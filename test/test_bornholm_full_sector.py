"""The extra run must accept retrofit gas bounds and reject simplified inputs."""

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import pypsa
import pytest

from scripts import run_bornholm_full_sector as full


@pytest.fixture
def network():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=2920, freq="3h"))
    n.snapshot_weightings.loc[:, :] = 3
    n.add("Bus", [f"AC{i}" for i in range(52)], carrier="AC")
    n.add("Bus", ["gas0", "gas1"], carrier="gas")
    n.add("Bus", ["H0", "H1"], carrier="H2")
    n.add("Bus", ["C0", "C1"], carrier="co2 stored")
    for name, carrier, bus0, bus1 in (
        ("new H2", "H2 pipeline", "H0", "H1"),
        ("converted H2", "H2 pipeline retrofitted", "H0", "H1"),
        ("CO2", "CO2 pipeline", "C0", "C1"),
        ("gas", "gas pipeline", "gas0", "gas1"),
    ):
        n.add(
            "Link",
            name,
            carrier=carrier,
            bus0=bus0,
            bus1=bus1,
            p_nom_extendable=True,
            p_nom=100,
            p_nom_min=0,
            p_nom_max=100,
        )
    for carrier in ("shipping oil", "shipping methanol", "kerosene for aviation"):
        n.add("Load", carrier, bus="AC0", carrier=carrier, p_set=10)
    n.add("GlobalConstraint", "CO2Limit", type="co2_atmosphere", constant=10000)
    n.meta = {
        "run": {"prefix": full.PREFIX},
        "sector": {
            k: True
            for k in (*full.RESTORED, "co2_spatial", "co2_network", "gas_network")
        },
        "electricity": {"transmission_limit": "vfixed"},
        "renewable": {
            k: {"resource_classes": v} for k, v in full.calibrated.CLASSES.items()
        },
    }
    return n


def test_retrofit_scope_and_actual_carbon_budget_are_recorded(network):
    report = full.check_scope(network)
    json.dumps(report, allow_nan=False)
    assert report["annual_weight_hours"] == 8760
    assert report["carbon_budget"]["CO2Limit"]["constant"] == 10000
    assert report["final_demand_TWh_by_carrier"]["shipping methanol"] == pytest.approx(
        0.0876
    )
    assert report["links_by_carrier"]["H2 pipeline retrofitted"] == 1


@pytest.mark.parametrize(
    "change",
    [
        "H2_network",
        "shipping",
        "CO2 pipeline",
        "aviation demand",
        "fixed gas",
        "weights",
        "carbon budget",
    ],
)
def test_stale_or_inconsistent_composed_inputs_are_rejected(network, change):
    if change in full.RESTORED:
        network.meta["sector"][change] = False
    elif change == "CO2 pipeline":
        network.remove("Link", "CO2")
    elif change == "aviation demand":
        network.loads.loc["kerosene for aviation", "p_set"] = 0
    elif change == "fixed gas":
        network.links.loc["gas", "p_nom_extendable"] = False
    elif change == "weights":
        network.snapshot_weightings.loc[:, "stores"] = 1
    else:
        network.remove("GlobalConstraint", "CO2Limit")
    with pytest.raises(ValueError):
        full.check_scope(network)


def test_effective_native_config_is_one_isolated_central_hybrid():
    context = full.study.context_for(full.RUN, configfiles=full.CONFIGS)
    reference = full.study.context_for(
        full.RUN, configfiles=full.final_campaign.CONFIGS
    )
    assert Path(context.input.network).resolve() == Path(full.TARGET).resolve()
    assert context.config["run"]["name"] == [full.RUN]
    assert context.config["run"]["shared_resources"]["policy"] is False
    assert context.config["clustering"]["temporal"]["averaging"] == "3h"
    for key in full.RESTORED:
        assert context.config["sector"][key] is True
        assert reference.config["sector"][key] is False
    assert context.config["sector"]["co2_network"] is True
    assert reference.config["sector"]["co2_network"] is True
    for key in (
        "bornholm_2035",
        "energy_islands",
        "renewable",
        "industry",
        "co2_budget",
    ):
        assert context.config[key] == reference.config[key]
    assert context.config["sector"]["shipping_methanol_share"][2035] == 0.35
    assert context.config["sector"]["shipping_oil_share"][2035] == 0.65
    assert context.config["sector"]["imports"]["enable"] is False
    assert context.config["data"]["swiss_energy_balances"] == {
        "source": "primary",
        "version": "2025-10-01",
    }


def test_shared_annual_runner_uses_new_target_and_preflight(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    output = Path("results") / full.PREFIX / "job" / full.RUN
    context = SimpleNamespace(input=SimpleNamespace(network=full.TARGET))

    # Stop deliberately at the custom preflight; the legacy fixed-gas check
    # must never run for this expanded sector configuration.
    def preflight(context):
        raise ValueError("expanded preflight reached")

    run = full.annual.AnnualRun(
        full.RUN, full.PREFIX, tuple(full.CONFIGS), full.TARGET, preflight, lambda n: {}
    )
    with (
        patch.object(full.study, "context_for", return_value=context),
        patch.object(
            full.study, "preflight", side_effect=AssertionError("legacy check")
        ),
        patch.object(full.annual, "configure", return_value={}),
        patch("scripts._helpers.configure_logging"),
        patch.object(full.annual.accuracy, "sha256", return_value="input-hash"),
        patch.object(full.annual, "code_provenance", return_value={}),
        pytest.raises(ValueError, match="expanded preflight reached"),
    ):
        full.annual.solve(output, annual_run=run)
    report = json.loads((output / "summary.json").read_text())
    assert report["scenario"] == full.RUN
    assert report["failed_phase"] == "preflight"
    assert report["source"] == full.TARGET
