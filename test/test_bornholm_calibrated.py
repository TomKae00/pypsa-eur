# SPDX-License-Identifier: MIT
"""Physical solves and provenance checks for the combined Bornholm correction."""

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import pytest
import yaml

from scripts import bornholm_model_assumptions as assumptions
from scripts import run_bornholm_calibrated as campaign


@pytest.fixture
def model(tmp_path):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=3, freq="h"))
    n.snapshot_weightings.loc[:, :] = np.array([1, 2, 3])[:, None]
    for name, country in (("DK2", "DK"), ("SE4", "SE"), ("MK", "MK")):
        n.add("Bus", name, carrier="AC", country=country, location=name)
    n.add("Bus", "EU coal", carrier="coal")
    n.add(
        "Generator",
        "EU coal",
        bus="EU coal",
        carrier="coal",
        p_nom_extendable=True,
        marginal_cost=7,
    )
    n.add(
        "Store",
        "EU coal Store",
        bus="EU coal",
        carrier="coal",
        e_nom_extendable=True,
        e_cyclic=True,
        capital_cost=0.1,
    )
    for k, profile in enumerate(([0, 0.3, 0.6], [0, 0.5, 0.8], [0, 0.7, 1])):
        n.add(
            "Generator",
            f"DK2 {k} onwind",
            bus="DK2",
            carrier="onwind",
            p_nom=0,
            p_max_pu=profile,
        )
    n.add(
        "Generator",
        "MK onwind",
        bus="MK",
        carrier="onwind",
        p_nom=0,
        p_max_pu=[0.1, 0.2, 0.3],
    )
    n.add(
        "Generator",
        "DK2 solar",
        bus="DK2",
        carrier="solar",
        p_nom=0,
        p_max_pu=[0.0, 0.4, 0.0],
    )
    n.add(
        "Generator",
        "BEI_WIND",
        bus="DK2",
        carrier="offwind-ac",
        p_nom=0,
        p_max_pu=[0.3, 0.6, 0.9],
    )
    factors = tmp_path / "factors.csv"
    factors.write_text("country,carrier,factor\nDK,onwind,0.6\n")
    limits = tmp_path / "limits.csv"
    limits.write_text(
        "planning_year,zone0,zone1,forward_MW,reverse_MW\n2035,DK2,SE4,50,20\n"
    )
    settings = {
        "enabled": True,
        "onshore_calibration": str(factors),
        "transfer_limits": str(limits),
        "coal_store_MWh": 1000,
        "excluded_countries": ["MK"],
        "expansion_policy": "fixed",
    }
    return n, settings


def prepare(n, settings):
    assumptions.apply_resource_and_coal_assumptions(n, settings)
    assumptions.install_transfer_assumptions(n, settings, 2035)


def solve(n, settings):
    status, condition = n.optimize(
        solver_name="highs",
        threads=1,
        log_to_console=False,
        extra_functionality=lambda net, snaps: (
            assumptions.add_zonal_transfer_constraints(
                net,
                snaps,
                settings,
            )
        ),
    )
    assert (status, condition) == ("ok", "optimal")
    return assumptions.check_solved_assumptions(n, settings)


def test_calibration_preserves_chronology_potentials_and_other_resources(model):
    n, settings = model
    n.generators.loc["DK2 0 onwind", "p_nom"] = 10
    before = n.generators_t.p_max_pu.copy()
    capacities = n.generators[["p_nom", "p_nom_min", "p_nom_max"]].copy()
    result = assumptions.apply_resource_and_coal_assumptions(n, settings)
    names = [f"DK2 {k} onwind" for k in range(3)]
    pd.testing.assert_frame_equal(n.generators_t.p_max_pu[names], before[names] * 0.6)
    pd.testing.assert_frame_equal(
        n.generators_t.p_max_pu.drop(columns=names), before.drop(columns=names)
    )
    pd.testing.assert_frame_equal(n.generators[capacities.columns], capacities)
    assert result["resources"]["countries"]["DK"][
        "reference_after_cf"
    ] == pytest.approx(0.24)
    assert not n.stores.at["EU coal Store", "e_nom_extendable"]
    assert n.stores.at["EU coal Store", "e_nom"] == 1000
    assert n.generators.at["EU coal", "p_nom_extendable"]
    again = n.generators_t.p_max_pu.copy()
    assumptions.apply_resource_and_coal_assumptions(n, settings)
    pd.testing.assert_frame_equal(again, n.generators_t.p_max_pu)


def test_upward_scaling_clips_and_records_peaks(model):
    n, settings = model
    Path(settings["onshore_calibration"]).write_text(
        "country,carrier,factor\nDK,onwind,2\n"
    )
    audit = assumptions.apply_resource_and_coal_assumptions(n, settings)
    assert n.generators_t.p_max_pu["DK2 2 onwind"].tolist() == [0, 1, 1]
    assert audit["resources"]["assets"]["DK2 2 onwind"][
        "clipped_fraction"
    ] == pytest.approx(2 / 3)


@pytest.mark.parametrize(
    "rows,match",
    [
        ("DK,onwind,-1\n", "positive and finite"),
        ("DK,onwind,0.6\nDK,onwind,0.7\n", "Duplicate"),
        ("MK,onwind,2\n", "excluded country"),
    ],
)
def test_invalid_calibration_rejected(model, rows, match):
    n, settings = model
    Path(settings["onshore_calibration"]).write_text("country,carrier,factor\n" + rows)
    with pytest.raises(ValueError, match=match):
        assumptions.apply_resource_and_coal_assumptions(n, settings)


def test_coal_requires_independent_fuel_supply(model):
    n, settings = model
    n.remove("Generator", "EU coal")
    with pytest.raises(ValueError, match="independent fuel supply"):
        assumptions.apply_resource_and_coal_assumptions(n, settings)


@pytest.mark.parametrize("reverse", [False, True])
def test_parallel_lines_respect_combined_directional_limit(model, reverse):
    n, settings = model
    n.add("Line", "a", bus0="DK2", bus1="SE4", x=0.1, s_nom=100)
    n.add("Line", "b", bus0="SE4", bus1="DK2", x=0.1, s_nom=100)
    cheap, expensive = ("SE4", "DK2") if reverse else ("DK2", "SE4")
    n.add("Generator", "cheap", bus=cheap, p_nom=1000, marginal_cost=1)
    n.add("Generator", "expensive", bus=expensive, p_nom=1000, marginal_cost=100)
    n.add("Load", "demand", bus=expensive, p_set=100)
    prepare(n, settings)
    report = solve(n, settings)
    expected = 20 if reverse else 50
    assert np.allclose(n.generators_t.p.cheap, expected)
    assert report["corridors"][int(reverse)]["maximum_net_export_MW"] == pytest.approx(
        expected
    )


@pytest.mark.parametrize("reverse", [False, True])
def test_lossy_reverse_dc_and_ac_share_the_border_limit(model, reverse):
    n, settings = model
    n.add("Line", "ac", bus0="SE4", bus1="DK2", x=0.1, s_nom=10)
    n.add(
        "Link",
        "dc",
        bus0="DK2",
        bus1="SE4",
        carrier="DC",
        p_nom=100,
        p_min_pu=0,
        efficiency=0.9,
    )
    n.add(
        "Link",
        "dc-reversed",
        bus0="SE4",
        bus1="DK2",
        carrier="DC",
        p_nom=100,
        p_min_pu=0,
        efficiency=0.9,
    )
    cheap, expensive = ("SE4", "DK2") if reverse else ("DK2", "SE4")
    n.add("Generator", "cheap", bus=cheap, p_nom=1000, marginal_cost=1)
    n.add("Generator", "expensive", bus=expensive, p_nom=1000, marginal_cost=100)
    n.add("Load", "demand", bus=expensive, p_set=100)
    prepare(n, settings)
    report = solve(n, settings)
    expected = 20 if reverse else 50
    assert np.allclose(n.generators_t.p.cheap, expected)
    assert report["corridors"][int(reverse)]["maximum_net_export_MW"] == pytest.approx(
        expected
    )
    assert n.links_t.p0["dc-reversed" if reverse else "dc"].gt(0).all()


@pytest.mark.parametrize(
    "policy,expected_supply,expected_capacity",
    [
        ("proportional", 150, 300),
        ("fixed", 50, 100),
    ],
)
def test_reinforcement_uses_increment_above_existing_envelope(
    model, policy, expected_supply, expected_capacity
):
    n, settings = model
    settings["expansion_policy"] = policy
    n.add(
        "Link",
        "dc",
        bus0="DK2",
        bus1="SE4",
        carrier="DC",
        p_nom=100,
        p_nom_min=100,
        p_nom_max=400,
        p_nom_extendable=True,
        p_min_pu=-1,
        capital_cost=0.1,
    )
    n.add("Generator", "cheap", bus="DK2", p_nom=1000, marginal_cost=1)
    n.add("Generator", "expensive", bus="SE4", p_nom=1000, marginal_cost=100)
    n.add("Load", "demand", bus="SE4", p_set=150)
    prepare(n, settings)
    report = solve(n, settings)
    assert np.allclose(n.generators_t.p.cheap, expected_supply)
    assert n.links.at["dc", "p_nom_opt"] == pytest.approx(expected_capacity)
    assert report["corridors"][0]["limit_MW"] == pytest.approx(expected_supply)


def test_coal_buffer_does_not_limit_annual_fuel(model):
    n, settings = model
    n.add("Link", "dc", bus0="DK2", bus1="SE4", carrier="DC", p_nom=100, p_min_pu=-1)
    n.add("Load", "fuel demand", bus="EU coal", p_set=2000)
    prepare(n, settings)
    solve(n, settings)
    assert n.snapshot_weightings.generators @ n.generators_t.p[
        "EU coal"
    ] == pytest.approx(12000)
    assert n.stores.at["EU coal Store", "e_nom_opt"] == 1000


def test_metadata_roundtrip_and_stale_input_rejection(model, tmp_path):
    n, settings = model
    n.add("Link", "dc", bus0="DK2", bus1="SE4", carrier="DC", p_nom=100, p_min_pu=-1)
    prepare(n, settings)
    path = tmp_path / "prepared.nc"
    n.export_to_netcdf(path)
    restored = pypsa.Network(path)
    assert (
        assumptions.validate_assumptions(restored, settings)
        == n.meta[assumptions.AUDIT_KEY]
    )
    with pytest.raises(ValueError, match="requires its assumption configuration"):
        assumptions.validate_assumptions(restored, {})
    Path(settings["onshore_calibration"]).write_text(
        "country,carrier,factor\nDK,onwind,0.7\n"
    )
    with pytest.raises(ValueError, match="stale"):
        assumptions.validate_assumptions(restored, settings)


def test_three_class_workflow_uses_explicit_tables_and_isolated_targets():
    from scripts._helpers import mock_snakemake

    for run in ("bornholm_zonal_baseline", "bornholm_zonal_open_expandable"):
        context = campaign.context_for(run)
        assert context.config["renewable"]["onwind"]["resource_classes"] == 3
        assert campaign.PREFIX in context.input.network
        compose = mock_snakemake(
            "compose_network", configfiles=campaign.CONFIGS, run=run, horizon="2035"
        )
        assert compose.input.onshore_calibration.endswith("onshore_calibration.csv")
        assert compose.input.zonal_transfer_limits.endswith(
            "reviewed_transfer_limits.csv"
        )
        assert (
            compose.input.profile_onwind
            != mock_snakemake(
                "compose_network",
                configfiles=campaign.CONFIGS[:-1],
                run=run,
                horizon="2035",
            ).input.profile_onwind
        )
        assert (
            context.params.solving["solver_options"]["gurobi-calibrated-base"][
                "threads"
            ]
            == 16
        )
    settings = yaml.safe_load(Path(campaign.CONFIGS[-1]).read_text())[
        "bornholm_model_assumptions"
    ]
    assert json.loads(json.dumps(assumptions.fingerprint(settings)))[0]


@pytest.mark.parametrize("fuel", [30.0, 17.3312, float("nan")])
def test_biomass_preflight_checks_primary_fuel_not_transport_or_conversion(fuel):
    n = pypsa.Network()
    for name in ("DK solid biomass", "DE solid biomass"):
        n.add("Bus", name, carrier="solid biomass")
        n.add(
            "Generator", name, bus=name, carrier="solid biomass", marginal_cost=fuel
        )
    n.add(
        "Generator",
        "DK solid biomass transported",
        bus="DK solid biomass",
        carrier="solid biomass",
        marginal_cost=35,
    )
    if fuel == 30.0:
        report = campaign.check_biomass_fuel(n)
        assert report == {"fuel_EUR_MWh": 30.0, "currency_year": 2025, "suppliers": 2}
    else:
        with pytest.raises(ValueError, match="rebuild composed inputs"):
            campaign.check_biomass_fuel(n)


def test_biomass_preflight_rejects_missing_primary_supply():
    with pytest.raises(ValueError, match="Missing primary solid biomass"):
        campaign.check_biomass_fuel(pypsa.Network())


def test_reviewed_shared_constraints_replace_legacy_caps_once(model, tmp_path):
    from scripts import transmission_assumptions as transfers

    n, settings = model
    table = Path(settings["transfer_limits"])
    table.write_text(
        "planning_year,zone0,zone1,forward_MW,reverse_MW\n2035,DK2,SE4,15,10\n"
    )
    n.add("Link", "dc", bus0="DK2", bus1="SE4", carrier="DC", p_nom=100, p_min_pu=-1)
    n.add("Generator", "cheap", bus="DK2", p_nom=100, marginal_cost=1)
    n.add("Generator", "local", bus="SE4", p_nom=100, marginal_cost=100)
    n.add("Load", "demand", bus="SE4", p_set=40)
    common = {"enabled": True, "limits": str(table), "expansion_policy": "fixed"}
    assumptions.apply_resource_and_coal_assumptions(n, settings)
    transfers.install(n, common, 2035)  # Planning year, although weather is 2013.
    assumptions.install_transfer_assumptions(n, settings, 2035)
    path = tmp_path / "reviewed.nc"
    n.export_to_netcdf(path)
    n = pypsa.Network(path)
    assert n.meta[transfers.AUDIT_KEY]["year"] == 2035
    status = n.optimize(
        solver_name="highs",
        include_objective_constant=False,
        extra_functionality=lambda net, sns: assumptions.add_zonal_transfer_constraints(
            net, sns, settings
        ),
    )
    assert status == ("ok", "optimal")
    np.testing.assert_allclose(n.generators_t.p.cheap, 15, atol=1e-5)
    assert len(assumptions.check_solved_assumptions(n, settings)["corridors"]) == 2


@pytest.mark.parametrize("fail_comparison", [False, True])
def test_check_stage_validates_existing_inputs_without_rebuilding(
    tmp_path, fail_comparison
):
    scheduler = tmp_path / "bin"
    scheduler.mkdir()
    for name in ["sbatch", "sacct", "scancel"]:
        path = scheduler / name
        path.write_text("#!/usr/bin/env bash\nexit 99\n")
        path.chmod(0o755)
    calls = tmp_path / "calls.txt"
    wrapper = tmp_path / "pixi-sophia"
    wrapper.write_text(
        '#!/usr/bin/env bash\nprintf "%s\\n" "$*" >> "$CALLS"\n'
        'if [[ "$FAIL_COMPARISON" == 1 && "$*" == *scripts.check_bornholm_networks* ]]; then exit 7; fi\n'
    )
    wrapper.chmod(0o755)
    env = dict(
        os.environ,
        SLURM_JOB_ID="123",
        SLURM_JOB_PARTITION="workq",
        SLURM_SUBMIT_DIR=str(tmp_path),
        PIXI_SOPHIA=str(wrapper),
        CALLS=str(calls),
        FAIL_COMPARISON=str(int(fail_comparison)),
    )
    env["PATH"] = str(scheduler) + os.pathsep + env["PATH"]
    script = (
        Path(__file__).resolve().parents[1] / "scripts/bornholm_calibrated_workq.sbatch"
    )
    result = subprocess.run(
        ["bash", str(script), "check"], env=env, capture_output=True, text=True
    )
    commands = calls.read_text()
    assert "scripts.check_bornholm_networks" in commands
    assert "snakemake" not in commands
    assert "--list-build-targets" not in commands
    if fail_comparison:
        assert result.returncode == 7
        assert "--preflight-all" not in commands
    else:
        assert result.returncode == 0, result.stderr
        assert "--preflight-all" in commands
