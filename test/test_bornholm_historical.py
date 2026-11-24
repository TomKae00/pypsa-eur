# SPDX-License-Identifier: MIT
import json
from pathlib import Path

import pandas as pd
import pypsa
import pytest
import yaml

from scripts.run_bornholm_historical import dispatch_solver_options, run
from scripts.run_bornholm_matched import CENTRAL_FOLLOWUPS


def test_central_followup_designs():
    assert CENTRAL_FOLLOWUPS == (
        ("central_hybrid_matched", 1000, 2000, 3000),
        ("central_germany_alone", 0, 2000, 2000),
        ("central_denmark_alone", 1200, 0, 1200),
        ("central_smaller_island", 1200, 2000, 2000),
    )


@pytest.mark.parametrize("sector_coupled", [False, True])
def test_historical_pair_preserves_the_fleet_and_accounts_for_project_cost(
    tmp_path, sector_coupled
):
    n = pypsa.Network()
    hours = 3 if sector_coupled else 1
    n.set_snapshots(pd.date_range("2024", "2025", freq=f"{hours}h", inclusive="left"))
    n.snapshot_weightings.loc[:, :] = hours
    n.add("Carrier", ["AC", "DC", "offwind-ac", "gas"])
    n.add("Bus", ["DKE1", "DE00"], carrier="AC")
    n.add(
        "Generator",
        "fixed DK offshore",
        bus="DKE1",
        carrier="offwind-ac",
        p_nom=100,
        p_max_pu=pd.Series(0.5, index=n.snapshots),
    )
    n.add(
        "Generator",
        "fixed DE gas",
        bus="DE00",
        carrier="gas",
        p_nom=4000,
        marginal_cost=pd.Series(60.0, index=n.snapshots),
    )
    n.add(
        "Generator",
        "fixed DK gas",
        bus="DKE1",
        carrier="gas",
        p_nom=2000,
        marginal_cost=60,
    )
    n.add("Load", "DE demand", bus="DE00", p_set=2000)
    n.add("Load", "DK demand", bus="DKE1", p_set=1500)
    if sector_coupled:
        n.meta = yaml.safe_load(Path("config/config.default.yaml").read_text())
        n.meta["sector"]["enabled"] = True
        n.meta["planning_horizons"] = [2024]
        n.add("Carrier", ["urban central heat", "heat pump"])
        n.add("Bus", "DKE1 urban central heat", carrier="urban central heat")
        n.add(
            "Link",
            "DK heat pump",
            carrier="heat pump",
            bus0="DKE1",
            bus1="DKE1 urban central heat",
            p_nom=100,
            efficiency=3,
        )
        n.add("Load", "DK heat demand", bus="DKE1 urban central heat", p_set=150)
    source = tmp_path / "historical.nc"
    n.export_to_netcdf(source)

    template = pypsa.Network()
    template.add("Bus", ["BEI_HUB", "DKE1", "DE00"], carrier=["DC", "AC", "AC"])
    template.add(
        "Generator",
        "BEI_WIND",
        bus="BEI_HUB",
        carrier="offwind-ac",
        p_nom_extendable=True,
        p_nom_min=3000,
        p_nom_max=3000,
        capital_cost=1000,
    )
    for name, landing, rating in (
        ("BEI_TO_DK2", "DKE1", 1200),
        ("BEI_TO_DE", "DE00", 2000),
    ):
        template.add(
            "Link",
            name,
            bus0="BEI_HUB",
            bus1=landing,
            carrier="DC",
            p_nom_extendable=True,
            p_nom_min=rating,
            p_nom_max=rating,
            capital_cost=500,
            efficiency=0.95,
        )
        template.add(
            "Link",
            name + "-reversed",
            bus0=landing,
            bus1="BEI_HUB",
            carrier="DC",
            p_nom_extendable=True,
            p_nom_min=rating,
            p_nom_max=rating,
            efficiency=0.95,
        )
    template.links["reversed"] = template.links.index.str.endswith("-reversed")
    template_path = tmp_path / "template.nc"
    template.export_to_netcdf(template_path)

    output = tmp_path / "comparison"
    run(source, template_path, output, 2024, "highs", 1, sector_coupled=sector_coupled)
    result = json.loads((output / "summary.json").read_text())
    savings = 1500 * 8784 * 0.95 * 60
    assert result["snapshots"] == 8784 // hours
    assert result["annual_weight_hours"] == 8784
    assert result["with_project"]["wind_mwh"] == pytest.approx(1500 * 8784)
    assert result["project_annualised_capital_cost_eur"] == pytest.approx(4_600_000)
    assert result["operating_savings_source_eur"] == pytest.approx(savings)
    assert result["net_savings_after_project_cost_eur2025"] == pytest.approx(
        savings * 1.021 - 4_600_000
    )
    project = pypsa.Network(output / "with_project.nc")
    pd.testing.assert_series_equal(
        project.generators.loc[n.generators.index, "p_nom"], n.generators.p_nom
    )
    assert not project.generators.p_nom_extendable.any()
    assert not project.links.p_nom_extendable.any()
    if sector_coupled:
        assert project.links.at["DK heat pump", "p_nom"] == 100
        assert (project.links_t.p0["DK heat pump"] == 50).all()


def test_historical_sector_uses_the_original_2035_gurobi_preset(monkeypatch):
    from scripts import run_bornholm_diagnostic as diagnostic

    monkeypatch.setattr(
        diagnostic, "memory_snapshot", lambda: {"MemTotal": 128, "MemAvailable": 128}
    )
    options = dispatch_solver_options("gurobi", 16, sector_coupled=True)
    assert options["threads"] == 16
    assert options["method"] == 2
    assert options["crossover"] == 0
    assert options["NumericFocus"] == 1
    assert options["BarConvTol"] == 1e-8
    assert options["SoftMemLimit"] == 96
