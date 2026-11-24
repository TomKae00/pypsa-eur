# SPDX-License-Identifier: MIT
import json

import numpy as np
import pandas as pd
import pypsa
import pytest

from scripts.run_bornholm_historical import AF26_BORNHOLM_CF, run, scale_project_wind
from scripts.run_bornholm_matched import CENTRAL_FOLLOWUPS


def test_central_followup_designs():
    assert CENTRAL_FOLLOWUPS == (
        ("central_hybrid_matched", 1000, 2000, 3000),
        ("central_germany_alone", 0, 2000, 2000),
        ("central_denmark_alone", 1200, 0, 1200),
        ("central_smaller_island", 1200, 2000, 2000),
    )


def test_project_wind_target_after_nameplate_clipping():
    proxy = pd.Series([0.0, 0.1, 0.6, 0.9])
    profile, multiplier = scale_project_wind(proxy, AF26_BORNHOLM_CF)
    assert multiplier > 1
    assert profile.mean() == pytest.approx(AF26_BORNHOLM_CF)
    assert profile.iloc[0] == 0
    assert profile.iloc[-1] == 1
    assert profile.is_monotonic_increasing
    with pytest.raises(ValueError, match="CF target exceeds"):
        scale_project_wind(proxy, 0.8)


def test_historical_pair_preserves_the_fleet_and_accounts_for_project_cost(tmp_path):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2024", "2025", freq="h", inclusive="left"))
    n.add("Carrier", ["AC", "DC", "offwind-ac", "gas"])
    n.add("Bus", ["DKE1", "DE00"], carrier="AC")
    existing_cf = 0.36
    n.add(
        "Generator",
        "fixed DK offshore",
        bus="DKE1",
        carrier="offwind-ac",
        p_nom=100,
        p_max_pu=pd.Series(existing_cf, index=n.snapshots),
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
    run(source, template_path, output, 2024, "highs", 1)
    result = json.loads((output / "summary.json").read_text())
    expected_generation = 3000 * 8784 * AF26_BORNHOLM_CF
    savings = expected_generation * 0.95 * 60
    assert result["snapshots"] == 8784
    assert result["wind_profile_proxy_capacity_factor"] == pytest.approx(existing_cf)
    assert result["with_project"]["wind_mwh"] == pytest.approx(expected_generation)
    assert result["with_project"]["wind_available_mwh"] == pytest.approx(
        expected_generation
    )
    assert result["with_project"]["wind_realised_capacity_factor"] == pytest.approx(
        AF26_BORNHOLM_CF
    )
    assert result["project_annualised_capital_cost_eur"] == pytest.approx(4_600_000)
    assert result["operating_savings_source_eur"] == pytest.approx(savings)
    assert result["net_savings_after_project_cost_eur2025"] == pytest.approx(
        savings * 1.021 - 4_600_000
    )
    project = pypsa.Network(output / "with_project.nc")
    without = pypsa.Network(output / "without_project.nc")
    assert "bornholm_historical_wind" not in without.meta
    wind_audit = project.meta["bornholm_historical_wind"]
    assert wind_audit["wind_profile_capacity_factor"] == pytest.approx(AF26_BORNHOLM_CF)
    np.testing.assert_allclose(
        project.generators_t.p_max_pu["fixed DK offshore"], existing_cf
    )
    pd.testing.assert_series_equal(
        project.generators.loc[n.generators.index, "p_nom"], n.generators.p_nom
    )
    assert not project.generators.p_nom_extendable.any()
    assert not project.links.p_nom_extendable.any()
