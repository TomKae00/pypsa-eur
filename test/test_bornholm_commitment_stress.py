# SPDX-License-Identifier: MIT
"""Focused checks for the commitment-stress interventions."""

import pandas as pd
import pypsa
import pytest

from scripts import bornholm_2035 as campaign
from scripts import run_bornholm_commitment_stress as stress


def test_same_demand_dose_in_either_country():
    cfg = campaign.read_config()
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035", periods=4, freq="3h"))
    n.snapshot_weightings.loc[:, :] = 2190.0
    for bus, country in (("DKW1", "DK"), ("DKE1", "DK"), ("DE00", "DE")):
        n.add("Bus", bus, country=country, carrier="AC")

    matrix = stress.cases()
    assert [row["index"] for row in matrix] == list(range(9))
    for name, bus in (("demand_dk2", "DKE1"), ("demand_de", "DE00")):
        selected = [row for row in matrix if row["stress"] == name]
        assert {row["layout"] for row in selected} == {
            "split_radial",
            "treaty_hybrid",
            "matched_hybrid",
        }
        stressed = stress.stressed_config(cfg, name)
        settings = campaign.definitions(stressed)["central_hybrid"]["bornholm_2035"]
        clone = n.copy()
        campaign.apply(clone, settings)
        assert len(clone.loads) == 1
        assert clone.loads.bus.iat[0] == bus
        assert clone.loads.p_set.iat[0] * 8760 / 1e6 == pytest.approx(
            stress.ADDED_DEMAND_TWH
        )
    assert not any(cfg["worlds"]["central"]["additional_data_centres_TWh"].values())


def test_offshore_floor_is_above_reference_and_uses_native_headroom(monkeypatch):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035", periods=4, freq="3h"))
    n.snapshot_weightings.loc[:, :] = 2190.0
    for bus, country in (("DKE1", "DK"), ("DE00", "DE")):
        n.add("Bus", bus, country=country, carrier="AC")
    project = [0.15, 0.35, 0.75, 0.5]
    n.add("Generator", "BEI_WIND_DK", bus="DKE1", p_max_pu=project)
    for bus in ("DKE1", "DE00"):
        for label, profile in (("aligned", project), ("other", [0.6, 0.5, 0.2, 0.4])):
            n.add(
                "Generator",
                f"{bus}-{label}",
                bus=bus,
                carrier="offwind-ac",
                p_nom_extendable=True,
                p_nom_min=10,
                p_nom_max=100,
                p_max_pu=profile,
            )
    solved = n.copy()
    solved.generators.loc[["DKE1-aligned", "DE00-aligned"], "p_nom_opt"] = 30
    solved.generators.loc[["DKE1-other", "DE00-other"], "p_nom_opt"] = 20

    with monkeypatch.context() as patch:
        patch.setattr(stress, "ADDED_OFFSHORE_MW", 85)
        plan = stress.offshore_minimums(n, solved)
    assert plan["minimums_mw"] == {
        "DKE1-aligned": 100.0,
        "DKE1-other": 35.0,
        "DE00-aligned": 100.0,
        "DE00-other": 35.0,
    }
    assert plan["added_sites"]["DKE1-aligned"]["extra_mw"] == 70.0
    assert plan["added_sites"]["DE00-other"]["extra_mw"] == 15.0
    assert n.generators.at["DE00-aligned", "p_nom_min"] == 10
    with pytest.raises(ValueError, match="Insufficient profiled offshore headroom"):
        stress.offshore_minimums(n, solved)
