# SPDX-License-Identifier: MIT
import numpy as np
import pandas as pd
import pypsa
import pytest

from scripts import hindcast_wind as wind
from scripts import historical_assets
from scripts.retrieve_danish_wind import municipality_area


def test_zone_specific_fit_preserves_weather_and_monthly_fleet(tmp_path):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2023-01-01", periods=48, freq="h"))
    raw = np.tile([0, 0.2, 0.5, 0.9], 12)
    rows = []
    for area, bus, factor in [("DK1", "DKW1", 0.5), ("DK2", "DKE1", 0.7)]:
        n.add("Bus", bus)
        n.add(
            "Generator",
            area + " onshore",
            bus=bus,
            carrier="onwind",
            p_nom=100,
            p_max_pu=raw,
        )
        n.add(
            "Generator",
            area + " offshore",
            bus=bus,
            carrier="offwind-ac",
            p_nom=100,
            p_max_pu=raw,
        )
        fleet = np.r_[np.full(24, 80.0), np.full(24, 90.0)]
        for t, g, c in zip(n.snapshots, raw * factor * fleet, fleet):
            rows.append(
                dict(
                    snapshot=t,
                    zone=area,
                    carrier="onwind",
                    generation_MW=g,
                    capacity_MW=c,
                    spot_price_EUR_MWh=50,
                )
            )
    p = tmp_path / "observed.csv.gz"
    pd.DataFrame(rows).to_csv(p, index=False)
    settings = dict(
        enabled=True,
        observations=str(p),
        minimum_fit_hours=10,
        factor_bounds=[0.3, 1.5],
    )
    report = wind.apply(n, settings, 2023)
    assert [r["factor"] for r in report["series"]] == pytest.approx([0.5, 0.7])
    np.testing.assert_allclose(n.generators_t.p_max_pu["DK1 offshore"], raw)
    assert n.generators.at["DK1 onshore", "p_nom"] == 90
    np.testing.assert_allclose(
        n.generators_t.p_max_pu["DK1 onshore"],
        raw * 0.5 * np.r_[np.full(24, 80 / 90), np.ones(24)],
    )
    assert wind.apply(n, settings, 2023) == report
    wind.validate(n, settings, 2023)
    with pytest.raises(ValueError, match="settings or observations changed"):
        wind.validate(n, {**settings, "minimum_fit_hours": 11}, 2023)


def test_nuclear_guard_preserves_total_capacity_and_other_countries():
    n = pypsa.Network()
    for b in ["SE03", "SE04", "FI00"]:
        n.add("Bus", b)
    n.add("Generator", "Ringhals", bus="SE04", carrier="nuclear", p_nom=2288)
    n.add("Generator", "Olkiluoto", bus="FI00", carrier="nuclear", p_nom=1600)
    before = n.generators.p_nom.sum()
    settings = {"enabled": True}
    historical_assets.apply(n, settings, 2023)
    assert n.generators.at["Ringhals", "bus"] == "SE03"
    assert n.generators.at["Olkiluoto", "bus"] == "FI00"
    assert n.generators.p_nom.sum() == before
    historical_assets.apply(n, settings, 2023)


def test_municipalities_cross_belt_explicitly():
    assert municipality_area(101) == "DK2"
    assert municipality_area(400) == "DK2"
    assert municipality_area(410) == "DK1"
    assert municipality_area(860) == "DK1"
    assert municipality_area(1) is None
