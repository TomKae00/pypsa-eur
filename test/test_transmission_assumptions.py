# SPDX-License-Identifier: MIT
"""Numerical dispatch tests for the operating constraints, including NetCDF reuse."""

import hashlib

import numpy as np
import pandas as pd
import pypsa
import pytest

from scripts import transmission_assumptions as t


def pair(tmp_path, *, year=2023, forward=30, reverse=10, extendable=False):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range(f"{year}-01-01", periods=3, freq="h"))
    n.add("Bus", "DKE1", carrier="AC")
    n.add("Bus", "SE04", carrier="AC")
    n.add(
        "Line",
        "AC",
        bus0="DKE1",
        bus1="SE04",
        x=0.1,
        s_nom=10 if extendable else 100,
        s_nom_min=10 if extendable else 0,
        s_nom_extendable=extendable,
        capital_cost=1,
    )
    n.add("Generator", "cheap", bus="DKE1", p_nom=100, marginal_cost=1)
    n.add("Generator", "local", bus="SE04", p_nom=100, marginal_cost=50)
    n.add("Load", "demand", bus="SE04", p_set=40)
    p = tmp_path / "limits.csv"
    pd.DataFrame(
        [
            dict(
                planning_year=year,
                zone0="DK2",
                zone1="SE4",
                forward_MW=forward,
                reverse_MW=reverse,
            )
        ]
    ).to_csv(p, index=False)
    settings = {"enabled": True, "limits": str(p), "expansion_policy": "fixed"}
    return n, settings


def solve(n, settings):
    t.install(n, settings, n.snapshots[0].year)
    status = n.optimize(
        solver_name="highs",
        include_objective_constant=False,
        extra_functionality=lambda net, sns: t.add_constraints(net, sns, settings),
    )
    assert status == ("ok", "optimal")
    t.check_solution(n, settings)


def test_hourly_export_roundtrip_and_dispatch(tmp_path):
    n, settings = pair(tmp_path)
    h = tmp_path / "hourly.csv"
    pd.DataFrame(
        dict(
            snapshot=n.snapshots,
            zone0="DK2",
            zone1="SE4",
            forward_MW=[10, 20, 0],
            reverse_MW=[5, 5, 5],
        )
    ).to_csv(h, index=False)
    settings["hourly_limits"] = str(h)
    t.install(n, settings, 2023)
    p = tmp_path / "prepared.nc"
    n.export_to_netcdf(p)
    m = pypsa.Network(p)
    solve(m, settings)
    np.testing.assert_allclose(m.lines_t.p0.AC, [10, 20, 0], atol=1e-5)
    np.testing.assert_allclose(m.generators_t.p.local, [30, 20, 40], atol=1e-5)
    m.set_snapshots(m.snapshots[:2])
    with pytest.raises(ValueError, match="snapshots changed"):
        t.validate(m, settings)


def test_lossy_reverse_pair_counted_at_sending_boundary(tmp_path):
    n, settings = pair(tmp_path, forward=50, reverse=17)
    n.lines.loc["AC", "s_nom"] = 30
    n.add(
        "Link",
        "cable",
        bus0="DKE1",
        bus1="SE04",
        carrier="DC",
        p_nom=30,
        efficiency=0.9,
    )
    n.add(
        "Link",
        "cable-reversed",
        bus0="SE04",
        bus1="DKE1",
        carrier="DC",
        p_nom=30,
        efficiency=0.9,
    )
    n.loads.loc["demand", "p_set"] = 80
    solve(n, settings)
    export = n.lines_t.p0.AC + n.links_t.p0.cable - 0.9 * n.links_t.p0["cable-reversed"]
    np.testing.assert_allclose(export, 50, atol=1e-5)
    assert (n.generators_t.p.local > 30).all()
    # Reverse the economic incentive; independently check the other direction.
    n.generators.loc["cheap", "marginal_cost"] = 100
    n.generators.loc["local", "marginal_cost"] = 1
    n.loads.loc["demand", "bus"] = "DKE1"
    solve(n, settings)
    receiving_export = (
        -n.lines_t.p0.AC + n.links_t.p0["cable-reversed"] - 0.9 * n.links_t.p0.cable
    )
    np.testing.assert_allclose(receiving_export, 17, atol=1e-5)


def test_investment_scales_operating_ceiling(tmp_path):
    n, settings = pair(tmp_path, year=2035, forward=5, reverse=5, extendable=True)
    settings["expansion_policy"] = "proportional"
    solve(n, settings)
    assert n.lines.at["AC", "s_nom_opt"] == pytest.approx(80)
    np.testing.assert_allclose(n.lines_t.p0.AC, 40, atol=1e-5)


def test_polish_gross_profile_cannot_cancel_export_against_import(tmp_path):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035-01-01", periods=1, freq="h"))
    for bus in ["PL00", "DE00", "CZ00"]:
        n.add("Bus", bus, carrier="AC")
    n.add("Line", "PLDE", bus0="PL00", bus1="DE00", x=0.1, s_nom=20)
    n.add("Line", "PLCZ", bus0="PL00", bus1="CZ00", x=0.1, s_nom=20)
    for name, bus, cost in [
        ("PLgen", "PL00", 2),
        ("CZgen", "CZ00", 1),
        ("DEgen", "DE00", 100),
    ]:
        n.add("Generator", name, bus=bus, p_nom=20, marginal_cost=cost)
    n.add("Load", "DEload", bus="DE00", p_set=8)
    p, j = tmp_path / "limits.csv", tmp_path / "joint.csv"
    pd.DataFrame(
        [
            dict(
                planning_year=2035, zone0="PL00", zone1=z, forward_MW=20, reverse_MW=20
            )
            for z in ["DE-LU", "CZ00"]
        ]
    ).to_csv(p, index=False)
    pd.DataFrame(
        [
            dict(
                planning_year=2035,
                constraint="profile",
                zone0="PL00",
                zone1=z,
                coefficient=1,
                positive_only=True,
                limit_MW=6,
            )
            for z in ["DE-LU", "CZ00"]
        ]
    ).to_csv(j, index=False)
    settings = {"enabled": True, "limits": str(p), "joint_limits": str(j)}
    solve(n, settings)
    assert n.lines_t.p0.PLDE.iloc[0] == pytest.approx(6)
    assert n.generators_t.p.DEgen.iloc[0] == pytest.approx(2)


def test_asset_date_derating_idempotency_and_stale_source(tmp_path):
    n, settings = pair(tmp_path)
    n.set_snapshots(pd.to_datetime(["2023-12-28", "2023-12-30"]))
    n.add(
        "Link",
        "public_asset",
        bus0="DKE1",
        bus1="SE04",
        carrier="DC",
        p_nom=1000,
        p_min_pu=-1,
    )
    p = tmp_path / "assets.csv"
    pd.DataFrame(
        [
            dict(
                year_min=2023,
                year_max=2023,
                asset="public_asset",
                nominal_MW=800,
                operating_MW=400,
                available_from_utc="2023-12-30T00:00:00Z",
            )
        ]
    ).to_csv(p, index=False)
    settings["asset_rules"] = str(p)
    t.apply_asset_rules(n, settings, 2023)
    t.apply_asset_rules(n, settings, 2023)
    assert n.links.at["public_asset", "p_nom"] == 800
    np.testing.assert_allclose(n.links_t.p_max_pu.public_asset, [0, 0.5])
    np.testing.assert_allclose(n.links_t.p_min_pu.public_asset, [0, -0.5])
    p.write_text(p.read_text() + "\n")
    with pytest.raises(ValueError, match="Asset rules changed"):
        t.apply_asset_rules(n, settings, 2023)


def test_reject_maxbex_as_simultaneous_limits(tmp_path):
    n, settings = pair(tmp_path)
    df = pd.read_csv(settings["limits"])
    df["basis"] = "JAO MaxBEX median"
    df.to_csv(settings["limits"], index=False)
    with pytest.raises(ValueError, match="MaxBEX is a diagnostic"):
        t.install(n, settings, 2023)


def test_native_loss_split_preserves_both_directional_hourly_limits():
    from scripts.prepare_sector_network import lossy_bidirectional_links

    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2023-01-01", periods=3, freq="h"))
    for bus in ["a", "b"]:
        n.add("Bus", bus)
    n.add(
        "Link",
        "dated cable",
        bus0="a",
        bus1="b",
        carrier="DC",
        p_nom=100,
        p_min_pu=-1,
        p_max_pu=1,
    )
    n.links_t.p_min_pu["dated cable"] = [0, -0.4, -0.4]
    n.links_t.p_max_pu["dated cable"] = [0, 0.8, 0.8]
    lossy_bidirectional_links(n, "DC", {"efficiency_static": 0.9})
    assert n.links.p_min_pu.eq(0).all()
    assert n.links_t.p_min_pu.empty
    np.testing.assert_allclose(
        n.links_t.p_max_pu["dated cable-reversed"], [0, 0.4, 0.4]
    )
    n.add("Generator", "a gen", bus="a", p_nom=200, marginal_cost=[1, 1, 100])
    n.add("Generator", "b gen", bus="b", p_nom=200, marginal_cost=[50, 50, 1])
    n.add("Load", "a load", bus="a", p_set=[0, 0, 50])
    n.add("Load", "b load", bus="b", p_set=[50, 100, 0])
    assert n.optimize(solver_name="highs", include_objective_constant=False) == (
        "ok",
        "optimal",
    )
    np.testing.assert_allclose(n.links_t.p0["dated cable"], [0, 80, 0], atol=1e-5)
    np.testing.assert_allclose(
        n.links_t.p0["dated cable-reversed"], [0, 0, 40], atol=1e-5
    )


@pytest.mark.parametrize("unit", ["ns", "us", "ms", "s"])
def test_snapshot_hash_is_independent_of_datetime_storage_unit(tmp_path, unit):
    n, settings = pair(tmp_path)
    n.set_snapshots(n.snapshots.as_unit(unit))
    audit = t.install(n, settings, 2023)
    assert audit["snapshot_hash_unit"] == "ns"
    original = audit["snapshots_sha256"]
    n.set_snapshots(n.snapshots.as_unit("ns"))
    assert t.validate(n, settings)["snapshots_sha256"] == original


@pytest.mark.parametrize("unit", ["ns", "us", "ms", "s"])
def test_legacy_three_hour_hash_survives_netcdf_roundtrip(tmp_path, unit):
    n, settings = pair(tmp_path, year=2035)
    n.set_snapshots(pd.date_range("2013-01-01", periods=3, freq="3h").as_unit(unit))
    audit = t.install(n, settings, 2035)
    audit.pop("snapshot_hash_unit")
    audit["snapshots_sha256"] = hashlib.sha256(
        t._utc(n.snapshots).asi8.tobytes()
    ).hexdigest()
    output = tmp_path / "legacy_3h.nc"
    n.export_to_netcdf(output)
    m = pypsa.Network(output)
    t.validate(m, settings)
    # A fractional-time change must not be hidden by rounding to a legacy unit.
    m.set_snapshots(m.snapshots.as_unit("ns") + pd.Timedelta(nanoseconds=1))
    with pytest.raises(ValueError, match="snapshots changed"):
        t.validate(m, settings)


@pytest.mark.parametrize("change", ["shift", "drop", "reverse"])
def test_legacy_snapshot_compatibility_rejects_real_time_changes(tmp_path, change):
    n, settings = pair(tmp_path)
    audit = t.install(n, settings, 2023)
    audit.pop("snapshot_hash_unit")
    audit["snapshots_sha256"] = t._snapshot_hash(n.snapshots, "us")
    changed = {
        "shift": n.snapshots + pd.Timedelta(hours=1),
        "drop": n.snapshots[:-1],
        "reverse": n.snapshots[::-1],
    }[change]
    n.set_snapshots(changed)
    with pytest.raises(ValueError, match="snapshots changed"):
        t.validate(n, settings)
