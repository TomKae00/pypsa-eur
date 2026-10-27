# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for selected helpers in scripts/add_electricity.py."""

import numpy as np
import pandas as pd
import pypsa
import pytest
import xarray as xr

from scripts.add_electricity import (
    attach_hydro,
    attach_load,
    fill_hydro_max_hours,
    load_and_aggregate_powerplants,
)


def test_attach_load(tmp_path):
    """Clustered demand is attached per bus with scaling applied."""

    times = pd.date_range("2000-01-01", periods=2, freq="h")
    buses = ["zone_1", "zone_2"]
    values = np.array([[1.0, 2.0], [3.0, 4.0]])

    data = xr.DataArray(
        values,
        coords={"time": times, "bus": buses},
        dims=["time", "bus"],
        name="electricity demand (MW)",
    )
    load_path = tmp_path / "electricity_demand.nc"
    data.to_netcdf(load_path)

    n = pypsa.Network()
    n.set_snapshots(times)
    n.add("Bus", buses)

    attach_load(n, load_path.as_posix(), scaling=2.0)

    assert sorted(n.loads.index) == buses
    assert sorted(n.loads_t.p_set.columns) == buses
    np.testing.assert_allclose(n.loads_t.p_set[buses].values, 2.0 * values)


@pytest.fixture
def hydro_stats():
    return pd.DataFrame(
        {"E_store[TWh]": [0.2], "p_nom_discharge[GW]": [0.01]}, index=["NO"]
    )


@pytest.mark.parametrize("missing_duration", [np.nan, 0])
def test_hydro_completion_preserves_reported_energy_and_row_order(
    hydro_stats, missing_duration
):
    plants = pd.DataFrame(
        {
            "country": ["NO"] * 3,
            "p_nom": [10.0, 30.0, 60.0],
            "max_hours": [1000, missing_duration, 0],
        },
        index=["known", "missing_a", "missing_b"],
    )
    hours = fill_hydro_max_hours(
        plants, hydro_stats, "energy_capacity_totals_by_country"
    )
    assert hours["known"] == 1000
    assert hours["missing_a"] == pytest.approx(190000 / 90)
    assert hours["missing_b"] == hours["missing_a"]
    assert (plants.p_nom * hours).sum() == pytest.approx(200000)
    reversed_hours = fill_hydro_max_hours(
        plants.iloc[::-1], hydro_stats, "energy_capacity_totals_by_country"
    )
    pd.testing.assert_series_equal(hours, reversed_hours.reindex(hours.index))
    if pd.isna(missing_duration):
        assert pd.isna(plants.loc["missing_a", "max_hours"])
    else:
        assert plants.loc["missing_a", "max_hours"] == missing_duration


@pytest.mark.parametrize("durations", [[1000.0], [1000.0, 2000.0]])
def test_hydro_complete_records_are_not_rescaled(hydro_stats, durations, caplog):
    plants = pd.DataFrame({"country": "NO", "p_nom": 10.0, "max_hours": durations})
    result = fill_hydro_max_hours(
        plants, hydro_stats, "energy_capacity_totals_by_country"
    )
    pd.testing.assert_series_equal(result, plants.max_hours)
    assert "without scaling" in caplog.text


def test_hydro_reported_storage_exceeding_estimate_is_preserved(hydro_stats, caplog):
    plants = pd.DataFrame(
        {"country": ["NO", "NO"], "p_nom": [100.0, 10.0], "max_hours": [3000.0, np.nan]}
    )
    result = fill_hydro_max_hours(
        plants, hydro_stats, "energy_capacity_totals_by_country"
    )
    np.testing.assert_allclose(result, [3000, 0])
    assert "exceeds the country estimate" in caplog.text


@pytest.mark.parametrize("power", [[0.0], [0.0, 10.0]])
def test_hydro_zero_capacity_does_not_produce_infinite_duration(hydro_stats, power):
    plants = pd.DataFrame({"country": "NO", "p_nom": power, "max_hours": np.nan})
    result = fill_hydro_max_hours(
        plants, hydro_stats, "energy_capacity_totals_by_country"
    )
    assert np.isfinite(result).all()
    expected_energy = 200000 if sum(power) else 0
    assert (plants.p_nom * result).sum() == pytest.approx(expected_energy)


def test_hydro_large_installation_and_missing_country_fallback(hydro_stats):
    plants = pd.DataFrame(
        {
            "country": ["NO", "NO", "SE"],
            "p_nom": [10.0] * 3,
            "max_hours": [1000.0, np.nan, np.nan],
        }
    )
    result = fill_hydro_max_hours(
        plants, hydro_stats, "estimate_by_large_installations"
    )
    np.testing.assert_allclose(result, [1000, 20000, 6])
    result = fill_hydro_max_hours(
        plants, hydro_stats, "energy_capacity_totals_by_country"
    )
    assert result.iloc[2] == 6


@pytest.mark.parametrize("missing_duration", [np.nan, 0])
@pytest.mark.parametrize("include_missing_capacity", [False, True])
def test_reservoir_energy_survives_native_plant_aggregation(
    tmp_path, hydro_stats, missing_duration, include_missing_capacity, caplog
):
    """Mixed known/missing durations in the same bus must survive completion."""
    plants = pd.DataFrame(
        {
            "Name": ["known", "missing_a", "missing_b"],
            "Fueltype": "Hydro",
            "Technology": "Reservoir",
            "Country": "NO",
            "Capacity": [10.0, 30.0, 60.0],
            "Duration": [1000, missing_duration, 0],
            "Efficiency": np.nan,
            "DateIn": 2000,
            "DateOut": np.nan,
            "bus": ["NO1", "NO1", "NO2"],
        }
    )
    if include_missing_capacity:
        unknown = plants.iloc[[0]].copy()
        unknown["Name"] = "unknown_capacity"
        unknown["Capacity"] = np.nan
        unknown["Duration"] = np.nan
        plants = pd.concat([plants, unknown], ignore_index=True)
    plants_path = tmp_path / "plants.csv"
    stats_path = tmp_path / "hydro.csv"
    plants.to_csv(plants_path)
    source_csv = plants_path.read_bytes()
    hydro_stats.to_csv(stats_path)
    costs = pd.DataFrame(
        {
            "VOM": [0.0],
            "FOM": [0.0],
            "efficiency": [0.9],
            "capital_cost": [0.0],
            "marginal_cost": [0.0],
            "fuel": [0.0],
            "lifetime": [80.0],
        },
        index=["hydro"],
    )
    kwargs = dict(
        ppl_fn=plants_path,
        costs=costs,
        hydro_capacities=stats_path,
        hydro_max_hours="energy_capacity_totals_by_country",
    )
    aggregated = load_and_aggregate_powerplants(**kwargs)
    individual = load_and_aggregate_powerplants(**kwargs, exclude_carriers=["hydro"])
    assert plants_path.read_bytes() == source_csv
    assert len(individual) == 3
    assert np.isfinite(individual.p_nom).all()
    assert not individual.index.str.contains("unknown_capacity").any()
    if include_missing_capacity:
        assert (
            "Excluding 1 reservoir records with missing discharge capacity"
            in caplog.text
        )
        assert "unknown_capacity" in caplog.text
        assert "omitted MW and storage are unknown" in caplog.text
    energy = aggregated.p_nom * aggregated.max_hours
    individual_energy = (
        (individual.p_nom * individual.max_hours).groupby(individual.bus).sum()
    )
    assert energy.sum() == pytest.approx(200000)
    assert energy["NO1 hydro"] == pytest.approx(10000 + 190000 / 3)
    np.testing.assert_allclose(
        energy.groupby(aggregated.bus).sum().sort_index(),
        individual_energy.sort_index(),
    )


@pytest.mark.parametrize("capacity", [np.nan, -1.0, np.inf, -np.inf])
def test_hydro_invalid_capacity_reports_source_record(hydro_stats, capacity):
    plants = pd.DataFrame(
        {
            "name": ["bad_reservoir"],
            "country": ["NO"],
            "bus": ["NO1"],
            "p_nom": [capacity],
            "max_hours": [np.nan],
        },
        index=["source_row_17"],
    )
    with pytest.raises(ValueError, match="finite and non-negative") as exc:
        fill_hydro_max_hours(plants, hydro_stats, "energy_capacity_totals_by_country")
    for detail in ["source_row_17", "bad_reservoir", "NO1"]:
        assert detail in str(exc.value)


def test_attach_hydro_keeps_completed_zero_duration(tmp_path):
    """A zero residual must not turn into fallback storage during attachment."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2023-01-01", periods=2, freq="h"))
    n.add("Bus", "NO1", country="NO")
    plants = pd.DataFrame(
        {
            "bus": "NO1",
            "country": "NO",
            "carrier": "hydro",
            "p_nom": [10.0, 10.0],
            "max_hours": [100.0, 0.0],
        },
        index=["reported", "zero_residual"],
    )
    costs = pd.DataFrame(
        {
            "CO2 intensity": [0.0],
            "efficiency": [0.9],
            "capital_cost": [0.0],
            "marginal_cost": [0.0],
        },
        index=["hydro"],
    )
    inflow_path = tmp_path / "inflow.nc"
    xr.DataArray(
        [[20.0], [40.0]],
        dims=("time", "countries"),
        coords={"time": n.snapshots.rename("time"), "countries": ["NO"]},
    ).to_netcdf(inflow_path)
    attach_hydro(n, costs, plants, inflow_path, ["hydro"])
    np.testing.assert_allclose(n.storage_units.max_hours, [100, 0])
    np.testing.assert_allclose(n.storage_units_t.inflow.sum(axis=1), [20, 40])
