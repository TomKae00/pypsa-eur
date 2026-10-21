# SPDX-License-Identifier: MIT
"""Prepare hourly fuel/CO2 inputs from a retrieved market-data collection.

All fuel outputs use EUR/MWh of net thermal energy; CO2 uses EUR/tCO2.
This step does not modify generation capacities, build a network or solve one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from hindcast import digest, write_json


def positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def read_sources(directory, config, years):
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if not set(years).issubset(manifest["requested_years"]):
        raise ValueError("Selected years are not covered by the source collection")
    records = {Path(item["path"]).name: item for item in manifest["outputs"]}
    result, provenance = {}, {}
    for name in ["gas", "coal", "co2", "oil", "fx"]:
        spec = config[name]
        path = directory / spec["file"]
        if path.name not in records or digest(path) != records[path.name]["sha256"]:
            raise ValueError(f"Source file does not match its collection manifest: {path}")
        frame = pd.read_csv(path)
        dates = pd.DatetimeIndex(pd.to_datetime(frame["date"], errors="raise", utc=True))
        if dates.hasnans or dates.has_duplicates or not dates.is_monotonic_increasing:
            raise ValueError(f"Source dates must be unique and increasing: {path}")
        if not dates.equals(dates.normalize()):
            raise ValueError(f"Expected daily observation labels: {path}")
        values = pd.to_numeric(frame[spec["column"]], errors="raise").to_numpy(dtype=float)
        if not np.isfinite(values[~np.isnan(values)]).all():
            raise ValueError(f"Non-finite source values: {path}")
        series = pd.Series(values, index=dates, name=name).dropna()
        if series.empty:
            raise ValueError(f"No usable observations: {path}")
        lag = config["time"]["effective_lag_days"][name]
        if isinstance(lag, bool) or not isinstance(lag, int) or lag < 0:
            raise ValueError(f"Effective lag for {name} must be a non-negative integer")
        result[name] = pd.DataFrame({
            "value": series.to_numpy(), "observation_date": series.index,
        }, index=series.index + pd.Timedelta(days=lag))
        result[name].index.name = "effective_date"
        provenance[name] = {
            "path": str(path), "sha256": digest(path), "column": spec["column"],
            "usable_observations": len(series), "missing_source_rows": int(np.isnan(values).sum()),
        }
    if (result["fx"]["value"] <= 0).any():
        raise ValueError("USD/EUR rates must be positive")
    provenance["collection_manifest"] = {
        "path": str(manifest_path), "sha256": digest(manifest_path),
    }
    return result, provenance


def previous(frame, target, max_age, name):
    """Look up only past/current observations; never interpolate or backfill."""
    positions = frame.index.get_indexer(target, method="pad")
    if (positions < 0).any():
        first = target[np.flatnonzero(positions < 0)[0]]
        raise ValueError(f"{name} has no seed observation on or before {first}")
    selected = frame.iloc[positions].copy()
    effective = frame.index.take(positions)
    ages = (target - effective) / pd.Timedelta(days=1)
    # Age limits apply to observation days, so a value held through day 10
    # expires at the start of day 11 rather than part-way through day 10.
    stale = np.floor(ages) > max_age
    if stale.any():
        first = np.flatnonzero(stale)[0]
        raise ValueError(
            f"{name} observation is too old at {target[first]} "
            f"(effective date {effective[first]}, limit {max_age} days)"
        )
    selected.index = target
    return selected, effective, np.asarray(ages)


def convert(sources, config):
    if config["time"]["holding"] != "previous_observation":
        raise ValueError("Only previous_observation holding is implemented")
    daily_age = positive(config["time"]["daily_max_age_days"], "Daily maximum age")
    weekly_age = positive(config["time"]["weekly_max_age_days"], "Weekly maximum age")
    gas = config["gas"]
    if gas["model_basis"] != "LHV" or gas["quoted_basis_assumption"] not in {"LHV", "HHV"}:
        raise ValueError("Gas basis must be HHV or LHV, with LHV output")
    ratio = positive(gas["lhv_over_hhv"], "Gas LHV/HHV ratio")
    if ratio > 1:
        raise ValueError("LHV/HHV cannot exceed one")
    gas_factor = 1 / ratio if gas["quoted_basis_assumption"] == "HHV" else 1.0
    coal_mwh_t = (
        positive(config["coal"]["net_calorific_value_kcal_kg"], "Coal NCV")
        * positive(config["coal"]["kcal_to_kj"], "kcal conversion") / 3600
    )
    oil_mwh_t = positive(config["oil"]["net_calorific_value_gj_t"], "Oil NCV") / 3.6
    lignite_factor = positive(config["lignite"]["coal_price_factor"], "Lignite factor")

    # Convert on the coal quote's effective date, then hold the converted quote.
    # A day with only a new FX observation does not create a new coal quote.
    fx, fx_dates, fx_ages = previous(sources["fx"], sources["coal"].index, daily_age, "Coal FX")
    converted = {name: frame.copy() for name, frame in sources.items() if name != "fx"}
    converted["gas"]["value"] *= gas_factor
    converted["coal"]["value"] /= fx["value"].to_numpy() * coal_mwh_t
    converted["coal"]["fx_effective_date"] = fx_dates
    converted["coal"]["fx_observation_date"] = fx["observation_date"].to_numpy()
    converted["oil"]["value"] /= oil_mwh_t
    factors = {
        "gas_price_multiplier": gas_factor, "gas_lhv_over_hhv_assumption": ratio,
        "coal_mwh_lhv_per_tonne": coal_mwh_t, "oil_mwh_lhv_per_tonne": oil_mwh_t,
        "lignite_coal_energy_price_ratio": lignite_factor,
        "maximum_fx_age_at_coal_quote_days": float(fx_ages.max()),
    }
    return converted, factors, daily_age, weekly_age


def annual(year, converted, config, factors, daily_age, weekly_age):
    index = pd.date_range(f"{year}-01-01", f"{year + 1}-01-01", freq="h", inclusive="left", tz="UTC")
    if config["drop_leap_day"]:
        index = index[~((index.month == 2) & (index.day == 29))]
    fuel = pd.DataFrame(index=index)
    dates = pd.DataFrame(index=index)
    summary = {}
    for name in ["gas", "coal", "oil", "co2"]:
        frame, effective, ages = previous(
            converted[name], index, weekly_age if name == "oil" else daily_age, name
        )
        if name == "co2":
            co2 = frame[["value"]].rename(columns={"value": "co2"})
        else:
            fuel[name] = frame["value"]
        dates[name + "_observation_date"] = frame["observation_date"]
        dates[name + "_effective_date"] = effective
        if name == "coal":
            dates["coal_fx_observation_date"] = frame["fx_observation_date"]
            dates["coal_fx_effective_date"] = frame["fx_effective_date"]
        summary[name] = {
            "first_observation_used": str(frame["observation_date"].iloc[0]),
            "last_observation_used": str(frame["observation_date"].iloc[-1]),
            "maximum_age_days": float(ages.max()),
            "distinct_observations_used": int(frame["observation_date"].nunique()),
        }
    fuel["lignite"] = fuel["coal"] * factors["lignite_coal_energy_price_ratio"]
    fuel = fuel[["gas", "coal", "lignite", "oil"]]
    for frame in [fuel, co2]:
        if not np.isfinite(frame.to_numpy()).all():
            raise ValueError(f"Missing or non-finite converted inputs in {year}")
    for frame in [fuel, co2, dates]:
        frame.index.name = "utc_timestamp"
    stats = pd.concat([fuel, co2], axis=1).agg(["min", "mean", "max"]).to_dict()
    return fuel, co2, dates, {"source_holding": summary, "price_statistics": stats}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/hindcast/cost-inputs.yaml"))
    parser.add_argument("--source-directory", type=Path)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--years", type=int, nargs="+")
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    years = sorted(set(args.years or config["years"]))
    source_directory = args.source_directory or Path(config["source_directory"])
    directory = args.directory or Path(config["output_directory"])
    paths = {
        year: {
            "fuel": directory / f"fuel_{year}_hourly.csv",
            "co2": directory / f"co2_{year}_hourly.csv",
            "source_dates": directory / f"cost_source_dates_{year}_hourly.csv",
            "audit": directory / f"cost-inputs_{year}.json",
        } for year in years
    }
    for year_paths in paths.values():
        for path in year_paths.values():
            if path.exists():
                parser.error(f"Refusing to overwrite {path}; select a new --directory")
    sources, provenance = read_sources(source_directory, config, years)
    converted, factors, daily_age, weekly_age = convert(sources, config)
    # Finish every selected year before creating outputs; a missing source or
    # excessively long gap must not leave a seemingly complete annual input set.
    prepared = {
        year: annual(year, converted, config, factors, daily_age, weekly_age)
        for year in years
    }
    directory.mkdir(parents=True, exist_ok=True)
    for year, (fuel, co2, dates, report) in prepared.items():
        year_paths = paths[year]
        for key, frame in [("fuel", fuel), ("co2", co2), ("source_dates", dates)]:
            frame.to_csv(year_paths[key], date_format="%Y-%m-%dT%H:%M:%SZ")
        write_json(year_paths["audit"], {
            "year": year, "hours": len(fuel),
            "units": {"fuel": "EUR/MWh_th, LHV", "co2": "EUR/tCO2"},
            "money_basis": "Nominal observation-year prices; no deflation",
            "time_convention": "Retrospective observation-day UTC blocks; not publication timestamps",
            "config": {"path": str(args.config), "sha256": digest(args.config), "contents": config},
            "input_files": provenance, "conversion_factors": factors,
            "outputs": {
                key: {"path": str(path), "sha256": digest(path)}
                for key, path in year_paths.items() if key != "audit"
            }, **report,
        })
        print(f"{year}: {len(fuel):,} hourly fuel and CO2 values -> {directory}", flush=True)


if __name__ == "__main__":
    main()
