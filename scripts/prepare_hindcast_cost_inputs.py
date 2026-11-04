# SPDX-License-Identifier: MIT
"""Prepare hourly fuel/CO2 inputs from a retrieved market-data collection.

Fuel outputs follow the price convention recovered from the 2023 archive;
their physical HHV/LHV basis is not independently established. CO2 uses EUR/tCO2.
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


def observations(dates, values, name, config):
    dates = pd.DatetimeIndex(pd.to_datetime(dates, errors="raise", utc=True))
    if dates.hasnans or dates.has_duplicates or not dates.is_monotonic_increasing:
        raise ValueError(f"Source dates must be unique and increasing: {name}")
    if not dates.equals(dates.normalize()):
        raise ValueError(f"Expected daily observation labels: {name}")
    values = np.asarray(pd.to_numeric(values, errors="raise"), dtype=float)
    if not np.isfinite(values[~np.isnan(values)]).all():
        raise ValueError(f"Non-finite source values: {name}")
    series = pd.Series(values, index=dates, name=name).dropna()
    if series.empty:
        raise ValueError(f"No usable observations: {name}")
    lag = config["time"]["effective_lag_days"][name]
    if isinstance(lag, bool) or not isinstance(lag, int) or lag < 0:
        raise ValueError(f"Effective lag for {name} must be a non-negative integer")
    frame = pd.DataFrame({
        "value": series.to_numpy(), "observation_date": series.index,
    }, index=series.index + pd.Timedelta(days=lag))
    frame.index.name = "effective_date"
    return frame


def weighted_oil(directory, manifest, config):
    spec = config["oil"]
    if spec["missing_consumption"] != "exclude_country" or spec["missing_price"] != "renormalize_available_weights":
        raise ValueError("Unsupported oil missing-data convention")
    candidates = [item for item in manifest["raw_sources"] if item["name"] == spec["workbook_source"]]
    if len(candidates) != 1:
        raise ValueError("Expected one original Commission workbook in the collection")
    record = candidates[0]
    path = directory / "raw" / Path(record["path"]).name
    if digest(path) != record["sha256"]:
        raise ValueError(f"Raw oil workbook differs from its collection manifest: {path}")
    countries = spec["countries"]
    if len(set(countries)) != len(countries):
        raise ValueError("Oil-weight countries must not be duplicated")
    price_columns = [country + spec["price_column_suffix"] for country in countries]
    weight_columns = [country + spec["consumption_column_suffix"] for country in countries]
    prices = pd.read_excel(path, sheet_name=spec["price_sheet"], header=0, skiprows=[1, 2])
    footer = prices.iloc[:, 0].astype(str).str.strip().eq("Notes:")
    if footer.any():
        prices = prices.iloc[:np.flatnonzero(footer.to_numpy())[0]]
    prices = prices.loc[prices.iloc[:, 0].notna()].copy()
    dates = pd.DatetimeIndex(pd.to_datetime(prices.iloc[:, 0], errors="raise", utc=True))
    if not set(price_columns).issubset(prices):
        raise ValueError("Oil workbook is missing configured national price columns")
    prices = prices[price_columns].apply(pd.to_numeric, errors="raise")
    prices.index = dates
    prices = prices.sort_index().loc[
        pd.Timestamp(manifest["download_start_date"], tz="UTC"):
        pd.Timestamp(manifest["download_end_date"], tz="UTC")
    ]
    consumption = pd.read_excel(path, sheet_name=spec["consumption_sheet"], header=0, skiprows=[1, 2])
    weight_year = spec["consumption_weight_year"]
    year_rows = consumption.loc[pd.to_numeric(consumption.iloc[:, 0], errors="coerce").eq(weight_year)]
    if len(year_rows) != 1 or not set(weight_columns).issubset(consumption):
        raise ValueError(f"Expected national oil consumption weights for {weight_year}")
    reported_weights = pd.to_numeric(year_rows.iloc[0][weight_columns], errors="raise").to_numpy(dtype=float)
    weights = np.nan_to_num(reported_weights, nan=0.0, posinf=np.inf, neginf=-np.inf)
    if not np.isfinite(weights).all() or (weights < 0).any() or weights.sum() <= 0:
        raise ValueError("Oil consumption weights must be finite, non-negative and non-zero in total")
    values = prices.to_numpy(dtype=float)
    if not np.isfinite(values[~np.isnan(values)]).all():
        raise ValueError("National oil prices contain non-finite values")
    available = np.isfinite(values) & (weights[None, :] > 0)
    denominator = np.where(available, weights[None, :], 0).sum(axis=1)
    numerator = np.where(available, values * weights[None, :], 0).sum(axis=1)
    aggregate = np.divide(numerator, denominator, out=np.full(len(prices), np.nan), where=denominator > 0)
    frame = observations(prices.index, aggregate, "oil", config)
    frame["countries_used"] = pd.Series(available.sum(axis=1), index=prices.index).reindex(
        pd.DatetimeIndex(frame["observation_date"])
    ).to_numpy()
    provenance = {
        "path": str(path), "sha256": digest(path), "source_url": record["url"],
        "retrieved_at_utc": record["retrieved_at_utc"], "price_sheet": spec["price_sheet"],
        "consumption_sheet": spec["consumption_sheet"], "weight_year": int(weight_year),
        "national_consumption_weights": {
            country: float(weight) for country, weight in zip(countries, weights)
            if weight > 0
        },
        "countries_without_reported_consumption": [
            country for country, weight in zip(countries, reported_weights) if np.isnan(weight)
        ],
        "usable_observations": len(frame), "missing_source_rows": int(np.isnan(aggregate).sum()),
        "aggregation": "Sum(price * consumption) / sum(consumption), over available positively weighted countries",
    }
    return frame, provenance


def read_sources(directory, config, years):
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if not set(years).issubset(manifest["requested_years"]):
        raise ValueError("Selected years are not covered by the source collection")
    records = {Path(item["path"]).name: item for item in manifest["outputs"]}
    result, provenance = {}, {}
    for name in ["gas", "coal", "co2"]:
        spec = config[name]
        path = directory / spec["file"]
        if path.name not in records or digest(path) != records[path.name]["sha256"]:
            raise ValueError(f"Source file does not match its collection manifest: {path}")
        frame = pd.read_csv(path)
        values = pd.to_numeric(frame[spec["column"]], errors="raise").to_numpy(dtype=float)
        result[name] = observations(frame["date"], values, name, config)
        provenance[name] = {
            "path": str(path), "sha256": digest(path), "column": spec["column"],
            "usable_observations": len(result[name]), "missing_source_rows": int(np.isnan(values).sum()),
        }
    result["oil"], provenance["oil"] = weighted_oil(directory, manifest, config)
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
    if config["procedure"] != "archive_2023_inferred":
        raise ValueError("Unsupported cost-input procedure")
    if config["time"]["holding"] != "previous_observation":
        raise ValueError("Only previous_observation holding is implemented")
    daily_age = positive(config["time"]["daily_max_age_days"], "Daily maximum age")
    weekly_age = positive(config["time"]["weekly_max_age_days"], "Weekly maximum age")
    gas_factor = positive(config["gas"]["price_multiplier"], "Gas multiplier")
    co2_factor = positive(config["co2"]["price_multiplier"], "CO2 multiplier")
    decimals = config["coal"]["source_price_decimals"]
    if isinstance(decimals, bool) or not isinstance(decimals, int) or decimals < 0:
        raise ValueError("Coal source-price decimals must be a non-negative integer")
    coal_factor = positive(config["coal"]["combined_price_factor"], "Coal combined conversion factor")
    oil_divisor = positive(config["oil"]["eur_per_tonne_divisor"], "Oil energy-price divisor")
    lignite_factor = positive(config["lignite"]["coal_price_factor"], "Lignite factor")

    # This is the empirically recovered combined coal factor. Applying ECB
    # rates here too would introduce another currency conversion.
    converted = {name: frame.copy() for name, frame in sources.items()}
    converted["gas"]["value"] *= gas_factor
    converted["coal"]["value"] = converted["coal"]["value"].round(decimals) * coal_factor
    converted["oil"]["value"] /= oil_divisor
    converted["co2"]["value"] *= co2_factor
    factors = {
        "gas_price_multiplier": gas_factor, "co2_price_multiplier": co2_factor,
        "coal_source_price_decimals": decimals, "coal_combined_price_factor": coal_factor,
        "additional_ecb_conversion": False, "oil_eur_per_tonne_divisor": oil_divisor,
        "oil_consumption_weight_year": config["oil"]["consumption_weight_year"],
        "lignite_coal_energy_price_ratio": lignite_factor,
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
        if name == "oil":
            dates["oil_countries_used"] = frame["countries_used"].astype(int)
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
    print(
        "Using the procedure inferred from the 2023 archive. "
        "CO2 provider and eight reference oil dates remain unresolved; "
        "parameters/weights in other years are explicit extensions.", flush=True,
    )
    for year, (fuel, co2, dates, report) in prepared.items():
        year_paths = paths[year]
        for key, frame in [("fuel", fuel), ("co2", co2), ("source_dates", dates)]:
            frame.to_csv(year_paths[key], date_format="%Y-%m-%dT%H:%M:%SZ")
        write_json(year_paths["audit"], {
            "year": year, "hours": len(fuel),
            "units": {
                "fuel": "EUR/MWh_th, implemented archive price convention",
                "fuel_basis_note": "Physical HHV/LHV basis not independently established for all quotes/factors",
                "co2": "EUR/tCO2",
            },
            "money_basis": "Nominal observation-year prices; no deflation",
            "time_convention": "Retrospective observation-day UTC blocks; not publication timestamps",
            "config": {"path": str(args.config), "sha256": digest(args.config), "contents": config},
            "input_files": provenance, "conversion_factors": factors,
            "alignment": {
                **config["alignment"],
                "year_has_verified_parameters": year in config["alignment"]["verified_parameter_years"],
                "co2_source_status": config["co2"]["alignment_status"],
                "oil_unresolved_intervals_in_this_year": [
                    interval for interval in config["alignment"]["oil_unresolved_reference_dates"]
                    if pd.Timestamp(interval["start"]).year == year
                ],
                "manual_oil_patches_applied": False,
            },
            "outputs": {
                key: {"path": str(path), "sha256": digest(path)}
                for key, path in year_paths.items() if key != "audit"
            }, **report,
        })
        print(f"{year}: {len(fuel):,} hourly fuel and CO2 values -> {directory}", flush=True)


if __name__ == "__main__":
    main()
