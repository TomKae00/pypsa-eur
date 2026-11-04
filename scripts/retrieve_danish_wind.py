# SPDX-License-Identifier: MIT
"""
Retrieve hourly DK1/DK2 wind production, prices and monthly fleet capacity.

Source: Energinet Energi Data Service, https://www.energidataservice.dk/
ProductionConsumptionSettlement includes the 2023 park-production reporting
change; no individual-turbine production series is spliced across that change.
Monthly reported capacity is an approximation to intra-month commissioning.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

API = "https://api.energidataservice.dk/dataset/"


def municipality_area(number):
    number = int(number)
    if 101 <= number <= 400:
        return "DK2"
    if 410 <= number <= 999:
        return "DK1"
    return None


def fetch(dataset, year, directory):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{dataset}.json"
    if path.exists():
        return json.loads(path.read_text()), hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
    params = {
        "start": str(pd.Timestamp(year, 1, 1) - pd.Timedelta(days=1))[:10],
        "end": str(pd.Timestamp(year + 1, 1, 1) + pd.Timedelta(days=1))[:10],
        "limit": 0,
    }
    if dataset != "CapacityPerMunicipality":
        params["filter"] = json.dumps({"PriceArea": ["DK1", "DK2"]})
    for attempt in range(5):
        response = requests.get(API + dataset, params=params, timeout=60)
        if response.status_code == 429:
            delay = float(response.headers.get("Retry-After", 10))
            if delay > 300 or attempt == 4:
                raise RuntimeError("Energinet rate limit; retry retrieval later")
            while delay > 0:
                time.sleep(min(30, delay))
                delay -= 30
            continue
        response.raise_for_status()
        data = response.json()
        if not data.get("records"):
            raise ValueError("Empty Energinet dataset: " + dataset)
        temporary = path.with_suffix(".tmp")
        temporary.write_bytes(response.content)
        temporary.replace(path)
        return data, hashlib.sha256(response.content).hexdigest()
    raise RuntimeError("Energinet retrieval failed")


def hourly(data, year, label):
    frame = pd.DataFrame(data["records"])
    frame["snapshot"] = pd.to_datetime(frame.HourUTC, utc=True)
    start, end = (
        pd.Timestamp(f"{year}-01-01", tz="UTC"),
        pd.Timestamp(f"{year + 1}-01-01", tz="UTC"),
    )
    frame = frame[
        (frame.snapshot >= start)
        & (frame.snapshot < end)
        & frame.PriceArea.isin(["DK1", "DK2"])
    ]
    expected = pd.date_range(start, end, inclusive="left", freq="h")
    for area in ["DK1", "DK2"]:
        idx = pd.DatetimeIndex(
            frame.loc[frame.PriceArea.eq(area), "snapshot"]
        ).sort_values()
        if not idx.equals(expected):
            raise ValueError(f"Incomplete or duplicate UTC hours: {label} {area}")
    return frame.set_index(["snapshot", "PriceArea"]).sort_index()


def retrieve(year, output, cache=None):
    output = Path(output)
    cache = Path(cache) if cache else output.parent / "sources" / str(year)
    sources, receipts = {}, {}
    for dataset in [
        "ProductionConsumptionSettlement",
        "CapacityPerMunicipality",
        "Elspotprices",
    ]:
        data, sha = fetch(dataset, year, cache)
        sources[dataset] = data
        receipts[dataset] = {"url": API + dataset, "raw_sha256": sha}
        time.sleep(2)
    production = hourly(sources["ProductionConsumptionSettlement"], year, "production")
    prices = hourly(sources["Elspotprices"], year, "prices")
    capacities = pd.DataFrame(sources["CapacityPerMunicipality"]["records"])
    capacities["month"] = pd.to_datetime(capacities.Month).dt.to_period("M").astype(str)
    capacities = capacities[capacities.month.str.startswith(str(year))].copy()
    if capacities.duplicated(["month", "MunicipalityNo"]).any():
        raise ValueError("Duplicate municipality capacity record")
    capacities["PriceArea"] = capacities.MunicipalityNo.map(municipality_area)
    wind_columns = ["OnshoreWindCapacity", "OffshoreWindCapacity"]
    if (
        capacities.loc[capacities.PriceArea.isna(), wind_columns]
        .fillna(0)
        .to_numpy()
        .sum()
        > 0
    ):
        raise ValueError(
            "Unmapped municipality with wind capacity; review DK1/DK2 crosswalk"
        )
    totals = capacities.groupby(["month", "PriceArea"])[wind_columns].sum(min_count=1)
    rows = []
    for area in ["DK1", "DK2"]:
        observed = production.xs(area, level="PriceArea")
        spot = prices.xs(area, level="PriceArea").SpotPriceEUR.reindex(observed.index)
        for carrier, prefix, capacity_column in [
            ("onwind", "OnshoreWind", "OnshoreWindCapacity"),
            ("offwind", "OffshoreWind", "OffshoreWindCapacity"),
        ]:
            fields = (
                [prefix + "Lt50kW_MWh", prefix + "Ge50kW_MWh"]
                if carrier == "onwind"
                else [prefix + "Lt100MW_MWh", prefix + "Ge100MW_MWh"]
            )
            generation = observed[fields].sum(axis=1, min_count=len(fields))
            monthly = totals.xs(area, level="PriceArea")[capacity_column]
            fleet = observed.index.strftime("%Y-%m").map(monthly)
            values = np.column_stack([generation, fleet, spot])
            if (
                not np.isfinite(values).all()
                or (values[:, :2] < 0).any()
                or (values[:, 1] <= 0).any()
            ):
                raise ValueError(f"Invalid wind observations: {area} {carrier}")
            rows.append(
                pd.DataFrame(
                    {
                        "snapshot": observed.index,
                        "zone": area,
                        "carrier": carrier,
                        "generation_MW": generation.to_numpy(),
                        "capacity_MW": fleet,
                        "spot_price_EUR_MWh": spot.to_numpy(),
                    }
                )
            )
    frame = pd.concat(rows, ignore_index=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False, compression={"method": "gzip", "mtime": 0})
    receipt = {
        "year": year,
        "sources": receipts,
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "rows": len(frame),
        "hours_per_series": len(production) // 2,
        "capacity_definition": "Monthly reported fleet, municipality DK1/DK2 allocation; no intra-month commissioning dates inferred",
    }
    output.with_suffix(".json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake = globals()["snakemake"]
        retrieve(int(snakemake.wildcards.year), snakemake.output.csv)
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--year", type=int, required=True)
        parser.add_argument("--output", type=Path, required=True)
        parser.add_argument("--cache", type=Path)
        args = parser.parse_args()
        print(json.dumps(retrieve(args.year, args.output, args.cache), indent=2))
