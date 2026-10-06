# SPDX-License-Identifier: MIT
"""
Retrieve directional ENTSO-E flows, preserving missing observations.

ENTSOE_API_KEY must be set in the environment. Example:
  python scripts/retrieve_hindcast_flows.py --year 2023
Physical flows and day-ahead commercial exchanges are saved separately.
"""

import argparse
import os
from pathlib import Path

import pandas as pd
import yaml
from hindcast import digest, utc_index, write_json


def hourly_complete(series, start, end):
    """Average regular subhourly MW only when every interval in an hour is present."""
    if not isinstance(series, pd.Series):
        raise ValueError("Expected a single ENTSO-E flow series")
    series = pd.to_numeric(series, errors="raise").copy()
    series.index = utc_index(series.index)
    series = series.loc[(series.index >= start) & (series.index < end)]
    if len(series) < 2:
        raise ValueError("Insufficient observations to determine flow resolution")
    deltas = series.index.to_series().diff().dropna()
    step = deltas.min()
    if step not in [pd.Timedelta(minutes=x) for x in [15, 30, 60]]:
        raise ValueError(f"Unsupported flow interval {step}")
    per_hour = int(pd.Timedelta(hours=1) / step)
    grid = pd.date_range(start, end, freq=step, inclusive="left")
    aligned = series.reindex(grid)
    mean = aligned.resample("h").mean()
    count = aligned.resample("h").count()
    return mean.where(count == per_hour)


def retrieve(client, cfg, year, kind):
    start = pd.Timestamp(f"{year}-01-01", tz="UTC")
    end = pd.Timestamp(f"{year + 1}-01-01", tz="UTC")
    index = pd.date_range(start, end, freq="h", inclusive="left")
    output = {}
    for border, pairs in cfg.items():
        contributions = []
        for origin, destination in pairs:
            directional = []
            for a, b in [(origin, destination), (destination, origin)]:
                if kind == "physical":
                    raw = client.query_crossborder_flows(a, b, start=start, end=end)
                else:
                    raw = client.query_scheduled_exchanges(
                        a, b, start=start, end=end, dayahead=True
                    )
                # NoMatchingDataError is deliberately NOT converted to a zero series.
                directional.append(hourly_complete(raw, start, end).reindex(index))
            contributions.append(directional[0] - directional[1])
        if not contributions:
            raise ValueError(f"No source pairs for {border}")
        parts = pd.concat(contributions, axis=1)
        output[border] = parts.sum(axis=1, min_count=len(contributions))
    return pd.DataFrame(output, index=index)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=2023)
    parser.add_argument(
        "--kind", choices=["physical", "scheduled_day_ahead"], default="physical"
    )
    parser.add_argument("--config", default="data/hindcast/flow-sources-zonal-2023.yaml")
    parser.add_argument("--output")
    args = parser.parse_args()
    key = os.environ.get("ENTSOE_API_KEY")
    if not key:
        parser.error(
            "Set ENTSOE_API_KEY in your environment; do not put it in a tracked file"
        )
    from entsoe import EntsoePandasClient

    cfg = yaml.safe_load(Path(args.config).read_text())
    output = Path(args.output or f"data/hindcast/flows_{args.kind}_{args.year}.csv")
    if output.exists():
        parser.error(f"{output} exists; choose another --output to refresh")
    frame = retrieve(EntsoePandasClient(api_key=key), cfg, args.year, args.kind)
    if frame.isna().all().any():
        raise ValueError("At least one border has no complete paired observations")
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index_label="utc_timestamp")
    write_json(
        output.with_suffix(".json"),
        {
            "kind": args.kind,
            "year": args.year,
            "sources": cfg,
            "coverage": frame.notna().mean().to_dict(),
            "sha256": digest(output),
            "unit": "MW",
            "direction": "first region -> second region",
            "missing_data": "no filling; missing directional query is an error",
        },
    )
    print(output)


if __name__ == "__main__":
    main()
