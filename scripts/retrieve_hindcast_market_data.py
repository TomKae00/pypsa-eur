# SPDX-License-Identifier: MIT
"""Retrieve fuel/CO2 source histories without filling or converting their prices.

Run from the repository root. Requires only the existing pandas, PyYAML and
openpyxl dependencies; no market-data account or API key is used. See
doc/hindcast.md for the source conventions and remaining preparation decisions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import date, datetime, timezone
from io import StringIO
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import yaml


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


class Sources:
    def __init__(self, directory):
        self.directory = directory
        directory.mkdir(parents=True, exist_ok=True)
        self.records = []

    def get(self, name, url, suffix, headers=None):
        fingerprint = hashlib.sha256(url.encode()).hexdigest()[:12]
        path = self.directory / f"{name}-{fingerprint}{suffix}"
        receipt = path.with_suffix(path.suffix + ".json")
        if path.exists() and receipt.exists():
            record = json.loads(receipt.read_text())
            if record["url"] != url or record["sha256"] != sha256(path):
                raise ValueError(f"Cached source or checksum changed: {path}")
        else:
            if path.exists():
                raise ValueError(f"Source exists without a receipt: {path}")
            request = Request(
                url,
                headers={"User-Agent": "Mozilla/5.0", **(headers or {})},
            )
            temporary = path.with_suffix(path.suffix + ".part")
            try:
                with urlopen(request, timeout=60) as response:
                    temporary.write_bytes(response.read())
                    record = {
                        "name": name,
                        "url": url,
                        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                        "last_modified": response.headers.get("Last-Modified"),
                    }
                record.update(
                    path=str(path), sha256=sha256(temporary), bytes=temporary.stat().st_size
                )
                # Write the receipt first: an interrupted rename can be resumed.
                write_json(receipt, record)
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        self.records.append(record)
        print(f"Source: {name} ({record['bytes']:,} bytes)", flush=True)
        return path


def source_series(dates, values, name, start, end):
    index = pd.DatetimeIndex(dates)
    if index.hasnans or index.has_duplicates:
        raise ValueError(f"Invalid or duplicate source dates: {name}")
    series = pd.Series(pd.to_numeric(values, errors="raise"), index=index, name=name)
    present = series.dropna()
    if not np.isfinite(present.to_numpy()).all():
        raise ValueError(f"Non-finite source prices: {name}")
    series = series.sort_index().loc[start:end]
    series.index.name = "date"
    if series.dropna().empty:
        raise ValueError(f"No source observations in the requested period: {name}")
    return series


def yahoo(config, sources, start, end):
    p1 = int(pd.Timestamp(start, tz="UTC").timestamp())
    p2 = int((pd.Timestamp(end, tz="UTC") + pd.Timedelta(days=1)).timestamp())
    url = config["endpoint"].format(symbol=quote(config["symbol"], safe=""))
    url += "?" + urlencode({"period1": p1, "period2": p2, "interval": "1d"})
    path = sources.get("gas-yahoo", url, ".json", {"Accept": "application/json"})
    chart = json.loads(path.read_text())["chart"]
    if chart.get("error") or not chart.get("result"):
        raise ValueError(f"Yahoo did not return a gas history: {chart.get('error')}")
    result = chart["result"][0]
    meta = result["meta"]
    if meta.get("symbol") != config["symbol"] or meta.get("currency") != "EUR":
        raise ValueError("Unexpected Yahoo gas symbol or currency")
    # Yahoo's timestamps are exchange-local midnight, including DST changes.
    dates = pd.to_datetime(result["timestamp"], unit="s", utc=True)
    dates = dates.tz_convert(meta["exchangeTimezoneName"]).tz_localize(None).normalize()
    series = source_series(
        dates, result["indicators"]["quote"][0][config["field"]],
        "gas_eur_mwh_quoted", start, end,
    )
    identity = {k: meta.get(k) for k in (
        "symbol", "currency", "shortName", "longName", "exchangeName",
        "instrumentType", "exchangeTimezoneName",
    )}
    return series.to_frame(), identity


def investing(config, sources, name, column, start, end):
    url = config["endpoint"].format(instrument_id=config["instrument_id"])
    url += "?" + urlencode({
        "start-date": start.isoformat(), "end-date": end.isoformat(),
        "time-frame": "Daily", "add-missing-rows": "false",
    })
    path = sources.get(name + "-investing", url, ".json", {
        "Accept": "application/json", "Referer": config["reference"],
        "domain-id": "www",
    })
    rows = json.loads(path.read_text())["data"]
    dates = pd.to_datetime([row["rowDateTimestamp"] for row in rows], utc=True)
    dates = dates.tz_localize(None).normalize()
    series = source_series(
        dates, [row[config["field"]] for row in rows], column, start, end
    )
    return series.to_frame(), {"instrument_id": config["instrument_id"]}


def oil(config, sources, start, end):
    path = sources.get("oil-commission", config["url"], ".xlsx")
    # The first row carries stable machine-readable product/country labels;
    # the next two rows describe the product and its physical unit.
    raw = pd.read_excel(path, sheet_name=config["sheet"], header=0, skiprows=[1, 2])
    # The workbook appends narrative notes and copyright text below the table.
    footer = raw.iloc[:, 0].astype(str).str.strip().eq("Notes:")
    if footer.any():
        raw = raw.iloc[:np.flatnonzero(footer.to_numpy())[0]]
    columns = list(config["columns"])
    if not set(columns).issubset(raw.columns):
        raise ValueError("Commission workbook is missing configured product columns")
    undated = raw.iloc[:, 0].isna()
    if raw.loc[undated, columns].notna().any().any():
        raise ValueError("Oil prices were supplied without an observation date")
    raw = raw.loc[~undated]
    dates = pd.to_datetime(raw.iloc[:, 0], errors="raise")
    result = {}
    for source, settings in config["columns"].items():
        if source not in raw:
            raise ValueError(f"Commission workbook is missing {source}")
        result[settings["name"]] = source_series(
            dates, raw[source].to_numpy(), settings["name"], start, end
        )
    return pd.DataFrame(result), {"sheet": config["sheet"]}


def fx(config, sources, start, end):
    url = config["endpoint"] + "?" + urlencode({
        "startPeriod": start.isoformat(), "endPeriod": end.isoformat(), "format": "csvdata",
    })
    path = sources.get("fx-ecb", url, ".csv", {"Accept": "text/csv"})
    raw = pd.read_csv(StringIO(path.read_text()))
    if not raw["KEY"].eq("EXR.D.USD.EUR.SP00.A").all():
        raise ValueError("ECB returned an unexpected exchange-rate series")
    series = source_series(
        pd.to_datetime(raw["TIME_PERIOD"]), raw["OBS_VALUE"].to_numpy(), "usd_per_eur", start, end
    )
    if (series.dropna() <= 0).any():
        raise ValueError("ECB exchange rates must be positive")
    return series.to_frame(), {"series": "EXR.D.USD.EUR.SP00.A", "quote": "USD per EUR"}


def describe(frame, years):
    summary = {}
    for year in years:
        part = frame.loc[str(year)] if year in frame.index.year else frame.iloc[:0]
        summary[str(year)] = {}
        for column in frame:
            values = part[column].dropna()
            summary[str(year)][column] = {
                "observations": len(values),
                "first_date": str(values.index.min().date()) if len(values) else None,
                "last_date": str(values.index.max().date()) if len(values) else None,
                "missing_source_rows": int(part[column].isna().sum()),
            }
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/hindcast/market-data.yaml")
    parser.add_argument("--years", nargs="+", type=int)
    parser.add_argument("--directory", type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    years = sorted(set(args.years or config["years"]))
    if not years or any(year < 1971 or year > datetime.now(timezone.utc).year - 1 for year in years):
        parser.error("Choose completed historical years from 1971 onwards")
    directory = args.directory or Path(f"data/hindcast/market-data-{years[0]}-{years[-1]}")
    if (directory / "manifest.json").exists():
        parser.error(f"Collection already exists: {directory}. Choose a new --directory to refresh it.")
    # Preserve December observations for later, explicit carry-forward at Jan 1.
    start, end = date(years[0] - 1, 12, 1), date(years[-1], 12, 31)
    sources = Sources(directory / "raw")
    jobs = [
        ("gas_ttf_daily.csv", "gas", lambda: yahoo(config["gas"], sources, start, end)),
        ("coal_rotterdam_daily.csv", "coal", lambda: investing(
            config["coal"], sources, "coal", "coal_usd_t", start, end)),
        ("co2_eua_daily.csv", "co2", lambda: investing(
            config["co2"], sources, "co2", "co2_eur_tco2", start, end)),
        ("oil_eu_weekly.csv", "oil", lambda: oil(config["oil"], sources, start, end)),
        ("fx_daily.csv", "fx", lambda: fx(config["fx"], sources, start, end)),
    ]
    outputs = []
    # Collect and parse all sources before exporting any normalized tables.
    # Failed downloads retain only source caches for a subsequent retry.
    for filename, name, retrieve in jobs:
        frame, identity = retrieve()
        outputs.append((filename, name, frame, identity))
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "requested_years": years,
        "download_start_date": start.isoformat(), "download_end_date": end.isoformat(),
        "config": {"path": args.config, "sha256": sha256(args.config), "contents": config},
        "status": "Source histories only; not generator marginal-cost inputs",
        "processing": {
            "date_semantics": "Trading/observation dates, not intraday publication timestamps",
            "weekend_or_holiday_filling": False, "interpolation": False,
            "currency_conversion": False, "energy_conversion": False,
            "lignite_proxy_applied": False, "co2_rounding": False,
            "money_basis": "Nominal prices in each observation year; no deflation",
        },
        "raw_sources": sources.records, "outputs": [],
    }
    for filename, name, frame, identity in outputs:
        path = directory / filename
        if path.exists():
            raise ValueError(f"Refusing to overwrite {path}")
        frame.to_csv(path, index_label="date", date_format="%Y-%m-%d")
        counts = describe(frame, years)
        manifest["outputs"].append({
            "path": str(path), "sha256": sha256(path),
            "source": name, "identity": identity, "coverage": counts,
        })
        print(f"Exported: {path}", flush=True)
        for year in years:
            print(f"  {year}: " + ", ".join(
                f"{col}={data['observations']} observations" for col, data in counts[str(year)].items()
            ), flush=True)
    write_json(directory / "manifest.json", manifest)
    print(f"Provenance: {directory / 'manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
