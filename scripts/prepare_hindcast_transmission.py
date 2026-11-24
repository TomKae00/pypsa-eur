# SPDX-License-Identifier: MIT
"""Rebuild the 2024/2025 input tables from the Nordic TSO and JAO publications."""

import argparse
import gzip
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from time import sleep
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

JAO = "https://publicationtool.jao.eu/nordic/api/data/maxBorderFlow"
NTC = "https://www.nordpoolgroup.com/globalassets/download-center/tso/max-ntc.pdf"
TTC = "https://www.fingrid.fi/globalassets/dokumentit/fi/kantaverkko/sahkonsiirto/appendix-1-9-to-soa-annex-os-v3.pdf"
SWITCH = pd.Timestamp(
    "2024-10-29T23:00:00Z"
)  # First delivery at 00:00 CET on 30 October.
QUARTER_HOUR = pd.Timestamp("2025-09-30T22:00:00Z")  # 1 October, 00:00 CEST.
YEARS = (2024, 2025)


def download(start, stop, cache):
    name = start.strftime("%Y%m%d-%H") + ".json"
    target = cache / name
    if not target.exists():
        query = urlencode(
            {
                "FromUtc": start.isoformat().replace("+00:00", "Z"),
                "ToUtc": stop.isoformat().replace("+00:00", "Z"),
                "skip": 0,
                "take": 500,
            }
        )
        request = Request(
            f"{JAO}?{query}", headers={"User-Agent": "PyPSA-Eur hindcast research"}
        )
        for attempt in range(4):
            try:
                with urlopen(request, timeout=30) as reply:
                    data = reply.read()
                break
            except (HTTPError, URLError, TimeoutError):
                if attempt == 3:
                    raise
                sleep(2 * (attempt + 1))
        validate = json.loads(data)
        if validate.get("rejected") or not validate.get("data"):
            raise ValueError(f"JAO returned no records for {start} to {stop}")
        target.write_bytes(data)
    return target.read_bytes()


def retrieve(cache):
    stop = pd.Timestamp("2026-01-01T00:00:00Z")
    boundaries = [*pd.date_range(SWITCH, stop, freq="2D"), stop]
    spans = list(zip(boundaries[:-1], boundaries[1:], strict=True))
    with ThreadPoolExecutor(max_workers=4) as pool:
        payloads = list(pool.map(lambda p: download(*p, cache), spans))
    records = [row for data in payloads for row in json.loads(data)["data"]]
    return (
        pd.DataFrame.from_records(records),
        {
            start.isoformat(): hashlib.sha256(data).hexdigest()
            for (start, _), data in zip(spans, payloads, strict=True)
        },
    )


def prepare(records, corridors, directory):
    dates = pd.DatetimeIndex(pd.to_datetime(records.dateTimeUtc, utc=True))
    finish = pd.Timestamp("2026-01-01T00:00:00Z")
    expected = pd.date_range(SWITCH, QUARTER_HOUR, freq="h", inclusive="left").append(
        pd.date_range(QUARTER_HOUR, finish, freq="15min", inclusive="left")
    )
    if not dates.sort_values().as_unit("ns").equals(expected.as_unit("ns")):
        missing = expected.difference(dates)
        extra = dates.difference(expected)
        raise ValueError(
            f"JAO snapshot coverage mismatch: {len(missing)} missing, "
            f"{len(extra)} unexpected, {dates.duplicated().sum()} duplicates"
        )
    dates = dates.floor("h")
    annual_rows = []
    for year in YEARS:
        hours = pd.date_range(
            f"{year}-01-01", f"{year + 1}-01-01", freq="h", inclusive="left", tz="UTC"
        )
        rows = []
        for c in corridors.itertuples(index=False):
            fields = [
                f"border_{c.jao_zone0}_{c.jao_zone1}",
                f"border_{c.jao_zone1}_{c.jao_zone0}",
            ]
            raw = records[fields].apply(pd.to_numeric, errors="raise").set_axis(dates)
            if not np.isfinite(raw.to_numpy()).all() or (raw < 0).any().any():
                raise ValueError(f"Invalid JAO direction values on {c.zone0}–{c.zone1}")
            # The model holds one mean MW value across the hour. Averaging the
            # quarter-hour ceilings preserves their total permissible hourly energy.
            bounds = raw.groupby(level=0).mean().reindex(hours)
            bounds.loc[hours < SWITCH, :] = [
                c.forward_pre_switch_MW,
                c.reverse_pre_switch_MW,
            ]
            if bounds.isna().any().any():
                raise ValueError(f"Incomplete hourly bounds on {c.zone0}–{c.zone1}")
            source = JAO
            basis = "Hourly JAO Max Border Flow; independent directional envelopes"
            if year == 2024:
                ttc = c.zone0 == "FI" and c.zone1 == "NO4"
                source = f"{TTC if ttc else NTC} | {JAO}"
                before = (
                    "published total transfer envelope"
                    if ttc
                    else "Nordic TSO maximum NTC"
                )
                basis += f"; {before} before 30 October local delivery"
            annual_rows.append((year, c.zone0, c.zone1, *bounds.max(), basis, source))
            rows.append(
                bounds.set_axis(["forward_MW", "reverse_MW"], axis=1).assign(
                    zone0=c.zone0, zone1=c.zone1
                )
            )
        hourly = pd.concat(rows).loc[:, ["zone0", "zone1", "forward_MW", "reverse_MW"]]
        with (directory / f"hourly_{year}.csv.gz").open("wb") as stream:
            with gzip.GzipFile(
                filename="", fileobj=stream, mode="wb", mtime=0
            ) as compressed:
                compressed.write(
                    hourly.to_csv(
                        index_label="snapshot",
                        float_format="%.2f",
                        date_format="%Y-%m-%dT%H:%M:%SZ",
                    ).encode()
                )
        print(f"{year}: {len(hours):,} UTC hours, {len(rows)} corridors")
    annual = pd.DataFrame(
        annual_rows,
        columns=[
            "planning_year",
            "zone0",
            "zone1",
            "forward_MW",
            "reverse_MW",
            "basis",
            "source",
        ],
    )
    annual.to_csv(directory / "limits.csv", index=False, float_format="%.2f")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--directory", type=Path, default=Path("data/hindcast/transmission")
    )
    parser.add_argument(
        "--cache", type=Path, required=True, help="Directory to cache raw JAO pages"
    )
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    records, hashes = retrieve(args.cache)
    prepare(
        records, pd.read_csv(args.directory / "nordic_corridors.csv"), args.directory
    )
    manifest = {
        "api": JAO,
        "mtu_records": len(records),
        "flow_based_from_utc": SWITCH.isoformat(),
        "quarter_hour_from_utc": QUARTER_HOUR.isoformat(),
        "raw_response_sha256_by_query_from_utc": hashes,
    }
    (args.directory / "source_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
