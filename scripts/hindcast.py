# SPDX-License-Identifier: MIT
"""
Reproduce historical dispatch inputs and benchmark raw PyPSA prices and flows.

Run from the PyPSA-Eur repository root. See doc/hindcast.md for assumptions.
This CLI is independent of the capacity-expansion Snakemake solve target.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

LOG = logging.getLogger(__name__)
ARCHIVE = "https://api.figshare.com/v2/articles/31248568/versions/1"
ANALYZER_COMMIT = "fd282cc8b1705b68aad0edd286d6b8f1c612ac5d"
ANALYZER_RAW = f"https://raw.githubusercontent.com/marco-saretta/pypsa-network-analyzer/{ANALYZER_COMMIT}"
CAPACITIES = {
    "generators": "p_nom",
    "links": "p_nom",
    "storage_units": "p_nom",
    "stores": "e_nom",
    "lines": "s_nom",
    "transformers": "s_nom",
}


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, default=str) + "\n")


def download(url, path, md5=None):
    """Atomic download; existing reference files are checked rather than overwritten."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        tmp = path.with_suffix(path.suffix + ".part")
        try:
            with (
                urllib.request.urlopen(url, timeout=60) as response,
                tmp.open("wb") as f,
            ):
                shutil.copyfileobj(response, f)
            if md5 and digest(tmp, "md5") != md5:
                raise ValueError(f"Archive checksum mismatch: {path}")
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)
    if md5 and digest(path, "md5") != md5:
        raise ValueError(f"Existing file has wrong checksum: {path}")
    return {"path": str(path), "url": url, "sha256": digest(path)}


def fetch_reference(year, directory):
    directory = Path(directory)
    with urllib.request.urlopen(ARCHIVE, timeout=60) as response:
        article = json.load(response)
    name = f"hindcast-dyn-wy{year}.nc"
    files = [f for f in article["files"] if f["name"] == name]
    if len(files) != 1:
        raise ValueError(f"Expected exactly one {name} in the pinned archive")
    f = files[0]
    records = [
        download(f["download_url"] + "?download=1", directory / name, f["computed_md5"])
    ]
    for name in ["electricity_prices.csv", "demand.csv"]:
        records.append(
            download(f"{ANALYZER_RAW}/data/benchmark/{name}", directory / name)
        )
    write_json(
        directory / f"sources-{year}.json",
        {"archive": ARCHIVE, "analyzer_commit": ANALYZER_COMMIT, "files": records},
    )


def utc_index(index):
    """Network naive timestamps are explicitly interpreted as UTC."""
    result = pd.DatetimeIndex(pd.to_datetime(index, utc=True))
    if result.has_duplicates:
        raise ValueError("Duplicate UTC timestamps")
    if not result.is_monotonic_increasing:
        raise ValueError("Timestamps must be sorted")
    return result


def read_hourly(path):
    df = pd.read_csv(path, index_col=0)
    df.index = utc_index(df.index)
    if not df.index.equals(df.index.floor("h")):
        raise ValueError(
            f"{path}: expected hourly timestamps; resample source data explicitly"
        )
    return df.apply(pd.to_numeric, errors="raise").replace([np.inf, -np.inf], np.nan)


def hourly_network(n):
    idx = utc_index(n.snapshots)
    # Missing Feb 29 is permitted, but typical-day or multi-hour networks are not.
    deltas = idx[1:] - idx[:-1]
    gap = deltas != pd.Timedelta(hours=1)
    if any(gap):
        allowed = (
            (idx[:-1].month == 2)
            & (idx[:-1].day == 28)
            & (idx[:-1].hour == 23)
            & (idx[1:].month == 3)
            & (idx[1:].day == 1)
            & (deltas == pd.Timedelta(hours=25))
        )
        if any(gap & ~allowed):
            raise ValueError(
                "Expected chronological hourly snapshots (optional leap-day omission)"
            )
    weights = n.snapshot_weightings
    if not np.allclose(weights.to_numpy(), 1.0):
        raise ValueError("Historical benchmark requires one-hour snapshot weights")
    return idx


def audit(n):
    idx = hourly_network(n)
    extension = {}
    for table, attr in CAPACITIES.items():
        df = getattr(n, table)
        extension[table] = df.index[df[f"{attr}_extendable"]].tolist()
    dynamic = n.generators_t.marginal_cost
    return {
        "snapshots": len(idx),
        "start_utc": str(idx[0]),
        "end_utc": str(idx[-1]),
        "buses": len(n.buses),
        "extendable": extension,
        "time_varying_cost_generators": int((dynamic.nunique() > 1).sum()),
        "global_constraints": n.global_constraints.reset_index().to_dict("records"),
        "network_meta": n.meta,
    }


def region_buses(n, regions):
    mapping = {}
    seen = set()
    for name, spec in regions.items():
        if "buses" in spec:
            buses = pd.Index(spec["buses"])
        elif "country" in spec:
            buses = n.buses.index[
                (n.buses.country == spec["country"])
                & n.buses.carrier.isin(["AC", "DC", ""])
            ]
        else:
            raise ValueError(f"{name}: provide explicit buses or country")
        if buses.empty or not buses.isin(n.buses.index).all():
            raise ValueError(f"{name}: empty or unknown bus mapping")
        if seen.intersection(buses):
            raise ValueError(f"{name}: overlapping region mappings")
        mapping[name] = list(buses)
        seen.update(buses)
    return mapping


def zonal_configuration(n, cfg):
    """Map the entire electricity network using explicit historical zone codes.

    A bus cannot be split into multiple bidding zones by postprocessing.
    """
    buses = n.buses.index[n.buses.carrier.isin(["AC", "DC", ""])]
    table = pd.read_csv(cfg["bus_zones"], dtype=str, keep_default_na=False)
    if not {"bus", "bidding_zone"}.issubset(table.columns):
        raise ValueError("bus_zones CSV requires bus,bidding_zone columns")
    if table.bus.duplicated().any() or table.bidding_zone.str.strip().eq("").any():
        raise ValueError("Each bus needs exactly one explicit bidding zone")
    if set(table.bus) != set(buses):
        raise ValueError(
            "Zone mapping must cover every electricity bus exactly: "
            f"missing={sorted(set(buses) - set(table.bus))}, "
            f"unknown={sorted(set(table.bus) - set(buses))}"
        )
    missing = set(cfg.get("expected_countries", [])) - set(n.buses.loc[buses, "country"])
    if missing:
        raise ValueError(f"European study network is missing countries: {sorted(missing)}")
    missing = set(cfg.get("required_zones", [])) - set(table.bidding_zone)
    if missing:
        raise ValueError(f"Required separate bidding zones are missing: {sorted(missing)}")
    regions = {
        zone: {"buses": group.bus.tolist(), "price_zones": [zone]}
        for zone, group in table.groupby("bidding_zone", sort=True)
    }
    zone_of = table.set_index("bus").bidding_zone.to_dict()
    borders = set()
    for frame in [n.lines, n.links.loc[n.links.carrier.eq("DC")]]:
        for branch in frame.itertuples():
            a, b = zone_of.get(branch.bus0), zone_of.get(branch.bus1)
            if a is not None and b is not None and a != b:
                borders.add(tuple(sorted([a, b])))
    return dict(cfg, regions=regions, borders=sorted(borders), spreads=sorted(borders))


def prepare_zonal(network, config, flow_sources):
    """Create an explicit bus mapping template, then derive all flow query pairs."""
    import pypsa

    n = pypsa.Network(network)
    cfg = yaml.safe_load(Path(config).read_text())
    path = Path(cfg["bus_zones"])
    if not path.exists():
        buses = n.buses.loc[n.buses.carrier.isin(["AC", "DC", ""]), ["country"]].copy()
        buses["bidding_zone"] = ""
        path.parent.mkdir(parents=True, exist_ok=True)
        buses.to_csv(path, index_label="bus")
        print(f"Created {path}. Fill bidding_zone with historical ENTSO-E area codes, then repeat this command.")
        return
    cfg = zonal_configuration(n, cfg)
    path = Path(flow_sources)
    path.parent.mkdir(parents=True, exist_ok=True)
    sources = {f"{a}->{b}": [[a, b]] for a, b in cfg["borders"]}
    path.write_text(yaml.safe_dump(sources, sort_keys=True))
    print(f"Wrote {len(sources)} inter-zone flow query pairs to {path}")


def model_prices(n, mapping):
    price = n.buses_t.marginal_price
    load = n.get_switchable_as_dense("Load", "p_set").T.groupby(n.loads.bus).sum().T
    out = {}
    for name, buses in mapping.items():
        p = price.reindex(index=n.snapshots, columns=buses)
        if len(buses) == 1:
            out[name] = p.iloc[:, 0]
        else:
            d = load.reindex(index=n.snapshots, columns=buses, fill_value=0)
            if (d < 0).any().any():
                raise ValueError(
                    "Negative loads require an explicit price-weighting choice"
                )
            valid = p.notna().all(axis=1) & d.notna().all(axis=1) & (d.sum(axis=1) > 0)
            out[name] = ((p * d).sum(axis=1) / d.sum(axis=1)).where(valid)
    result = pd.DataFrame(out)
    result.index = hourly_network(n)
    return result


def reference_prices(prices, loads, regions):
    """Require all declared zones. Missing weights never silently renormalise a country."""
    result = {}
    for name, spec in regions.items():
        zones = spec["price_zones"]
        p = prices[zones]  # fail on unknown columns
        if len(zones) == 1:
            result[name] = p.iloc[:, 0]
        else:
            if loads is None:
                raise ValueError(
                    f"{name}: load data required for multi-zone aggregation"
                )
            d = loads[zones].reindex(p.index)
            valid = p.notna().all(axis=1) & d.notna().all(axis=1)
            valid &= (d >= 0).all(axis=1) & (d.sum(axis=1) > 0)
            result[name] = ((p * d).sum(axis=1) / d.sum(axis=1)).where(valid)
    return pd.DataFrame(result)


def border_flows(n, mapping, borders):
    """
    Net MW leaving the FIRST region, using both branch ends to respect losses.

    A->B: p0 if bus0 is A; p1 if bus1 is A. p1 is normally negative
    for a branch importing into A. Only electrical Lines and DC Links count.
    """
    out = {}
    for a, b in borders:
        series = []
        for table in ["lines", "links"]:
            df = getattr(n, table)
            if table == "links":
                df = df[df.carrier.eq("DC")]
            forward = df.index[df.bus0.isin(mapping[a]) & df.bus1.isin(mapping[b])]
            reverse = df.index[df.bus1.isin(mapping[a]) & df.bus0.isin(mapping[b])]
            ts = getattr(n, table + "_t")
            for ids, port in [(forward, "p0"), (reverse, "p1")]:
                if len(ids):
                    values = ts[port].reindex(index=n.snapshots, columns=ids)
                    series.append(values.sum(axis=1, min_count=len(ids)))
        if not series:
            raise ValueError(
                f"No electrical branches found for {a}->{b}; check mappings"
            )
        out[f"{a}->{b}"] = pd.concat(series, axis=1).sum(axis=1, min_count=len(series))
    result = pd.DataFrame(out, index=n.snapshots)
    result.index = hourly_network(n)
    return result


def load_shedding(n):
    ids = n.generators.index[
        n.generators.carrier.str.contains("load|shedding", case=False, na=False)
    ]
    # Historical PyPSA-Eur files may use sign=0.001 for kW load-shedding units.
    p = n.generators_t.p.reindex(index=n.snapshots, columns=ids)
    p = p.mul(n.generators.loc[ids, "sign"], axis=1).clip(lower=0)
    result = p.T.groupby(n.generators.loc[ids, "bus"]).sum().T
    result.index = hourly_network(n)
    return result


def paired_metrics(sim, obs, frequency="hourly", min_coverage=1.0):
    """Align before aggregation; both means use exactly the same observed hours."""
    pair = pd.concat(
        [sim.rename("model"), obs.reindex(sim.index).rename("observed")], axis=1
    )
    pair = pair.replace([np.inf, -np.inf], np.nan)
    valid = pair.notna().all(axis=1)
    n_expected, n_valid = len(pair), int(valid.sum())
    matched = pair.where(valid)
    if frequency != "hourly":
        rule = {"daily": "D", "weekly": "W-MON"}[frequency]
        coverage = valid.resample(rule).mean()
        matched = matched.resample(rule).mean().loc[coverage >= min_coverage]
    matched = matched.dropna()
    if matched.empty:
        raise ValueError("No valid matched observations at requested resolution")
    error = matched.model - matched.observed
    denom = matched.model.abs() + matched.observed.abs()
    smape = np.divide(
        200 * error.abs(), denom, out=np.zeros(len(error)), where=denom.to_numpy() != 0
    )
    metrics = {
        "n_expected_hours": n_expected,
        "n_matched_hours": n_valid,
        "coverage": n_valid / n_expected,
        "n_scored": len(matched),
        "model_mean": matched.model.mean(),
        "observed_mean": matched.observed.mean(),
        "bias": error.mean(),
        "mae": error.abs().mean(),
        "rmse": np.sqrt((error**2).mean()),
        "smape_percent": np.mean(smape),
        "correlation": matched.model.corr(matched.observed)
        if len(matched) > 1
        else np.nan,
    }
    return metrics, matched


def plot_pair(pair, path, unit):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(11, 6), layout="constrained")
    pair.resample("W-MON").mean().plot(ax=axes[0], color=["#007f8b", "#df8627"])
    axes[0].set(title=path.stem + " — weekly means", ylabel=unit, xlabel="UTC")
    for col, color in zip(pair.columns, ["#007f8b", "#df8627"]):
        axes[1].plot(np.sort(pair[col].dropna())[::-1], label=col, color=color)
    axes[1].set(title="Duration curves — matched hours", xlabel="Hours", ylabel=unit)
    axes[1].legend()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def benchmark(network, cfg, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    import pypsa

    n = pypsa.Network(network)
    report = audit(n)
    if cfg.get("bus_zones"):
        cfg = zonal_configuration(n, cfg)
        report["bus_zones_sha256"] = digest(cfg["bus_zones"])
    mapping = region_buses(n, cfg["regions"])
    model = model_prices(n, mapping)
    year = cfg["year"]
    model = model.loc[model.index.year == year]
    if model.empty:
        raise ValueError(f"Network has no snapshots in {year}")
    obs = reference_prices(
        read_hourly(cfg["prices"]),
        read_hourly(cfg["loads"]) if cfg.get("loads") else None,
        cfg["regions"],
    )
    report.update(
        {
            "network_sha256": digest(network),
            "config": cfg,
            "region_buses": mapping,
            "prices_sha256": digest(cfg["prices"]),
            "loads_sha256": digest(cfg["loads"]) if cfg.get("loads") else None,
            "leap_day_policy": "Use network timestamps; never shift subsequent hours",
            "price_filter": "none",
            "missing_data": "paired observations only",
        }
    )
    rows = []

    def compare(kind, key, sim, ref, unit):
        for freq in ["hourly", "daily", "weekly"]:
            try:
                metrics, pair = paired_metrics(
                    sim, ref, freq, cfg.get("min_aggregate_coverage", 1.0)
                )
            except ValueError as exc:
                if str(exc) != "No valid matched observations at requested resolution":
                    raise
                matched = sim.notna() & ref.reindex(sim.index).notna()
                rows.append({
                    "kind": kind, "series": key, "resolution": freq, "unit": unit,
                    "status": "NOT SCORED: insufficient paired observations",
                    "n_expected_hours": len(sim), "n_matched_hours": int(matched.sum()),
                    "coverage": float(matched.mean()), "n_scored": 0,
                })
                continue
            if kind == "flow" and freq == "hourly":
                metrics.update(
                    {
                        "direction_agreement": (
                            np.sign(pair.model) == np.sign(pair.observed)
                        ).mean(),
                        "model_net_twh_matched": pair.model.sum() / 1e6,
                        "observed_net_twh_matched": pair.observed.sum() / 1e6,
                        "model_export_twh_matched": pair.model.clip(lower=0).sum()
                        / 1e6,
                        "observed_export_twh_matched": pair.observed.clip(lower=0).sum()
                        / 1e6,
                        "model_import_twh_matched": -pair.model.clip(upper=0).sum()
                        / 1e6,
                        "observed_import_twh_matched": -pair.observed.clip(
                            upper=0
                        ).sum()
                        / 1e6,
                    }
                )
            rows.append(
                {
                    "kind": kind,
                    "series": key,
                    "resolution": freq,
                    "unit": unit,
                    "status": "scored",
                    **metrics,
                }
            )
            if freq == "hourly":
                stem = f"{kind}_{key}".replace("->", "_to_").replace("/", "_")
                pair.to_csv(output / f"{stem}_paired.csv", index_label="utc_timestamp")
                plot_pair(pair, output / f"{stem}.png", unit)

    for key in model:
        compare("price", key, model[key], obs[key], "EUR/MWh")
    for a, b in cfg.get("spreads", []):
        compare("spread", f"{a}-{b}", model[a] - model[b], obs[a] - obs[b], "EUR/MWh")
    model.to_csv(output / "model_prices_raw.csv", index_label="utc_timestamp")
    shedding = load_shedding(n).reindex(model.index)
    shedding.to_csv(output / "load_shedding_mw.csv", index_label="utc_timestamp")
    report["load_shedding_mwh"] = float(shedding.sum().sum())
    report["load_shedding_hours"] = int((shedding.sum(axis=1) > 1e-6).sum())
    if cfg.get("borders"):
        flows = border_flows(n, mapping, cfg["borders"]).reindex(model.index)
        flows.to_csv(output / "model_border_flows_mw.csv", index_label="utc_timestamp")
        if cfg.get("flows"):
            if cfg.get("flow_kind") not in ["physical", "scheduled_day_ahead"]:
                raise ValueError("Set flow_kind to physical or scheduled_day_ahead")
            observed_flows = read_hourly(cfg["flows"])
            for key in flows:
                compare("flow", key, flows[key], observed_flows[key], "MW")
            report["flow_validation"] = cfg["flow_kind"]
            report["flows_sha256"] = digest(cfg["flows"])
        else:
            report["flow_validation"] = (
                "NOT PERFORMED: no historical flow file configured"
            )
            LOG.warning(report["flow_validation"])
    pd.DataFrame(rows).to_csv(output / "scores.csv", index=False)
    write_json(output / "audit.json", report)
    print(
        f"Price benchmark written to {output}; {report.get('flow_validation', 'no borders requested')}"
    )


def solve(
    network, output, solver, threads, marginal_costs=None, transmission_losses=None
):
    import pypsa

    n = pypsa.Network(network)
    report = audit(n)
    if any(report["extendable"].values()):
        raise ValueError(
            "Hindcast must use fixed historical capacities; extendable assets found"
        )
    if not n.investment_periods.empty:
        raise ValueError(
            "Multi-investment-period networks are not historical dispatch cases"
        )
    for df in [n.generators, n.links]:
        if "committable" in df and df.committable.any():
            raise ValueError(
                "This replication expects continuous dispatch without unit commitment"
            )
    if marginal_costs:
        costs = read_hourly(marginal_costs)
        if set(costs.columns) != set(n.generators.index):
            raise ValueError(
                "Marginal-cost CSV must contain every generator exactly once"
            )
        values = costs.reindex(hourly_network(n))
        if values.isna().any().any():
            raise ValueError(
                "Marginal costs missing for one or more snapshots/generators"
            )
        values.index = n.snapshots
        n.generators_t.marginal_cost = values
    elif report["time_varying_cost_generators"] == 0:
        raise ValueError(
            "No dynamic generator costs: supply verified hourly --marginal-costs"
        )
    if transmission_losses is None:
        transmission_losses = (
            n.meta.get("solving", {}).get("options", {}).get("transmission_losses", 0)
        )
    output = Path(output)
    if output.exists() or output.resolve() == Path(network).resolve():
        raise ValueError(
            "Choose a new output path; reference networks are never overwritten"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    # Preserve the source network's load-shedding components and cost/sign.
    reference_prices_stored = n.buses_t.marginal_price.copy()
    reference_objective = float(n.objective) if n.is_solved else None
    # Fresh optimisation discards stored dispatch/duals as a source of results.
    options = (
        {"threads": threads, "solver": "ipm", "run_crossover": "off"}
        if solver == "highs"
        else {"Threads": threads, "Method": 2, "Crossover": 0}
    )
    status, condition = n.optimize(
        solver_name=solver,
        solver_options=options,
        transmission_losses=transmission_losses,
    )
    if status != "ok" or condition != "optimal":
        raise RuntimeError(
            f"Solve failed: {status}/{condition}; no solved network exported"
        )
    n.meta = dict(
        n.meta,
        hindcast_reproduction={
            "source": str(network),
            "source_sha256": digest(network),
            "pypsa_version": pypsa.__version__,
            "solver": solver,
            "transmission_losses": transmission_losses,
            "marginal_costs_sha256": digest(marginal_costs) if marginal_costs else None,
            "scope": "Fixed-capacity full-horizon re-solve; external callbacks not reconstructed",
        },
    )
    n.export_to_netcdf(output)
    check = audit(n)
    check["objective"] = float(n.objective)
    check["source_objective"] = reference_objective
    if not reference_prices_stored.empty:
        errors = (n.buses_t.marginal_price - reference_prices_stored).abs()
        check["price_mae_against_source"] = errors.mean().to_dict()
        check["price_max_abs_error_against_source"] = errors.max().to_dict()
    write_json(output.with_suffix(".audit.json"), check)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    fetch = sub.add_parser("fetch-reference")
    fetch.add_argument("--year", type=int, choices=range(2020, 2025), default=2023)
    fetch.add_argument("--directory", default="data/hindcast/reference")
    inspect = sub.add_parser("inspect")
    inspect.add_argument("network")
    inspect.add_argument("--output", default="results/hindcast-inspect.json")
    prepare = sub.add_parser("prepare-zonal")
    prepare.add_argument("network")
    prepare.add_argument("--config", default="config/hindcast/benchmark.yaml")
    prepare.add_argument("--flow-sources", default="data/hindcast/flow-sources-zonal-2023.yaml")
    run = sub.add_parser("solve")
    run.add_argument("network")
    run.add_argument("--output", required=True)
    run.add_argument("--solver", choices=["highs", "gurobi"], default="highs")
    run.add_argument("--threads", type=int, default=1)
    run.add_argument("--marginal-costs")
    run.add_argument("--transmission-losses", type=int, default=None)
    bench = sub.add_parser("benchmark")
    bench.add_argument("--config", default="config/hindcast/benchmark.yaml")
    bench.add_argument("--network")
    bench.add_argument("--output")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    if args.command == "fetch-reference":
        fetch_reference(args.year, args.directory)
    elif args.command == "prepare-zonal":
        prepare_zonal(args.network, args.config, args.flow_sources)
    elif args.command == "inspect":
        import pypsa

        write_json(args.output, audit(pypsa.Network(args.network)))
    elif args.command == "solve":
        solve(
            args.network,
            args.output,
            args.solver,
            args.threads,
            args.marginal_costs,
            args.transmission_losses,
        )
    else:
        cfg = yaml.safe_load(Path(args.config).read_text())
        benchmark(args.network or cfg["network"], cfg, args.output or cfg["output"])


if __name__ == "__main__":
    main()
