# SPDX-License-Identifier: MIT
"""Dispatch a fixed historical fleet with and without the treaty Bornholm project."""

import argparse
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

from scripts import transmission_assumptions
from scripts.dispatch_checks import check as check_dispatch
from scripts.run_bornholm_matched import project_legs, require
from scripts.run_bornholm_welfare import digest, write_json

HUB = "BEI_HUB"
WIND = "BEI_WIND"
# Eurostat euro-area annual average HICP. Used only for the illustrative capital-cost comparison.
EUR2025_PER_SOURCE_EUR = {2023: 1.024 * 1.021, 2024: 1.021, 2025: 1.0}
LANDINGS = {"BEI_TO_DK2": "DKE1", "BEI_TO_DE": "DE00"}
TEMPLATE = Path(
    "resources/bornholm_2035_final_v1/central_hybrid/networks/composed_2035.nc"
)


def fixed_fleet(n):
    return {
        component: n.components[component].static[attribute].copy()
        for component, attribute in (
            ("Generator", "p_nom"),
            ("Link", "p_nom"),
            ("Line", "s_nom"),
            ("Transformer", "s_nom"),
            ("Store", "e_nom"),
            ("StorageUnit", "p_nom"),
        )
    }


def validate_historical(n, year):
    expected = pd.date_range(str(year), str(year + 1), freq="h", inclusive="left")
    require(
        isinstance(n.snapshots, pd.DatetimeIndex)
        and n.snapshots.equals(expected)
        and np.allclose(n.snapshot_weightings, 1),
        f"Expected the full hourly {year} hindcast with unit weights",
    )
    require(
        n.investment_periods.empty and not n.has_scenarios,
        "Expected one historical period",
    )
    for component, values in fixed_fleet(n).items():
        static = n.components[component].static
        require(
            not static[f"{values.name}_extendable"].any(),
            f"Historical {component} capacities must already be fixed",
        )
        require(
            not static.get("committable", pd.Series(False, index=static.index)).any(),
            f"Unexpected unit commitment on historical {component}",
        )
    require(
        not n.generators_t.marginal_cost.empty,
        "Historical dynamic fuel costs are missing",
    )
    require(
        not any(
            table.index.astype(str).str.startswith("BEI_").any()
            for table in (n.buses, n.generators, n.links)
        ),
        "Historical source already contains Bornholm project assets",
    )
    audit = transmission_assumptions.validate(n)
    if audit is not None:
        require(
            audit["year"] == year, "Historical transmission limits use another year"
        )


def historical_wind_profile(n):
    g = n.generators
    offshore = g.index[
        g.bus.eq("DKE1") & g.carrier.str.startswith("offwind") & g.p_nom.gt(0)
    ].intersection(n.components["Generator"].active_assets)
    require(
        len(offshore) > 0,
        "No existing DK2 offshore wind fleet for the historical profile",
    )
    availability = n.generators_t.p_max_pu
    profile = pd.DataFrame(
        np.broadcast_to(
            g.loc[offshore, "p_max_pu"].to_numpy(), (len(n.snapshots), len(offshore))
        ).copy(),
        index=n.snapshots,
        columns=offshore,
    )
    dynamic = offshore.intersection(availability.columns)
    require(len(dynamic) > 0, "Historical DK2 offshore wind has no weather series")
    profile.loc[:, dynamic] = availability.loc[n.snapshots, dynamic]
    require(
        np.isfinite(profile).all().all()
        and profile.ge(0).all().all()
        and profile.le(1).all().all(),
        "Invalid historical wind availability",
    )
    capacity = g.loc[offshore, "p_nom"]
    return profile @ (capacity / capacity.sum()), capacity.to_dict()


def add_project(n, template, profile):
    for bus in LANDINGS.values():
        require(
            bus in n.buses.index and n.buses.at[bus, "carrier"] == "AC",
            f"Missing historical AC landing {bus}",
        )
    project_legs(template, {"project": "hybrid"}, 1200)
    hub = template.buses.loc[HUB]
    # Give the offshore bus its own location, outside background corridor audits.
    n.add("Bus", HUB, x=hub.x, y=hub.y, carrier=hub.carrier, location=HUB)
    wind = template.generators.loc[WIND]
    require(
        wind.p_nom_min == wind.p_nom_max == 3000, "Expected the 3 GW treaty template"
    )
    if wind.carrier not in n.carriers.index:
        n.add("Carrier", wind.carrier)
    n.add(
        "Generator",
        WIND,
        bus=HUB,
        carrier=wind.carrier,
        p_nom=3000,
        p_max_pu=profile,
        marginal_cost=wind.marginal_cost,
        capital_cost=wind.capital_cost,
        lifetime=wind.lifetime,
        efficiency=wind.efficiency,
    )
    annual_cost = float(3000 * wind.capital_cost)
    for forward in LANDINGS:
        for name in (forward, forward + "-reversed"):
            link = template.links.loc[name]
            require(
                not link.get("bus2", ""),
                f"Unexpected third terminal {name}",
            )
            for attr in ("p_min_pu", "p_max_pu", "efficiency", "marginal_cost"):
                require(
                    name not in template.links_t[attr],
                    f"Unexpected varying template {attr}: {name}",
                )
            n.add(
                "Link",
                name,
                bus0=link.bus0,
                bus1=link.bus1,
                carrier="DC",
                p_nom=link.p_nom_min,
                p_min_pu=link.p_min_pu,
                p_max_pu=link.p_max_pu,
                efficiency=link.efficiency,
                marginal_cost=link.marginal_cost,
                length=link.length,
                capital_cost=link.capital_cost,
                lifetime=link.lifetime,
            )
            annual_cost += float(link.p_nom_min * link.capital_cost)
    transmission_assumptions.validate(n)
    return annual_cost


def solve(n, directory, case, solver, threads):
    solver_options = (
        {"Threads": threads, "Method": 2, "Crossover": -1, "NumericFocus": 2}
        if solver == "gurobi"
        else {"threads": threads, "solver": "ipm", "run_crossover": "off"}
    )
    status, termination = n.optimize(
        extra_functionality=transmission_assumptions.add_constraints,
        solver_name=solver,
        solver_options=solver_options,
        transmission_losses=n.meta.get("solving", {})
        .get("options", {})
        .get("transmission_losses", 0),
        include_objective_constant=False,
        log_fn=str(directory / f"{case}.log"),
    )
    require(
        (status, termination) == ("ok", "optimal"), f"{case}: {status}/{termination}"
    )
    require(math.isfinite(n.objective), f"{case}: invalid dispatch cost")
    physical_quality = check_dispatch(n)
    transmission_assumptions.check_solution(n)
    w = n.snapshot_weightings.generators
    shedding = n.generators.index[
        n.generators.carrier.eq("load") & n.generators.sign.gt(0)
    ]
    unserved_mwh = (
        n.generators_t.p.reindex(columns=shedding)
        .mul(n.generators.sign[shedding])
        .sum(axis=1)
        @ w
    )

    path = directory / f"{case}.nc"
    temporary = path.with_suffix(".tmp.nc")
    n.export_to_netcdf(temporary)
    temporary.replace(path)
    generation = (
        n.generators_t.p[WIND]
        if WIND in n.generators.index
        else pd.Series(0.0, index=n.snapshots)
    )
    output_mwh = float(generation @ w)
    capture = (
        float((generation * n.buses_t.marginal_price[HUB]) @ w / output_mwh)
        if output_mwh > 0
        else None
    )
    return dict(
        operating_cost_eur=float(n.objective),
        wind_mwh=output_mwh,
        capture_eur_per_mwh=capture,
        physical_quality=physical_quality,
        load_shedding_mwh=float(unserved_mwh),
        network=str(path),
    )


def run(
    source,
    template_path,
    directory,
    year,
    solver,
    threads,
    *,
    operating_to_eur2025=None,
    check=False,
):
    factor = operating_to_eur2025
    if factor is None:
        factor = EUR2025_PER_SOURCE_EUR.get(year)
    require(
        factor is None or (math.isfinite(factor) and factor > 0),
        "Invalid EUR2025 conversion",
    )
    historical = pypsa.Network(source)
    validate_historical(historical, year)
    original_fleet = fixed_fleet(historical)
    profile, proxy = historical_wind_profile(historical)
    with_project = historical.copy()
    annual_cost = add_project(with_project, pypsa.Network(template_path), profile)
    updated_fleet = fixed_fleet(with_project)
    for component, nominal in original_fleet.items():
        if not nominal.empty:
            pd.testing.assert_series_equal(
                updated_fleet[component].loc[nominal.index], nominal
            )
    if check:
        print(
            f"{year}: {len(profile)} hours; DK2 offshore profile capacity factor {profile.mean():.2%}"
        )
        print(f"Treaty layout annualised capital cost: EUR {annual_cost:,.0f}")
        return
    directory.mkdir(parents=True, exist_ok=False)
    report = dict(
        weather_and_fleet_year=year,
        snapshots=len(profile),
        source=str(source),
        source_sha256=digest(source),
        treaty_template=str(template_path),
        treaty_template_sha256=digest(template_path),
        wind_profile_proxy_mw=proxy,
        wind_profile_capacity_factor=float(profile.mean()),
        project_annualised_capital_cost_eur=annual_cost,
        project_cost_basis="Bornholm 2035 Central treaty template, EUR2025",
        solver=solver,
        operating_cost_to_eur2025_factor=factor,
        currency_note="Euro-area HICP approximation for historical operating costs, unless overridden.",
    )
    report["without_project"] = solve(
        historical, directory, "without_project", solver, threads
    )
    del historical
    write_json(directory / "summary.json", report)
    report["with_project"] = solve(
        with_project, directory, "with_project", solver, threads
    )
    savings = (
        report["without_project"]["operating_cost_eur"]
        - report["with_project"]["operating_cost_eur"]
    )
    savings_2025 = None if factor is None else savings * factor
    report.update(
        operating_savings_source_eur=savings,
        operating_savings_eur2025=savings_2025,
        net_savings_after_project_cost_eur2025=None
        if savings_2025 is None
        else savings_2025 - annual_cost,
        break_even_annualised_project_cost_eur2025=savings_2025,
        load_shedding_difference_mwh=(
            report["with_project"]["load_shedding_mwh"]
            - report["without_project"]["load_shedding_mwh"]
        ),
    )
    write_json(directory / "summary.json", report)
    print(f"Operating savings: EUR {savings:,.0f} in the historical cost basis")
    if savings_2025 is not None:
        print(
            f"Net of project costs: EUR {savings_2025 - annual_cost:,.0f} in 2025 prices"
        )
    if any(
        report[case]["load_shedding_mwh"] > 1
        for case in ("with_project", "without_project")
    ):
        print(
            "Load shedding occurred. Inspect both totals before interpreting savings."
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--year", default=2024, type=int)
    parser.add_argument("--template", type=Path, default=TEMPLATE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--solver", choices=("gurobi", "highs"), default="gurobi")
    parser.add_argument(
        "--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))
    )
    parser.add_argument("--operating-to-eur2025", type=float)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    directory = args.output
    if not args.check and directory is None:
        job = os.environ.get("SLURM_JOB_ID")
        if not job:
            parser.error("Provide --output when running outside Slurm")
        directory = Path(f"results/bornholm_historical/{args.year}/{job}")
    run(
        args.source,
        args.template,
        directory,
        args.year,
        args.solver,
        args.threads,
        operating_to_eur2025=args.operating_to_eur2025,
        check=args.check,
    )


if __name__ == "__main__":
    main()
