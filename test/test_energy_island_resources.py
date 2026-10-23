# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Native resource pipeline and project composition integration regressions."""

import runpy
from pathlib import Path
from types import SimpleNamespace

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import pytest
import xarray as xr
import yaml
from shapely.geometry import box
from snakemake.iocontainers import InputFiles, Log, OutputFiles, Params, Wildcards

from scripts._helpers import offshore_connection_cost
from scripts.add_energy_islands import main


def run_script(script, inputs, outputs, params):
    snake = SimpleNamespace(
        input=InputFiles(fromdict=inputs),
        output=OutputFiles(fromdict=outputs),
        params=Params(fromdict=params),
        wildcards=Wildcards(),
        threads=1,
        log=Log(
            fromdict={
                "python": str(Path(next(iter(outputs.values()))).with_suffix(".log"))
            }
        ),
        rule=script,
        config={"run": {}, "atlite": {"show_progress": False}},
    )
    runpy.run_path(
        f"scripts/{script}.py", run_name="__main__", init_globals={"snakemake": snake}
    )


@pytest.fixture
def project_profile(tmp_path):
    """Run real atlite eligibility and native profile code on synthetic weather."""
    cutout = atlite.Cutout(
        str(tmp_path / "weather.nc"),
        module="era5",
        x=slice(14, 15),
        y=slice(54, 55),
        time=pd.date_range("2013-01-01", periods=4, freq="h"),
    )
    shape = (4, *cutout.shape)
    # Deliberately different cell weather: a copied regional proxy is detectable.
    speed = np.broadcast_to(np.linspace(4, 12, cutout.shape[1]), shape).copy()
    speed += np.arange(4)[:, None, None]
    cutout.data["wnd100m"] = (("time", "y", "x"), speed)
    cutout.data["roughness"] = (("time", "y", "x"), np.full(shape, 0.0002))
    cutout.data.attrs["prepared_features"] = ["wind"]
    cutout.data.to_netcdf(cutout.path)
    area = gpd.GeoDataFrame(
        {"name": ["BEI_WIND"]}, geometry=[box(14.1, 54.1, 14.45, 54.45)], crs=4326
    )
    area_path = str(tmp_path / "area.geojson")
    area.to_file(area_path)
    hubs = pd.DataFrame(
        {"hub_id": ["BEI_HUB"], "country": ["DK"], "x": [15.0], "y": [55.0]}
    )
    wind = pd.read_csv("data/energy_islands/bornholm_wind.csv")
    for name, table in [("hubs", hubs), ("wind", wind)]:
        table.to_csv(tmp_path / f"{name}.csv", index=False)
    availability = str(tmp_path / "availability.nc")
    simple = {"offwind-ac": {"natura": False, "excluder_resolution": 200}}
    run_script(
        "determine_availability_matrix",
        {"regions": area_path, "cutout": str(cutout.path)},
        {"nc": availability},
        {
            "technology": "offwind-ac",
            "renewable": simple,
            "project_area": True,
            "plot_availability_matrix": False,
        },
    )
    with open("config/config.default.yaml") as f:
        renewable = yaml.safe_load(f)["renewable"]
    profile = str(tmp_path / "profile.nc")
    run_script(
        "build_renewable_profiles",
        {
            "availability_matrix": availability,
            "cutout": str(cutout.path),
            "resource_regions": area_path,
            "distance_regions": area_path,
            "energy_island_hubs": str(tmp_path / "hubs.csv"),
            "energy_island_wind": str(tmp_path / "wind.csv"),
        },
        {"profile": profile, "class_regions": str(tmp_path / "classes.geojson")},
        {
            "technology": "offwind-ac",
            "renewable": renewable,
            "snapshots": {
                "start": "2013-01-01",
                "end": "2013-01-01 04:00",
                "inclusive": "left",
            },
            "drop_leap_day": False,
        },
    )
    return tmp_path, cutout, area_path, simple, profile


def test_native_exclusion_and_profile(project_profile):
    tmp_path, cutout, area_path, simple, profile = project_profile
    # Native regional availability before and after subtracting the project.
    region_path = str(tmp_path / "region.geojson")
    gpd.GeoDataFrame(
        {"name": ["DK"]}, geometry=[box(14.05, 54.05, 14.95, 54.95)], crs=4326
    ).to_file(region_path)
    arrays = []
    for excluded in [False, True]:
        path = str(tmp_path / f"generic_{excluded}.nc")
        run_script(
            "determine_availability_matrix",
            {
                "regions": region_path,
                "cutout": str(cutout.path),
                "project_exclusion": area_path if excluded else [],
            },
            {"nc": path},
            {
                "technology": "offwind-ac",
                "renewable": simple,
                "plot_availability_matrix": False,
            },
        )
        arrays.append(xr.load_dataarray(path).sel(bus="DK"))
    project = xr.load_dataarray(tmp_path / "availability.nc").sel(bus="BEI_WIND")
    np.testing.assert_allclose(arrays[0] - arrays[1], project, atol=0.005)
    assert float(arrays[1].sum()) < float(arrays[0].sum())
    ds = xr.load_dataset(profile)
    assert ds.profile.sizes["time"] == 4
    assert (
        np.isfinite(ds.profile).all() and ((ds.profile >= 0) & (ds.profile <= 1)).all()
    )
    assert float(ds.average_distance.item()) > 40  # distance to hub, not polygon (zero)
    assert float(ds.p_nom_max.sum()) > 0


def setup_project(project_profile, case):
    tmp_path, cutout, _, _, profile = project_profile
    n = pypsa.Network()
    n.set_snapshots(cutout.data.time.to_index())
    n.add("Carrier", ["AC", "DC", "offwind-ac", "gas"])
    n.add("Bus", ["DK", "DE"], country=["DK", "DE"], carrier="AC")
    links = pd.read_csv("data/energy_islands/bornholm_links.csv")
    links["target_bus"] = links.target_country
    links.to_csv(tmp_path / "links.csv", index=False)
    paths = {
        f"energy_island_{key}": str(tmp_path / f"{key}.csv")
        for key in ["hubs", "wind", "links"]
    }
    paths["energy_island_profile"] = profile
    costs = pd.DataFrame(
        {
            "capital_cost": [100000.0, 10.0, 20.0, 100.0, 50.0, 10.0, 15.0],
            "marginal_cost": 0.0,
            "efficiency": 1.0,
            "lifetime": 30.0,
        },
        index=[
            "offwind",
            "HVDC overhead",
            "HVDC submarine",
            "HVDC inverter pair",
            "offwind-ac-station",
            "offwind-ac-connection-submarine",
            "offwind-ac-connection-underground",
        ],
    )
    config = {
        "enabled": True,
        "case": case,
        "profile_mode": "project",
        "exclude_project_area": True,
        "subtract_generic_potential": False,
    }
    return n, paths, config, costs


@pytest.mark.parametrize(
    "case,capacities", [("radial", [3000]), ("hybrid", [1200, 2000])]
)
@pytest.mark.parametrize(
    "budget_type,budget_column",
    [
        ("transmission_volume_expansion_limit", "length"),
        ("transmission_expansion_cost_limit", "capital_cost"),
    ],
)
def test_native_aggregation_costs_and_budget(
    project_profile, case, capacities, budget_type, budget_column
):
    from scripts.prepare_network import apply_temporal_aggregation

    n, paths, config, costs = setup_project(project_profile, case)
    n.add(
        "Generator",
        "generic",
        bus="DK",
        carrier="offwind-ac",
        p_nom_extendable=True,
        p_nom_max=800,
        capital_cost=1e9,
    )
    main(n, paths, config, costs, 2035, stage="generation")
    hourly = n.generators_t.p_max_pu["BEI_WIND"].copy()
    assert n.links.empty
    assert n.generators.at["generic", "p_nom_max"] == 800
    distance = n.generators.at["BEI_WIND", "collection_distance_km"]
    expected_cost = (
        100000 + 50 + offshore_connection_cost(costs, "offwind-ac", distance)
    )
    assert n.generators.at["BEI_WIND", "capital_cost"] == pytest.approx(expected_cost)
    weights = Path(paths["energy_island_profile"]).with_name("weights.csv")
    n.snapshot_weightings.resample("2h").sum().to_csv(weights)
    apply_temporal_aggregation(
        n,
        SimpleNamespace(snapshot_weightings=weights),
        SimpleNamespace(
            clustering_temporal={
                "averaging": "2h",
                "segmentation": False,
                "representative": False,
            }
        ),
    )
    pd.testing.assert_series_equal(
        n.generators_t.p_max_pu["BEI_WIND"],
        hourly.resample("2h").mean(),
        check_names=False,
        check_freq=False,
    )
    n.add(
        "GlobalConstraint",
        "grid_limit",
        type=budget_type,
        carrier_attribute="AC, DC",
        sense="<=",
        constant=0,
    )
    main(n, paths, config, costs, 2035, stage="links")
    assert list(n.links.p_nom_min) == capacities
    assert list(n.links.p_nom_max) == capacities
    assert n.global_constraints.at["grid_limit", "constant"] == pytest.approx(
        float(n.links.p_nom_min @ n.links[budget_column])
    )
    n.add("Load", "load", bus="DE", p_set=1000)
    n.add("Generator", "gas", bus="DE", carrier="gas", p_nom=2000, marginal_cost=100)
    # Match the native compose -> solve boundary. Temporal aggregation replaces
    # network internals; the solve rule loads a fresh Network from the export.
    composed = Path(paths["energy_island_profile"]).with_name(f"{case}_composed.nc")
    n.export_to_netcdf(composed)
    n = pypsa.Network(composed)
    status, condition = n.optimize(solver_name="highs", threads=1)
    assert (status, condition) == ("ok", "optimal")
    np.testing.assert_allclose(n.links.p_nom_opt, capacities)
    capex = (
        n.generators.p_nom_opt @ n.generators.capital_cost
        + n.links.p_nom_opt @ n.links.capital_cost
    )
    opex = n.snapshot_weightings.objective @ n.generators_t.p.mul(
        n.generators.marginal_cost
    ).sum(axis=1)
    assert n.objective == pytest.approx(capex + opex)
    n.export_to_netcdf(
        Path(paths["energy_island_profile"]).with_name(f"{case}_solved.nc")
    )


def test_no_double_subtraction(project_profile):
    n, paths, config, costs = setup_project(project_profile, "hybrid")
    config["subtract_generic_potential"] = True
    with pytest.raises(ValueError, match="cannot both"):
        main(n, paths, config, costs, 2035)
    assert "BEI_HUB" not in n.buses.index


def test_scenario_resource_settings():
    for path in ["config/config.bornholm.yaml", "config/config.bornholm.smoke.yaml"]:
        cfg = yaml.safe_load(Path(path).read_text())
        scenarios = yaml.safe_load(Path(cfg["run"]["scenarios"]["file"]).read_text())
        for scenario in scenarios.values():
            island = {**cfg["energy_islands"], **scenario["energy_islands"]}
            assert (
                island["exclude_project_area"]
                and not island["subtract_generic_potential"]
            )
            assert Path(island["project_area"]).is_file()
            assert Path(island["wind"]).is_file()
            assert Path(island["links"]).is_file()
