# SPDX-License-Identifier: MIT
"""Planned-grid selection, project mapping and independent investment checks."""

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import pytest
from shapely.geometry import LineString, box

from scripts.add_energy_islands import main as add_energy_islands
from scripts.add_transmission_projects_and_dlr import select_grid_year
from scripts.build_transmission_projects import connect_new_lines, get_project_files
from scripts.prepare_bornholm_grid import STATUSES, select_projects
from scripts.prepare_network import set_transmission_limit
from scripts.prepare_sector_network import lossy_bidirectional_links
from scripts.solve_network import add_lossy_bidirectional_link_constraints
from test.test_energy_island_resources import project_profile, setup_project  # noqa: F401


def test_committed_2035_subset_is_valid_and_excludes_bornholm():
    root = Path("data/transmission_projects/tyndp2026")
    manifest = json.loads((root / "audit/manifest.json").read_text())
    tables = get_project_files(root)
    assert {k: len(v) for k, v in tables.items()} == {
        k: v for k, v in manifest["selected_counts"].items() if v
    }
    for frame in tables.values():
        assert frame.build_year.le(2035).all()
        assert set(frame.project_status) <= STATUSES
        assert not frame.index.str.startswith("TYNDP2026_1106_").any()
    assert tables["new_links"].p_nom.gt(0).all()
    harmony = tables["new_lines"].loc["TYNDP2026_170_1034"]
    assert harmony.v_nom == 220
    assert harmony.num_parallel == 2
    assert harmony.build_year == 2030
    assert 100 < harmony.length < 200
    assert "TYNDP2026_170_1034" not in tables["new_links"].index
    for project in manifest["technical_review_exclusions"]:
        assert project not in tables["new_links"].index


def test_selection_removes_future_speculative_and_project_rows():
    frame = (
        pd.read_csv("data/transmission_projects/tyndp2026/new_links.csv", index_col=0)
        .iloc[:4]
        .copy()
    )
    frame.index = ["future", "speculative", "TYNDP2026_1106_1793", "accepted"]
    frame.build_year = [2036, 2030, 2030, 2035]
    frame.project_status = [
        "under_construction",
        "under_consideration",
        "under_construction",
        "in_permitting",
    ]
    selected, reasons = select_projects(frame)
    assert selected.index.tolist() == ["accepted"]
    assert reasons.at["TYNDP2026_1106_1793"] == "bornholm_modelled_separately"
    frame.at["accepted", "p_nom"] = np.nan
    with pytest.raises(ValueError, match="positive finite rating"):
        select_projects(frame)


def test_offshore_endpoints_get_distinct_buses():
    n = pypsa.Network()
    n.add("Bus", "mainland", x=20, y=1, carrier="AC")
    lines = pd.DataFrame({"geometry": [LineString([(5, 1), (7, 1)])]}, index=["cable"])
    offshore = gpd.GeoDataFrame(
        {"country": ["DK"]}, geometry=[box(0, 0, 10, 10)], crs=4326
    )
    connected, buses = connect_new_lines(
        lines, n, pd.DataFrame(), offshore, distance_upper_bound=0.4
    )
    assert len(buses) == 2
    assert connected.at["cable", "bus0"] != connected.at["cable", "bus1"]
    assert set(buses.country) == {"DK"}


def test_project_collapsing_to_one_bus_is_not_a_line():
    n = pypsa.Network()
    n.add("Bus", "one", x=0, y=0, carrier="AC")
    lines = pd.DataFrame(
        {"geometry": [LineString([(0, 0), (0.01, 0.01)])]}, index=["loop"]
    )
    connected, _ = connect_new_lines(lines, n, pd.DataFrame())
    assert connected.empty


def test_grid_year_replaces_upgraded_line_and_excludes_future_assets():
    n = pypsa.Network()
    n.add("Bus", ["a", "b"], carrier="AC")
    n.add(
        "Line",
        ["old", "replacement", "future"],
        bus0="a",
        bus1="b",
        carrier="AC",
        s_nom=[100, 200, 500],
        build_year=[1990, 2030, 2040],
        lifetime=[40, np.inf, np.inf],
    )
    n.add(
        "Link", "hydrogen", bus0="a", bus1="b", carrier="H2 pipeline", build_year=2040
    )
    select_grid_year(n, 2035)
    assert n.lines.index.tolist() == ["replacement"]
    assert n.lines.s_nom.sum() == 200
    assert "hydrogen" in n.links.index


@pytest.mark.parametrize("case", ["radial_de", "optimised", "open"])
def test_fixed_background_leaves_project_investment_free(
    project_profile,  # noqa: F811
    case,
    tmp_path,
):
    n, inputs, config, costs = setup_project(project_profile, case)
    # Annual cable costs must be compared with an annual dispatch benefit.
    n.snapshot_weightings.loc[:, :] = 8760 / len(n.snapshots)
    n.add("Bus", ["SE", "PL"], carrier="AC", country=["SE", "PL"])
    costs.loc["HVAC overhead"] = [10, 0, 1, 30]
    n.add(
        "Line",
        "background_ac",
        bus0="DK",
        bus1="DE",
        s_nom=100,
        x=1,
        length=50,
        s_nom_extendable=True,
    )
    n.add(
        "Link",
        "background_dc",
        bus0="DK",
        bus1="DE",
        carrier="DC",
        p_nom=100,
        p_min_pu=-1,
        length=100,
        underwater_fraction=1,
        p_nom_extendable=True,
    )
    n.add(
        "Link",
        "hydrogen",
        bus0="DK",
        bus1="DE",
        carrier="H2 pipeline",
        p_nom_extendable=True,
    )
    set_transmission_limit(n, "v", "fixed", costs)
    assert not n.lines.s_nom_extendable.any()
    assert not n.links.at["background_dc", "p_nom_extendable"]
    assert n.links.at["hydrogen", "p_nom_extendable"]
    n.remove("Link", "hydrogen")
    add_energy_islands(n, inputs, config, costs, 2035)
    assert n.global_constraints.empty
    candidates = n.links.index[n.links.bus0.eq("BEI_HUB")]
    if case != "radial_de":
        assert n.links.loc[candidates, "p_nom_min"].eq(0).all()
    lossy_bidirectional_links(n, "DC", {"efficiency_static": 0.98})
    n.add("Load", "DE demand", bus="DE", p_set=1000)
    n.add("Generator", "gas", bus="DE", carrier="gas", p_nom=2000, marginal_cost=100)
    path = tmp_path / "composed.nc"
    n.export_to_netcdf(path)
    n = pypsa.Network(path)
    assert n.optimize(
        solver_name="highs",
        threads=1,
        extra_functionality=lambda n, snapshots: (
            add_lossy_bidirectional_link_constraints(n)
        ),
    ) == ("ok", "optimal")
    assert n.lines.at["background_ac", "s_nom_opt"] == 100
    assert n.links.at["background_dc", "p_nom_opt"] == 100
    assert n.links.loc[candidates, "p_nom_opt"].sum() > 0
    for name in candidates:
        assert n.links.at[name, "p_nom_opt"] == pytest.approx(
            n.links.at[name + "-reversed", "p_nom_opt"]
        )
    if case == "open":
        # There is no load/benefit in SE or PL in this controlled example.
        assert (
            n.links.loc[["BEI_TO_SE4", "BEI_TO_PL"], "p_nom_opt"].abs().lt(1e-6).all()
        )
