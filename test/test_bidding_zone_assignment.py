# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-License-Identifier: MIT
"""Zone boundaries must govern both plant capacities and transmission endpoints."""

from types import SimpleNamespace

import geopandas as gpd
import pypsa
import pytest
from shapely.geometry import box

from scripts.base_network import build_admin_shapes
from scripts.build_powerplants import map_to_country_bus
from scripts.cluster_network import busmap_for_admin_regions, clustering_for_n_clusters


def test_plants_and_corridors_follow_zone_boundary_inside_county(tmp_path):
    # A county crossing the zone border was previously assigned wholly to SE04.
    nuts = gpd.GeoDataFrame(
        {"id": ["SEcounty"], "country": ["SE"], "bidding_zone": ["SE04"]},
        geometry=[box(12, 56, 14, 58)],
        crs=4326,
    )
    zones = gpd.GeoDataFrame(
        {"zone_name": ["SE03", "SE04"], "country": ["SE", "SE"]},
        geometry=[box(12, 57, 14, 58), box(12, 56, 14, 57)],
        crs=4326,
    )
    nuts_path, zones_path = tmp_path / "nuts.geojson", tmp_path / "zones.geojson"
    nuts.to_file(nuts_path)
    zones.to_file(zones_path)
    n = pypsa.Network()
    n.add(
        "Bus",
        ["north_a", "north_b", "south"],
        country="SE",
        carrier="AC",
        v_nom=380,
        x=[12.5, 13.5, 13],
        y=[57.5, 57.5, 56.5],
    )
    n.add(
        "Line",
        ["internal", "border_a", "border_b"],
        bus0=["north_a", "north_a", "north_b"],
        bus1=["north_b", "south", "south"],
        s_nom=[9000, 1000, 2000],
        s_max_pu=0.7,
        r=0.01,
        x=0.1,
        length=100,
    )
    n.add(
        "Link",
        ["internal_dc", "border_dc"],
        bus0=["north_a", "north_b"],
        bus1=["north_b", "south"],
        carrier="DC",
        p_nom=[500, 600],
        p_min_pu=-1,
        length=100,
    )
    shapes = build_admin_shapes(
        n,
        nuts_path,
        "administrative",
        {"level": "bz"},
        ["SE"],
        bidding_zones=zones_path,
    )
    admin_path = tmp_path / "admin.geojson"
    shapes[["country", "geometry"]].to_file(admin_path)
    busmap = busmap_for_admin_regions(
        n,
        admin_path,
        SimpleNamespace(countries=["SE"], administrative={"level": "bz"}),
    )
    assert busmap.to_dict() == {"north_a": "SE03", "north_b": "SE03", "south": "SE04"}

    regions = shapes[["geometry"]].rename_axis("name")
    plants = gpd.GeoDataFrame(
        {"Country": "SE", "Capacity": [1000.0, 2000.0, 3000.0]},
        geometry=gpd.points_from_xy(n.buses.x, n.buses.y),
        crs=4326,
        index=n.buses.index.rename("plant"),
    )
    assigned = map_to_country_bus(plants, regions)
    assert assigned.bus.to_dict() == busmap.to_dict()
    assert assigned.groupby("bus").Capacity.sum().to_dict() == {
        "SE03": 3000,
        "SE04": 3000,
    }

    clustered = clustering_for_n_clusters(n, busmap)
    nc = clustered.n
    # Internal assets must not inflate the interzonal capacity.
    assert "internal" not in clustered.linemap.index
    assert nc.lines.s_nom.sum() == pytest.approx(3000)
    assert (nc.lines.s_nom * nc.lines.s_max_pu).sum() == pytest.approx(2100)
    assert set(nc.links.index) == {"border_dc"}
    assert nc.links.at["border_dc", "p_nom"] == 600
    assert set(nc.links.loc["border_dc", ["bus0", "bus1"]]) == {"SE03", "SE04"}
    assert nc.lines.bus0.ne(nc.lines.bus1).all()
