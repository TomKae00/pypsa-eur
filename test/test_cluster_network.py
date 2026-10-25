# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for helper utilities in scripts/cluster_network.py."""

import geopandas as gpd
import pandas as pd
import pypsa
from shapely.geometry import Polygon

from scripts.cluster_network import (
    busmap_from_shapes,
    copperplate_buses,
    sanitize_busmap,
)


def test_sanitize_busmap_strips_values():
    """Whitespace in busmap values is removed while indices stay untouched."""

    raw = pd.Series([" zone_a ", "zone_b  "], index=[" bus1 ", "bus2  "], name="busmap")

    cleaned = sanitize_busmap(raw)

    assert cleaned.index.tolist() == [" bus1 ", "bus2  "]
    assert cleaned.tolist() == ["zone_a", "zone_b"]
    assert cleaned.name == "busmap"


def test_busmap_from_shapes_returns_trimmed_labels():
    """Busmap inferred from shapes produces sanitized labels."""

    n = pypsa.Network()
    n.add("Bus", "bus1", x=0.0, y=0.0, carrier="AC", country="AA")
    n.add("Bus", "bus2", x=2.0, y=0.0, carrier="AC", country="AA")

    shapes = gpd.GeoDataFrame(
        {
            "name": [" region_1 ", "region_2  "],
            "geometry": [
                Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)]),
                Polygon([(1, -1), (3, -1), (3, 1), (1, 1)]),
            ],
        },
        crs="EPSG:4326",
    )

    busmap = busmap_from_shapes(n, shapes)

    assert busmap.index.tolist() == ["bus1", "bus2"]
    assert busmap.tolist() == ["region_1", "region_2"]


def test_copperplate_export_preserves_construction_flags(tmp_path):
    """Market links must coexist with native boolean branch metadata in NetCDF."""
    n = pypsa.Network()
    n.add("Carrier", ["AC", "DC"])
    n.add("Bus", ["DE00", "LUG1", "PL00"], carrier="AC", country=["DE", "LU", "PL"])
    n.add(
        "Link",
        ["existing", "planned"],
        bus0=["DE00", "LUG1"],
        bus1="PL00",
        carrier="DC",
        p_nom=[1000, 2000],
        under_construction=[False, True],
        underwater_fraction=[0.5, 1.0],
    )
    original = n.links.copy()
    copperplate_buses(n, [["DE00", "LUG1"]])

    path = tmp_path / "clustered.nc"
    n.export_to_netcdf(path)
    restored = pypsa.Network(path)
    pd.testing.assert_frame_equal(
        restored.links.loc[
            original.index, ["bus0", "bus1", "p_nom", "underwater_fraction"]
        ],
        original[["bus0", "bus1", "p_nom", "underwater_fraction"]],
    )
    assert restored.links.loc[original.index, "under_construction"].tolist() == [
        False,
        True,
    ]
    copper = restored.links.loc["copper_DE00_LUG1_0"]
    assert not copper.under_construction
    assert copper.carrier == "copper"
    assert copper.p_nom == float("inf")
    assert copper.efficiency == 1.0
    assert copper.p_min_pu == -1.0
    assert restored.buses.country.to_dict() == {
        "DE00": "DE",
        "LUG1": "LU",
        "PL00": "PL",
    }
