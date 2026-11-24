# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-License-Identifier: MIT

"""Demand allocation with landless grid buses and mixed native data sources."""

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

from scripts.build_electricity_demand_base import upsample_load


@pytest.fixture
def demand_inputs(tmp_path):
    names = ["NL_land", "NL_hv", "NO_land", "GB_land", "GBNI_land"]
    regions = gpd.GeoDataFrame(
        {"name": names, "country": ["NL", "NL", "NO", "GB", "GB"]},
        geometry=[box(2 * i, 0, 2 * i + 1, 1) for i in range(5)],
        crs=3035,
    )
    regions_fn = tmp_path / "regions.geojson"
    regions.to_file(regions_fn)
    n = pypsa.Network()
    n.add("Bus", names, country=regions.country.values, substation_lv=True)
    n.buses.loc["NL_hv", "substation_lv"] = False
    n.add("Bus", "NL_no_cell", country="NL", substation_lv=True)
    load = pd.DataFrame(
        {"NL": [10.0, 20.0], "NO": [30.0, 40.0], "GB": [50.0, 60.0]},
        index=pd.date_range("2023-01-01", periods=2, freq="h"),
    )
    load_fn = tmp_path / "load.csv"
    load.to_csv(load_fn)
    raster_fn = tmp_path / "energy_atlas.tif"
    with rasterio.open(
        raster_fn,
        "w",
        driver="GTiff",
        height=1,
        width=10,
        count=1,
        dtype="float32",
        crs="EPSG:3035",
        transform=from_origin(0, 1, 1, 1),
        nodata=-1,
    ) as raster:
        raster.write(
            np.array([[2, -1, 4, -1, 6, -1, 8, -1, 10, -1]], dtype="float32"), 1
        )
    nuts = regions.assign(pop=1.0, gdp=1.0)
    nuts_fn = tmp_path / "nuts.geojson"
    nuts.to_file(nuts_fn)
    # GB statistics cover Great Britain; NI uses the native NUTS fallback.
    gb_fn = tmp_path / "gb.geojson"
    gpd.GeoDataFrame(
        {"LAD24CD": ["GBcode"]}, geometry=[regions.geometry.iloc[3]], crs=3035
    ).to_file(gb_fn)
    excel_fn = tmp_path / "gb.xlsx"
    with pd.ExcelWriter(excel_fn) as writer:
        pd.DataFrame(
            {
                "Local authority": ["GB district"],
                "Code": ["GBcode"],
                "Total consumption\n(GWh):\nAll meters": [1.0],
            }
        ).to_excel(writer, sheet_name="2019", startrow=4, index=False)
    kwargs = dict(
        n=n,
        regions_fn=regions_fn,
        load_fn=load_fn,
        raster_fn=raster_fn,
        gb_excel_fn=excel_fn,
        gb_geojson_fn=gb_fn,
        nuts3_fn=nuts_fn,
        distribution_key={"gdp": 0.6, "population": 0.4},
    )
    return kwargs, regions, load


@pytest.mark.parametrize("substation_only", [False, True])
def test_allocation_preserves_hourly_country_totals(demand_inputs, substation_only):
    kwargs, _, load = demand_inputs
    n = kwargs["n"]
    original_buses = n.buses.copy(deep=True)
    result = upsample_load(**kwargs, substation_only=substation_only)
    assert np.isfinite(result).all()
    assert "NL_no_cell" not in result.bus
    assert ("NL_hv" in result.bus) != substation_only
    for country in load:
        buses = result.bus.to_index().intersection(
            n.buses.index[n.buses.country.eq(country)]
        )
        np.testing.assert_allclose(result.sel(bus=buses).sum("bus"), load[country])
    np.testing.assert_allclose(result.sel(bus="GB_land"), load.GB * 2 / 3)
    np.testing.assert_allclose(result.sel(bus="GBNI_land"), load.GB / 3)
    pd.testing.assert_frame_equal(n.buses, original_buses)


def test_missing_input_geometry_identifies_the_bus(demand_inputs):
    kwargs, regions, _ = demand_inputs
    regions.loc[0, "geometry"] = None
    regions.to_file(kwargs["regions_fn"])
    with pytest.raises(ValueError, match="Missing or empty land geometries.*NL_land"):
        upsample_load(**kwargs)


def test_country_without_eligible_land_is_not_silently_dropped(demand_inputs):
    kwargs, regions, _ = demand_inputs
    regions.loc[regions.country.ne("NO")].to_file(kwargs["regions_fn"])
    with pytest.raises(ValueError, match="No eligible land cells.*NO"):
        upsample_load(**kwargs)


def test_zero_distribution_data_is_not_replaced_by_arbitrary_weights(demand_inputs):
    kwargs, _, _ = demand_inputs
    nuts = gpd.read_file(kwargs["nuts3_fn"])
    nuts.loc[nuts.country.eq("NO"), ["pop", "gdp"]] = 0.0
    nuts.to_file(kwargs["nuts3_fn"])
    with pytest.raises(
        ValueError, match="Invalid or missing demand distribution weights for NO"
    ):
        upsample_load(**kwargs)
