# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-License-Identifier: MIT
"""Retain accepted SciGRID values when selectively correcting gas pipes."""

import pandas as pd
import pytest
from pypsa.geo import haversine_pts
from shapely.geometry import LineString

from scripts.build_gas_network import prepare_dataset
from scripts.cluster_gas_network import aggregate_parallel_pipes, reindex_pipes


def source_pipe(capacity_mw=5000.0, length_km=100.0, **overrides):
    row = {
        "name": "source-pipe",
        "geometry": LineString([(0.0, 0.0), (1.0, 0.0)]),
        "diameter_mm": 600.0,
        "diameter_mm_method": "raw",
        "is_H_gas": 1,
        "is_bothDirection": 0,
        "length_km": length_km,
        "max_cap_M_m3_per_d": capacity_mw / 437.5,
        "max_pressure_bar": 70.0,
        "start_year": 2000,
    }
    return pd.DataFrame([{**row, **overrides}])


@pytest.mark.parametrize("capacity", [1000.0, 5000.0, 8000.0])
def test_preserves_plausible_source_capacity(capacity):
    actual = prepare_dataset(source_pipe(capacity))
    assert actual.p_nom.iloc[0] == pytest.approx(capacity)
    assert actual.length.iloc[0] == 100.0


@pytest.mark.parametrize("capacity", [0.0, 100.0, 100000.0])
def test_still_replaces_capacities_selected_for_correction(capacity):
    actual = prepare_dataset(source_pipe(capacity))
    assert actual.p_nom.iloc[0] == pytest.approx(5000.0)


def test_preserves_nordstream_pressure_exception():
    actual = prepare_dataset(
        source_pipe(100000.0, max_pressure_bar=250.0)
    )
    assert actual.p_nom.iloc[0] == pytest.approx(100000.0)


@pytest.mark.parametrize("length", [1.0, 20.0, 1000.0])
def test_still_replaces_lengths_selected_for_correction(length):
    actual = prepare_dataset(source_pipe(length_km=length))
    expected = 1.5 * haversine_pts([0.0, 0.0], [1.0, 0.0])
    assert actual.length.iloc[0] == pytest.approx(expected)


def test_hungary_serbia_source_capacity_survives_clustering():
    # SciGRID IGGIELGN 1.1.2, DOI 10.5281/zenodo.4767098.
    # This exact segment appears in the uploaded HU00 -> RS00 Link tags.
    source = source_pipe(
        name="PLM_1711_Seg_0_Seg_0_Seg_1",
        geometry=LineString([(19.90758, 46.22199), (19.96264, 46.01847)]),
        diameter_mm=900.0,
        diameter_mm_method="Median(diameter_mm)",
        length_km=23.838341,
        max_cap_M_m3_per_d=12.492560438325587,
        start_year=None,
    )
    actual = prepare_dataset(source)
    assert actual.p_nom.iloc[0] == pytest.approx(5465.495191767444)
    assert actual.length.iloc[0] == pytest.approx(23.838341)
    # Preserve the native conservative estimate for inferred diameter.
    assert actual.p_nom_diameter.iloc[0] == pytest.approx(1500.0)

    actual = actual.drop(columns=["geometry", "point0", "point1", "H_gas"])
    actual["bus0"] = "HU00"
    actual["bus1"] = "RS00"
    reindex_pipes(actual)
    clustered = aggregate_parallel_pipes(actual)
    assert clustered.loc["gas pipeline HU00 -> RS00", "p_nom"] == pytest.approx(
        5465.495191767444
    )


def test_parallel_accepted_capacities_are_conserved():
    source = pd.concat([source_pipe(2000.0), source_pipe(3000.0)], ignore_index=True)
    actual = prepare_dataset(source)
    actual = actual.drop(columns=["geometry", "point0", "point1", "H_gas"])
    actual["bus0"] = "HU00"
    actual["bus1"] = "RS00"
    reindex_pipes(actual)
    assert aggregate_parallel_pipes(actual).p_nom.iloc[0] == pytest.approx(5000.0)
