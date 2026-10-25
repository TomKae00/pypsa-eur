# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Regression for Instrat CSVs whose date column is not the first column."""

import runpy
from types import SimpleNamespace

import pandas as pd
import pytest

from scripts import _helpers


@pytest.mark.parametrize(
    "columns",
    [
        ["date", "volume", "indeks", "price"],
        ["volume", "date", "indeks", "price"],
    ],
)
def test_hourly_co2_prices_ignore_source_column_order(tmp_path, monkeypatch, columns):
    source = tmp_path / "prices_eu_ets_all.csv"
    output = tmp_path / "co2_price.csv"
    pd.DataFrame(
        {
            "volume": [300000, 400000, 500000],
            "date": [
                "2013-01-01T11:02:00Z",
                "2013-01-02T11:01:00Z",
                "2013-01-03T11:00:00Z",
            ],
            "indeks": ["t2pa"] * 3,
            "price": [10.0, 20.0, 30.0],
        }
    )[columns].to_csv(source, index=False)
    monkeypatch.setattr(_helpers, "configure_logging", lambda _: None)
    snake = SimpleNamespace(
        input={"csv": source},
        output={"csv": output},
        params={"rolling_window": 3},
        config={"run": {}},
    )
    runpy.run_path(
        "scripts/build_co2_prices.py",
        run_name="__main__",
        init_globals={"snakemake": snake},
    )
    actual = pd.read_csv(output, index_col="date", parse_dates=True)["price"]
    expected = pd.Series(
        [15.0] * 24 + [20.0] * 24 + [25.0],
        index=pd.date_range("2013-01-01 11:00", periods=49, freq="h", name="date"),
        name="price",
    )
    pd.testing.assert_series_equal(actual, expected, check_freq=False)
