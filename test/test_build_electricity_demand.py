# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-License-Identifier: MIT

"""Exercise the demand entry point with incomplete synthetic-year coverage."""

import runpy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts import _helpers


class AttrDict(dict):
    __getattr__ = dict.__getitem__


@pytest.mark.parametrize(
    "synthetic_year,observed_gap", [(2023, False), (2024, True), (2023, True)]
)
def test_synthetic_coverage_and_observed_precedence(
    tmp_path, monkeypatch, caplog, synthetic_year, observed_gap
):
    snapshots = pd.date_range("2024-01-01", periods=2, freq="h")
    observed = pd.DataFrame({"DE": [10.0, 20.0], "SE": [30.0, 40.0]}, index=snapshots)
    if observed_gap:
        observed.loc[snapshots[1], "SE"] = np.nan
    inputs = AttrDict()
    for source in ["opsd", "entsoe", "neso"]:
        inputs[source] = str(tmp_path / f"{source}.csv")
        observed.to_csv(inputs[source], date_format="%Y-%m-%dT%H:%M:%SZ")
    synthetic = pd.DataFrame(
        {"DE": [999.0, 999.0], "SE": [80.0, 90.0]},
        index=pd.date_range(f"{synthetic_year}-01-01", periods=2, freq="h"),
    )
    inputs.synthetic = str(tmp_path / "synthetic.csv")
    synthetic.to_csv(inputs.synthetic)
    output = tmp_path / "demand.csv"
    job = SimpleNamespace(
        input=inputs,
        output=[str(output)],
        params=AttrDict(
            snapshots={},
            drop_leap_day=True,
            countries=["DE", "SE"],
            load={
                "fill_gaps": {
                    "enable": False,
                    "interpolate_limit": 6,
                    "time_shift_for_large_gaps": "1w",
                },
                "manual_adjustments": False,
                "supplement_synthetic": True,
                "fixed_year": 2024,
            },
        ),
    )
    monkeypatch.setattr(_helpers, "configure_logging", lambda _: None)
    monkeypatch.setattr(_helpers, "set_scenario_config", lambda _: None)
    monkeypatch.setattr(_helpers, "get_snapshots", lambda *_: snapshots)

    def run():
        runpy.run_path(
            str(
                Path(__file__).resolve().parents[1]
                / "scripts/build_electricity_demand.py"
            ),
            run_name="__main__",
            init_globals={"snakemake": job},
        )

    if observed_gap and synthetic_year != 2024:
        with pytest.raises(AssertionError, match="Load data contains nans"):
            run()
        assert not output.exists()
    else:
        run()
        result = pd.read_csv(output, index_col=0, parse_dates=True)
        expected = observed.copy()
        if observed_gap:
            expected.loc[snapshots[1], "SE"] = 90
        np.testing.assert_allclose(result, expected)
        assert result.index.equals(snapshots)
    if synthetic_year != 2024:
        assert "Synthetic demand does not cover 2 requested hours" in caplog.text
