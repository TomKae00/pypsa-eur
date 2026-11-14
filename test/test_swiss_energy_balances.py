"""Preserve native energy categories across the original and revised BFE layouts."""

import pandas as pd
import pytest

from scripts.build_swiss_energy_balances import (
    build_energy_balances,
    residential_cooking,
    select_row,
)


def write_workbook(path, *, revised):
    """Minimal real Excel workbook with known units and both source layouts."""
    sheets = {
        "Tabelle1": {"int. Flugverkehr": 36.0},
        "Tabelle17": {
            "Total" if revised else "Total Endenergieverbrauch": 72.0,
            "Raumwärme": 36.0,
            "Warmwasser": 18.0,
            **(
                {"Kochen": 10.8, "Geschirrspüler": 3.6}
                if revised
                else {"Kochen / Geschirrspüler": 14.4}
            ),
        },
        "Tabelle18": {
            "Total": 36.0,
            "Raumwärme": 7.2,
            "Warmwasser": 3.6,
            "Kochherde": 1.8,
        },
        "Tabelle26": {
            "Total" if revised else "Total Endenergie": 72.0,
            "Raumwärme": 36.0,
            "Warmwasser": 18.0,
            "Prozesswärme": 7.2,
        },
        "Tabelle28": {
            "Total" if revised else "Total Elektrizität": 36.0,
            "Raumwärme": 7.2,
            "Warmwasser": 3.6,
            "Prozesswärme": 1.8,
        },
        "Tabelle35": {
            "Schiene": 7.2,
            "Strasse": 36.0,
            "Luft (Inland)": 3.6,
            "Wasser": 1.8,
        },
        "Tabelle38": {"Schiene*" if revised else "Schiene": 7.2, "Strasse": 3.6},
    }
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet, rows in sheets.items():
            index_name = (
                "Verkehrsträger"
                if sheet in ("Tabelle35", "Tabelle38")
                else "Verwendungszweck"
            )
            frame = pd.DataFrame({year: pd.Series(rows) for year in (2000, 2023, 2024)})
            frame.to_excel(writer, sheet_name=sheet, startrow=5, index_label=index_name)


def test_excel_layouts_produce_identical_native_categories_and_twh(tmp_path):
    old, revised = tmp_path / "original.xlsx", tmp_path / "revised.xlsx"
    write_workbook(old, revised=False)
    write_workbook(revised, revised=True)
    expected = build_energy_balances(old)
    actual = build_energy_balances(revised)
    pd.testing.assert_frame_equal(actual, expected)
    assert actual.shape == (3, 24)
    assert actual.index.names == ["country", "year"]
    assert actual.loc[("CH", 2023), "total residential"] == pytest.approx(20)
    assert actual.loc[("CH", 2023), "total residential cooking"] == pytest.approx(4)
    assert actual.loc[("CH", 2023), "electricity rail"] == pytest.approx(2)
    assert actual["total international navigation"].eq(0).all()


def test_unrecognised_or_ambiguous_total_labels_fail_explicitly():
    for labels in (["Unrecognised total"], ["Total", "Total Endenergie"]):
        table = pd.DataFrame({2023: [1.0] * len(labels)}, index=labels)
        with pytest.raises(ValueError, match="Expected exactly one"):
            select_row(table, "Total Endenergie", "Total")


def test_split_cooking_requires_both_components_and_preserves_missing_values():
    table = pd.DataFrame(
        {2023: [10.8, float("nan")]}, index=["Kochen", "Geschirrspüler"]
    )
    assert pd.isna(residential_cooking(table).loc[2023])
    with pytest.raises(KeyError):
        residential_cooking(table.drop(index="Geschirrspüler"))
