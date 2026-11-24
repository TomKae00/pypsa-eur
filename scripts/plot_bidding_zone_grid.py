# SPDX-License-Identifier: MIT
"""Compare the 2025 and 2035 electricity grid from prepared PyPSA networks."""

import argparse
import json
import logging
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
from pyproj import Transformer
from shapely import make_valid, transform, union_all
from shapely.geometry import box, shape

from scripts import transmission_assumptions as transmission

LAND_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    "v5.1.2/geojson/ne_50m_admin_0_countries.geojson"
)
LABELS = {
    "GB00": "GB",
    "GBZE": "Shetland",
    "DE-LU": "DE/LU",
    "IE/NI": "SEM",
    "ESMA": "Balearics",
    "FR15": "Corsica",
    "ITN1": "IT-N",
    "ITCN": "IT-CN",
    "ITCS": "IT-CS",
    "ITS1": "IT-S",
    "ITCA": "IT-CAL",
    "ITSA": "IT-SAR",
    "ITSI": "IT-SIC",
    "BEI_HUB": "Bornholm",
}
LABEL_OFFSETS = {
    "GBZE": (5, 5),
    "GB00": (4, -9),
    "IE/NI": (-6, -8),
    "FR15": (7, 0),
    "ESMA": (1, -14),
    "ITCA": (-5, 7),
    "ITSI": (-7, -9),
    "ITS1": (5, 0),
    "ITN1": (-6, 9),
    "NL00": (-7, 5),
    "BE00": (-8, -5),
    "ME00": (-8, 2),
    "XK00": (-7, -8),
    "AL00": (-7, 1),
    "MK00": (5, -8),
    "SI00": (6, -6),
    "HR00": (-8, 4),
    "NO2": (-9, -8),
    "NO1": (6, -7),
    "NO5": (-9, 1),
    "DK1": (-9, -10),
    "DK2": (5, -9),
    "SE4": (5, 4),
    "BEI_HUB": (8, -6),
}
# Representative positions for schematic connectors. The zone polygons show geography.
ANCHORS = {"DE-LU": (10.15, 51.3), "IE/NI": (-7.3, 53.5), "ITCA": (16.52, 38.93)}
BEND = {
    ("GB00", "NO2"): 0.16,
    ("DK1", "GB00"): 0.08,
    ("GB00", "NL00"): -0.13,
    ("BE00", "GB00"): -0.11,
    ("FR00", "GB00"): 0.13,
    ("DE-LU", "NO2"): 0.12,
    ("DK1", "NL00"): 0.20,
    ("NL00", "NO2"): -0.04,
    ("DE-LU", "DK2"): -0.12,
    ("DE-LU", "SE4"): -0.19,
    ("PL00", "SE4"): 0.09,
    ("LT", "SE4"): -0.08,
    ("FI", "SE3"): 0.03,
    ("BEI_HUB", "DE-LU"): -0.10,
    ("BEI_HUB", "DK2"): 0.08,
    ("DK1", "NO2"): 0.13,
    ("DK1", "SE3"): 0.16,
    ("NO2", "NO5"): -0.16,
    ("EE", "FI"): -0.10,
}
NUMBER_OFFSETS = {
    ("SE2", "SE3"): (21, -8),
    ("SE3", "SE4"): (30, -5),
    ("DK2", "SE4"): (-2, 20),
    ("BEI_HUB", "DK2"): (12, -10),
    ("BEI_HUB", "DE-LU"): (17, -3),
    ("DE-LU", "DK2"): (-20, -4),
    ("DE-LU", "DK1"): (-27, -1),
    ("DK1", "DK2"): (-6, 20),
    ("DK1", "NO2"): (-20, 8),
    ("GB00", "NO2"): (-12, 6),
    ("FI", "SE1"): (11, 16),
    ("FI", "SE3"): (8, 9),
    ("LT", "SE4"): (15, 5),
    ("NO2", "NO5"): (-19, 12),
    ("DK1", "GB00"): (-10, -12),
    ("EE", "FI"): (15, 4),
}
BALTIC = {
    "DK1",
    "DK2",
    "DE-LU",
    "BEI_HUB",
    "FI",
    "EE",
    "LV",
    "LT",
    "NO1",
    "NO2",
    "NO3",
    "NO4",
    "NO5",
    "SE1",
    "SE2",
    "SE3",
    "SE4",
    "PL00",
    "GB00",
    "NL00",
}
INK = "#132c40"
MUTED = "#4a6172"
SEA = "#e9f2f7"
STATE = "#f4f1e8"
OTHER = "#f8f8f5"
BORDER = "#bcc9cc"
REVIEWED = "#176277"
NATIVE = "#708794"
PROJECT = "#b65a26"
UP = "#13756b"
DOWN = "#ab5029"
PROJECT_ZONE = "BEI_HUB"
PROJECTOR = Transformer.from_crs(4326, 3035, always_xy=True)


def zone_label(zone):
    return LABELS.get(zone, zone.removesuffix("00"))


def projected(lon, lat):
    x, y = PROJECTOR.transform(lon, lat)
    return np.array([x, y]) / 1000


def read_geojson(path, url):
    if path:
        return json.loads(Path(path).read_text())
    with urlopen(url, timeout=45) as response:
        return json.load(response)


def basemap(zones_geojson, land_geojson):
    areas = defaultdict(list)
    for feature in zones_geojson["features"]:
        key = feature["properties"]["zone_name"]
        zone = "DK2" if key == "DKBH" else transmission.ALIASES.get(key, key)
        areas[zone].append(make_valid(shape(feature["geometry"])))
    clipped = {}
    for key, geometries in areas.items():
        geom = union_all(geometries)
        if not geom.is_empty and geom.bounds[0] < 45 and geom.bounds[2] > -22:
            clipped[key] = geom
    neighbors = []
    for feature in land_geojson["features"]:
        geom = make_valid(shape(feature["geometry"]))
        if geom.bounds[0] < 46 and geom.bounds[2] > -22 and geom.bounds[3] > 30:
            neighbors.append(geom)
    return clipped, neighbors


def grid(n, year):
    meta = n.meta.get(transmission.AUDIT_KEY, {})
    if meta.get("year") != year:
        raise ValueError(f"Expected {year} transmission assumptions on the network")
    corridors = {}
    for item in meta["corridors"]:
        a, b = item["zone0"], item["zone1"]
        values = np.asarray(item.get("hourly_limits_MW", item["limits_MW"]))
        corridors[tuple(sorted([a, b]))] = (
            values if a < b else values[..., ::-1],
            item,
        )
    bus_zones = transmission.zones(n)
    electricity = n.buses.loc[n.buses.carrier.isin(["AC", "DC", ""])]
    anchors = electricity.groupby(bus_zones).agg({"x": "mean", "y": "mean"})
    for z, xy in ANCHORS.items():
        if z in anchors.index:
            anchors.loc[z] = xy
    accum = defaultdict(lambda: np.zeros((len(n.snapshots), 2)))
    carriers = defaultdict(set)
    for component, frame, nominal in (
        ("Line", n.lines, "s_nom"),
        ("Link", n.links, "p_nom"),
    ):
        selected = frame.bus0.isin(electricity.index) & frame.bus1.isin(
            electricity.index
        )
        if component == "Link":
            selected &= frame.carrier.isin(["AC", "DC", "B2B"])
        for name, branch in frame.loc[selected].iterrows():
            a, b = bus_zones[branch.bus0], bus_zones[branch.bus1]
            if a == b:
                continue
            capacity = branch[nominal]
            if branch[nominal + "_extendable"]:
                capacity = (
                    branch[nominal + "_opt"]
                    if n.is_solved
                    else branch[nominal + "_min"]
                )
            if capacity <= 0:
                continue
            if component == "Line":
                pulse = transmission.dense(n, component, "s_max_pu", [name]).to_numpy()
                flow = capacity * np.repeat(pulse, 2, axis=1)
            else:
                forward = transmission.dense(
                    n, component, "p_max_pu", [name]
                ).to_numpy()
                reverse = -transmission.dense(
                    n, component, "p_min_pu", [name]
                ).to_numpy()
                efficiency = transmission.dense(
                    n, component, "efficiency", [name]
                ).to_numpy()
                flow = capacity * np.maximum(
                    0, np.concatenate([forward, reverse * efficiency], axis=1)
                )
            pair = tuple(sorted([a, b]))
            accum[pair] += flow if a < b else flow[:, ::-1]
            carriers[pair].add("AC" if component == "Line" else branch.carrier)

    rows = []
    for (a, b), native in sorted(accum.items()):
        is_project = a.startswith("BEI_") or b.startswith("BEI_")
        has_limit = (a, b) in corridors
        if has_limit:
            limits, record = corridors[a, b]
            plotted = np.minimum(native, limits)
            source = record["source"]
            published = np.quantile(np.broadcast_to(limits, native.shape), 0.95, axis=0)
        else:
            plotted = native
            source = (
                "Bornholm project configuration"
                if is_project
                else "PyPSA-Eur physical branch inputs"
            )
            published = (np.nan, np.nan)
        p95 = np.quantile(plotted, 0.95, axis=0)
        native_p95 = np.quantile(native, 0.95, axis=0)
        rows.append(
            {
                "year": year,
                "zone0": a,
                "zone1": b,
                "out_GW": p95[0] / 1000,
                "back_GW": p95[1] / 1000,
                "native_out_GW": native_p95[0] / 1000,
                "native_back_GW": native_p95[1] / 1000,
                "published_out_GW": published[0] / 1000,
                "published_back_GW": published[1] / 1000,
                "type": "+".join(sorted(carriers[a, b])),
                "basis": "Project"
                if is_project
                else "Reviewed"
                if has_limit
                else "Native",
                "source": source,
            }
        )
    return anchors, pd.DataFrame(rows).set_index(["zone0", "zone1"])


def geometry_path(geom):
    if geom.is_empty:
        return
    if geom.geom_type == "Polygon":
        from shapely.geometry.polygon import orient

        geom = orient(geom, sign=1)
        paths = []
        for ring in (geom.exterior, *geom.interiors):
            coords = np.asarray(ring.coords)
            codes = [
                MplPath.MOVETO,
                *([MplPath.LINETO] * (len(coords) - 2)),
                MplPath.CLOSEPOLY,
            ]
            paths.append(MplPath(coords, codes))
        yield MplPath.make_compound_path(*paths)
    elif hasattr(geom, "geoms"):
        for polygon in geom.geoms:
            yield from geometry_path(polygon)


def make_paths(areas, neighbors):
    def project(coordinates):
        return projected(*coordinates.T).T

    window = box(-28, 29, 47, 74)

    def convert(geom):
        geom = geom.intersection(window)
        return list(
            geometry_path(
                transform(geom, project).simplify(2.4, preserve_topology=True)
            )
        )

    return (
        {k: convert(v) for k, v in areas.items()},
        [p for g in neighbors for p in convert(g)],
    )


def format_map(ax, paths, extent, active):
    ax.set_facecolor(SEA)
    ax.set(xlim=extent[:2], ylim=extent[2:], aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#cedbe1")
        spine.set_linewidth(0.75)
    zones, neighbors = paths
    for p in neighbors:
        ax.add_patch(
            PathPatch(p, facecolor=OTHER, edgecolor="#d7dddb", linewidth=0.38, zorder=0)
        )
    for zone, shapes in zones.items():
        if zone in active:
            for p in shapes:
                ax.add_patch(
                    PathPatch(
                        p, facecolor=STATE, edgecolor=BORDER, linewidth=0.55, zorder=1
                    )
                )


def line_geometry(p0, p1, pair):
    p0, p1 = np.asarray(p0), np.asarray(p1)
    delta = p1 - p0
    center = (p0 + p1) / 2 + BEND.get(pair, 0) * np.array([-delta[1], delta[0]])
    a, b = p0, p1
    if pair == ("BEI_HUB", "DE-LU"):
        center = projected(13.62, 54.14) * 0.65 + center * 0.35
    return MplPath([a, center, b], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]), (
        a + 2 * center + b
    ) / 4


def width(gw):
    return 0.34 + 1.38 * np.sqrt(gw)


def plot_map(ax, paths, extent, anchors, rows, mode):
    format_map(ax, paths, extent, set(anchors.index))
    positions = {zone: projected(x, y) for zone, (x, y) in anchors.iterrows()}
    direct_labels = mode == "numbers"
    ordering = rows.assign(best=rows[["out_GW", "back_GW"]].max(axis=1)).sort_values(
        ["basis", "best"], ascending=[True, True]
    )
    for (a, b), row in ordering.iterrows():
        pair = (a, b)
        if not row.best > 0:
            continue
        path, middle = line_geometry(positions[a], positions[b], pair)
        if mode == "change":
            value = row.delta_GW
            color = UP if value > 0 else DOWN
            visible = abs(value) >= 0.25
            lw = width(abs(value)) if visible else 0.42
            color = color if visible else "#a8b7bf"
            zorder = 6 if visible else 3
        else:
            color = {"Reviewed": REVIEWED, "Native": NATIVE, "Project": PROJECT}[
                row.basis
            ]
            lw = (
                (1.35 if row.basis == "Project" else 0.9)
                if direct_labels
                else width(row.best)
            )
            visible = not direct_labels or pair in NUMBER_OFFSETS
            zorder = (
                6 if row.basis == "Project" else 5 if row.basis == "Reviewed" else 4
            )
        alpha = 0.94 if visible else 0.52
        ax.add_patch(
            PathPatch(
                path,
                facecolor="none",
                edgecolor=color,
                linewidth=lw,
                zorder=zorder,
                capstyle="round",
                alpha=alpha,
            )
        )
        if direct_labels and pair in NUMBER_OFFSETS and a in BALTIC and b in BALTIC:
            dx, dy = NUMBER_OFFSETS[pair]
            ax.annotate(
                f"{row.best:.2f}".rstrip("0").rstrip("."),
                xy=middle,
                xytext=(dx, dy),
                textcoords="offset points",
                ha="center",
                va="center",
                color=color,
                fontsize=9.3,
                weight="bold",
                zorder=11,
                annotation_clip=True,
                bbox={
                    "boxstyle": "round,pad=.19,rounding_size=.30",
                    "fc": "white",
                    "ec": "#cdd9db",
                    "lw": 0.55,
                    "alpha": 0.98,
                },
                arrowprops={
                    "arrowstyle": "-",
                    "color": color,
                    "lw": 0.58,
                    "alpha": 0.75,
                    "shrinkA": 2,
                    "shrinkB": 1,
                },
            )
    for z, (x, y) in positions.items():
        if not (extent[0] <= x <= extent[1] and extent[2] <= y <= extent[3]):
            continue
        project = z.startswith("BEI_")
        ax.scatter(
            x,
            y,
            s=39 if project else 11,
            c=PROJECT if project else INK,
            marker="D" if project else "o",
            edgecolors="white",
            linewidths=0.55,
            zorder=8,
        )
        if mode == "numbers" and z not in BALTIC | {"GBZE"}:
            continue
        dx, dy = LABEL_OFFSETS.get(z, (5, 5))
        ax.annotate(
            zone_label(z),
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left" if dx > 0 else "right" if dx < 0 else "center",
            va="center",
            weight="bold",
            fontsize=8.4 if direct_labels else 6.6,
            color=PROJECT if project else INK,
            zorder=9,
            annotation_clip=True,
            path_effects=[pe.withStroke(linewidth=2.1, foreground=STATE)],
        )


EUROPE = (2480, 6350, 1280, 5600)
NORDICS = (3270, 5650, 3020, 5500)


def header(fig, title, subtitle):
    fig.text(0.047, 0.948, title, color=INK, fontsize=24, weight="bold", va="top")
    fig.text(0.048, 0.890, subtitle, color=MUTED, fontsize=11, va="top")
    fig.add_artist(
        Line2D(
            [0.048, 0.952],
            [0.85, 0.85],
            color="#d4e0e2",
            lw=0.85,
            transform=fig.transFigure,
        )
    )


def footer(fig, lines):
    for k, line in enumerate(lines):
        fig.text(
            0.048, 0.061 - k * 0.024, line, color=MUTED, fontsize=8.25, va="center"
        )


def pair_page(paths, cases, mode, status):
    number = mode == "numbers"
    fig = plt.figure(figsize=(16, 9.5), facecolor="white")
    if number:
        header(
            fig,
            "Nordic and Baltic transfer limits",
            "Direct labels show the larger of the two directional limits, in GW.",
        )
    else:
        header(
            fig,
            "The European electricity grid in the model",
            "Bidding zones and modelled transfer limits for the 2025 hindcast and the 2035 Bornholm treaty hybrid.",
        )
    for i, (year, (anchors, rows)) in enumerate(cases.items()):
        left = 0.049 + 0.464 * i
        fig.text(
            left, 0.821, str(year), fontsize=17, weight="bold", color=INK, va="center"
        )
        caption = (
            "Historical fleet • hourly limits"
            if year == 2025
            else "Selected TYNDP 2026 projects • treaty hybrid"
        )
        fig.text(left + 0.068, 0.821, caption, fontsize=10, color=MUTED, va="center")
        ax = fig.add_axes([left, 0.143, 0.438, 0.645])
        plot_map(ax, paths, NORDICS if number else EUROPE, anchors, rows, mode)

    handles = [
        Line2D([], [], color=color, lw=2, label=label)
        for color, label in (
            (REVIEWED, "Reviewed corridor"),
            (NATIVE, "Native branch limit"),
            (PROJECT, "Bornholm cable"),
        )
    ]
    if not number:
        handles.extend(
            Line2D([], [], color=INK, lw=width(gw), label=f"{gw:g} GW")
            for gw in (1, 5, 10, 20)
        )
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.095),
        ncol=len(handles),
        frameon=False,
        labelcolor=INK,
        fontsize=9.2,
        handlelength=1.7,
        columnspacing=2.0,
    )
    footer(
        fig,
        [
            "2025 shows the 95th percentile of hourly model limits. 2035 uses fixed ERAA and grid-operator inputs where available. "
            "Branch ratings also bound each corridor.",
            status
            + "  •  Schematic routes  •  Zone and coast geometry: Electricity Maps, ENTSO-E and Natural Earth",
        ],
    )
    return fig


def change_page(paths, cases, status):
    fig = plt.figure(figsize=(16, 9.5), facecolor="white")
    header(
        fig,
        "Change in reviewed transfer limits",
        "2035 minus 2025 in the larger direction. The table gives both directions in GW.",
    )
    a25, y25 = cases[2025]
    a35, y35 = cases[2035]
    paired = y25[["out_GW", "back_GW"]].join(
        y35[["out_GW", "back_GW"]], how="outer", lsuffix="_25", rsuffix="_35"
    )
    delta = paired[["out_GW_35", "back_GW_35"]].max(axis=1).fillna(0) - paired[
        ["out_GW_25", "back_GW_25"]
    ].max(axis=1).fillna(0)
    old_basis = y25["basis"].reindex(paired.index)
    new_basis = y35["basis"].reindex(paired.index)
    compared = new_basis.eq("Reviewed") & old_basis.eq("Reviewed")
    new_links = old_basis.isna() & new_basis.isin(["Reviewed", "Project"])
    changes = y35.combine_first(y25).assign(
        delta_GW=delta.where(compared | new_links, 0)
    )
    ax = fig.add_axes([0.048, 0.135, 0.53, 0.67])
    plot_map(ax, paths, EUROPE, a35.combine_first(a25), changes, "change")

    table_ax = fig.add_axes([0.595, 0.174, 0.355, 0.632])
    table_ax.axis("off")
    keys = [
        ("SE2", "SE3"),
        ("SE3", "SE4"),
        ("DK2", "SE4"),
        ("DE-LU", "DK1"),
        ("DE-LU", "DK2"),
        ("DK1", "DK2"),
        ("DK1", "NO2"),
        ("FI", "SE1"),
        ("NO2", "NO5"),
        ("BEI_HUB", "DK2"),
        ("BEI_HUB", "DE-LU"),
    ]
    for x, label, ha in [
        (0.025, "Corridor", "left"),
        (0.54, "2025  → / ←", "center"),
        (0.84, "2035  → / ←", "center"),
    ]:
        table_ax.text(
            x, 0.955, label, ha=ha, va="center", color=INK, fontsize=11, weight="bold"
        )
    for i, key in enumerate(keys):
        if key not in paired.index:
            continue
        pos = 0.879 - i * 0.074
        table_ax.axhspan(
            pos - 0.037, pos + 0.037, color="#edf4f5" if i % 2 == 0 else "#f7fafb", lw=0
        )
        table_ax.text(
            0.025,
            pos,
            f"{zone_label(key[0])}–{zone_label(key[1])}",
            va="center",
            color=PROJECT if key[0] == PROJECT_ZONE else INK,
            fontsize=10.2,
        )
        for year, col in (("25", 0.54), ("35", 0.84)):
            forward, backward = paired.loc[key, [f"out_GW_{year}", f"back_GW_{year}"]]

            shown = "—" if pd.isna(forward) else f"{forward:.3f} / {backward:.3f}"
            table_ax.text(
                col, pos, shown, va="center", ha="center", color=INK, fontsize=10.4
            )
    fig.legend(
        handles=[
            Line2D([], [], color=UP, lw=2.6, label="Higher in 2035"),
            Line2D([], [], color=DOWN, lw=2.6, label="Lower in 2035"),
            Line2D([], [], color="#a8b7bf", lw=1.2, label="Other corridors"),
        ],
        loc="lower left",
        bbox_to_anchor=(0.075, 0.082),
        ncol=3,
        frameon=False,
        labelcolor=INK,
        fontsize=9.5,
    )
    footer(
        fig,
        [
            "Differences reflect the grid and the operating assumptions. "
            "Pale lines show changes under 0.25 GW or corridors outside the reviewed comparison.",
            status
            + "  •  Schematic routes  •  Zone and coast geometry: Electricity Maps, ENTSO-E and Natural Earth",
        ],
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network-2025", type=Path, required=True)
    parser.add_argument("--network-2035", type=Path, required=True)
    parser.add_argument(
        "--zones", type=Path, required=True, help="PyPSA-Eur bidding_zones.geojson"
    )
    parser.add_argument(
        "--land", type=Path, help="Local copy of the Natural Earth countries"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--source-note", default="Sources: 2025 and 2035 PyPSA network files"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cases = {
        y: grid(pypsa.Network(path), y)
        for y, path in ((2025, args.network_2025), (2035, args.network_2035))
    }
    paths = make_paths(
        *basemap(json.loads(args.zones.read_text()), read_geojson(args.land, LAND_URL))
    )
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    mpl.rcParams.update(
        {
            "font.family": ["Open Sans", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "figure.autolayout": False,
        }
    )
    with PdfPages(args.output / "grid_capacity_comparison.pdf") as pdf:
        for mode, name in (
            ("width", "europe"),
            ("numbers", "nordics"),
            ("change", "changes"),
        ):
            fig = (
                change_page(paths, cases, args.source_note)
                if mode == "change"
                else pair_page(paths, cases, mode, args.source_note)
            )
            pdf.savefig(fig)
            if mode != "change":
                fig.savefig(args.output / f"grid_capacity_{name}.png", dpi=200)
            plt.close(fig)
    pd.concat([rows for _, rows in cases.values()]).assign(
        network_basis=args.source_note
    ).to_csv(args.output / "grid_capacity_corridors.csv", float_format="%.4f")


if __name__ == "__main__":
    main()
