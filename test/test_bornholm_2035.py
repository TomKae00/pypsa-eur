"""Physical regressions for the final 2035 matrix and separate EV participation."""

import os
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytest
import xarray as xr
import yaml

from scripts import bornholm_2035 as campaign
from scripts.add_energy_islands import main as add_project
from scripts.prepare_sector_network import add_EVs, lossy_bidirectional_links


def ev_network(v2g_share):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=6, freq="h"))
    nodes = pd.Index(["DE00"])
    n.add("Bus", nodes, carrier="AC", country="DE", location=nodes)
    profile = pd.DataFrame(10.0, index=n.snapshots, columns=nodes)
    options = yaml.safe_load(open("config/config.default.yaml"))["sector"]
    options["bev_v2g_availability"] = v2g_share
    add_EVs(
        n,
        profile * 0 + 0.8,
        profile * 0 + 0.2,
        profile,
        0.45,
        pd.Series(1000, index=nodes),
        profile * 0 + 18,
        SimpleNamespace(nodes=nodes),
        options,
    )
    return n


@pytest.mark.parametrize("share", [0.0, 0.25, 0.5])
def test_v2g_changes_discharge_access_but_preserves_mobility_charging_and_energy(share):
    n = ev_network(share)
    baseline = ev_network(0.0)
    assert n.stores.e_nom.sum() == pytest.approx(1000 * 0.45 * 0.05 * 0.5)
    assert n.links.loc[
        n.links.carrier.eq("BEV charger"), "p_nom"
    ].sum() == pytest.approx(1000 * 0.45 * 0.011)
    assert n.links.loc[n.links.carrier.eq("V2G"), "p_nom"].sum() == pytest.approx(
        1000 * 0.45 * 0.011 * share
    )
    np.testing.assert_allclose(
        n.loads_t.p_set.sum(axis=1), baseline.loads_t.p_set.sum(axis=1)
    )
    if share:
        battery = n.links.loc[n.links.carrier.eq("V2G"), "bus0"].iloc[0]
        energy = n.stores.loc[n.stores.bus.eq(battery), "e_nom"].sum()
        assert energy == pytest.approx(1000 * 0.45 * 0.05 * share)
        # DSM-only energy has no direct discharge connection to the grid.
        other_batteries = set(n.stores.bus) - {battery}
        assert (
            not n.links.loc[n.links.carrier.eq("V2G"), "bus0"]
            .isin(other_batteries)
            .any()
        )


def test_unspecified_v2g_share_retains_original_names_and_pooled_fleet():
    n = ev_network(None)
    assert set(n.stores.index) == {"DE00 EV battery"}
    assert "DE00 V2G" in n.links.index
    assert n.links.at["DE00 V2G", "p_nom"] == pytest.approx(2.475)


def test_v2g_cannot_exceed_managed_fleet():
    with pytest.raises(ValueError, match="DSM share"):
        ev_network(0.6)


@pytest.mark.parametrize("share", [0.0, 0.25, 0.5])
def test_ev_cohorts_supply_mobility_in_a_physical_dispatch(share):
    n = ev_network(share)
    n.add(
        "Generator",
        "grid",
        bus="DE00",
        p_nom=100,
        marginal_cost=[10, 10, 100, 100, 10, 10],
    )
    result = n.optimize(solver_name="highs", threads=1, log_to_console=False)
    assert result == ("ok", "optimal")
    assert np.isfinite(n.objective)
    n.meta["sector"] = {"bev_charge_efficiency": 0.9}
    use = campaign.electricity_use(n)
    assert len(use) == 1
    assert use[0]["electricity_TWh"] == pytest.approx(
        n.loads_t.p_set.sum().sum() / 0.9 / 1e6
    )
    if share == 0:
        assert not n.links.carrier.eq("V2G").any()


def test_matrix_has_only_two_additional_single_country_radial_tests():
    config = campaign.read_config()
    cases = campaign.scenarios(config)
    assert len(cases) == 15
    assert [c["index"] for c in cases] == list(range(15))
    assert yaml.safe_load(campaign.SCENARIOS.read_text()) == campaign.definitions(
        config
    )
    radial = [c for c in cases if c["project"] in ("radial_dk", "radial_de")]
    assert len(radial) == 2
    assert all(c["world"] == "central" and c["v2g_share"] == 0.25 for c in radial)


def test_additional_data_centre_energy_preserves_other_service_demands():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=4, freq="h"))
    n.snapshot_weightings.loc[:, :] = 2190
    n.add("Bus", ["DKW1", "DKE1", "DE00"], carrier="AC", country=["DK", "DK", "DE"])
    n.add("Bus", "DE heat", carrier="heat", location="DE00")
    n.add("Load", "heat", bus="DE heat", carrier="heat", p_set=[1, 2, 3, 4])
    before = n.loads_t.p_set.copy()
    campaign.apply(
        n,
        {
            "enabled": True,
            "additional_data_centres_TWh": {"DK1": 2, "DK2": 1, "DE": 40},
        },
    )
    assert n.loads_t.p_set["heat"].equals(before["heat"])
    report = pd.DataFrame(campaign.demand_report(n))
    extra = report[report.carrier.eq("additional data centres")].set_index("country")
    assert extra.at["DK", "final_demand_TWh"] == pytest.approx(3)
    assert extra.at["DE", "final_demand_TWh"] == pytest.approx(40)
    with pytest.raises(ValueError, match="already applied"):
        campaign.apply(n, {"enabled": True, "additional_data_centres_TWh": {}})


@pytest.mark.parametrize("losses", [False, True])
def test_split_radial_cannot_transfer_the_danish_wind_block_to_germany(
    tmp_path, losses
):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=4, freq="h"))
    n.add("Carrier", ["AC", "DC", "offwind-ac"])
    n.add("Bus", ["DKE1", "DE00"], carrier="AC", country=["DK", "DE"])
    n.add("Generator", "backup", bus="DE00", p_nom=5000, marginal_cost=100)
    n.add("Load", "German demand", bus="DE00", p_set=2000)
    ds = xr.Dataset(
        {
            "profile": (("time", "bus", "bin"), np.full((4, 1, 1), 0.5)),
            "average_distance": (("bus", "bin"), [[10.0]]),
        },
        coords={"time": n.snapshots.to_numpy(), "bus": ["BEI_WIND"], "bin": [0]},
    )
    profile = tmp_path / "profile.nc"
    ds.to_netcdf(profile)
    paths = {
        f"energy_island_{name}": f"data/energy_islands/bornholm_{name}.csv"
        for name in ("hubs", "wind", "links")
    }
    paths["energy_island_profile"] = str(profile)
    costs = pd.DataFrame(
        {
            "capital_cost": 1.0,
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
    add_project(
        n,
        paths,
        {
            "enabled": True,
            "case": "split_radial",
            "profile_mode": "project",
            "exclude_project_area": True,
            "subtract_generic_potential": False,
        },
        costs,
        2035,
    )
    if losses:
        lossy_bidirectional_links(
            n, "DC", {"efficiency_static": 0.98, "efficiency_per_1000km": 0.977}
        )
    campaign.project_check(n, "split_radial")
    result = n.optimize(solver_name="highs", threads=1, log_to_console=False)
    assert result == ("ok", "optimal")
    np.testing.assert_allclose(n.generators_t.p["BEI_WIND_DK"], 0, atol=1e-6)
    np.testing.assert_allclose(n.generators_t.p["BEI_WIND_DE"], 1000, atol=1e-6)
    delivered = 1000 * n.links.at["BEI_TO_DE_SPLIT", "efficiency"]
    np.testing.assert_allclose(n.generators_t.p["backup"], 2000 - delivered, atol=1e-6)
    restored = tmp_path / "network.nc"
    n.export_to_netcdf(restored)
    campaign.project_check(pypsa.Network(restored), "split_radial")


@pytest.mark.parametrize("prepare_exit", [0, 19])
def test_controller_never_submits_solves_before_all_input_checks_pass(
    tmp_path, prepare_exit
):
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    for name in ("bornholm_2035_workq.sbatch", "bornholm_workq_bridge.sh"):
        shutil.copy2(Path("scripts") / name, scripts / name)
    tools = tmp_path / "bin"
    tools.mkdir()
    wrapper = tools / "pixi"
    wrapper.write_text("""#!/usr/bin/env bash
set -eu
case "$*" in
  *--list-build-targets*) for i in {0..14}; do echo "resources/case_${i}.nc"; done ;;
  *--prepare*) exit "${FAKE_PREPARE_EXIT}" ;;
esac
""")
    sbatch = tools / "sbatch"
    sbatch.write_text("""#!/usr/bin/env bash
printf '%s\\n' "$*" >> "$FAKE_SUBMISSIONS"
echo 800
""")
    cancel = tools / "scancel"
    cancel.write_text("#!/usr/bin/env bash\nexit 0\n")
    for path in (wrapper, sbatch, cancel):
        path.chmod(0o755)
    journal = tmp_path / "submissions.txt"
    env = dict(
        os.environ,
        PATH=f"{tools}:{os.environ['PATH']}",
        PIXI_SOPHIA=str(wrapper),
        SLURM_SUBMIT_DIR=str(tmp_path),
        SLURM_JOB_ID="123",
        SLURM_JOB_PARTITION="workq",
        FAKE_PREPARE_EXIT=str(prepare_exit),
        FAKE_SUBMISSIONS=str(journal),
    )
    result = subprocess.run(
        ["bash", str(scripts / "bornholm_2035_workq.sbatch"), "all"],
        env=env,
        text=True,
        capture_output=True,
        timeout=20,
    )
    if prepare_exit:
        assert result.returncode == prepare_exit
        assert not journal.exists()
    else:
        assert result.returncode == 0, result.stderr
        calls = journal.read_text().splitlines()
        assert len(calls) == 2
        assert "--array=0-14%6" in calls[0]
        assert "--dependency=afterany:800" in calls[1]
