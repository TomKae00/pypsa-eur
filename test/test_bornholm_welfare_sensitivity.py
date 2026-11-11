# SPDX-License-Identifier: MIT
"""Meaningful small-network and submission tests; full annual solves run on Sophia."""

import json
import subprocess
from pathlib import Path

import pandas as pd
import pypsa
import pytest

from scripts import bornholm_welfare_sensitivity as sensitivity
from scripts import run_bornholm_welfare as runner


def test_timestamp_hash_preserves_precision_without_version_specific_keyword():
    from scripts.transmission_assumptions import _snapshot_hash

    index = pd.DatetimeIndex(["2035-01-01 00:00:00.000000001"]).as_unit("ns")
    _snapshot_hash(index)
    with pytest.raises(ValueError, match="lose precision"):
        _snapshot_hash(index, "us")


@pytest.fixture
def config():
    return runner.read_config()


@pytest.fixture
def network():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=2, freq="3h"))
    n.snapshot_weightings.loc[:, :] = 3.0
    for name, country in (
        ("DKW1", "DK"),
        ("DKE1", "DK"),
        ("DE00", "DE"),
        ("SE04", "SE"),
    ):
        n.add("Bus", name, carrier="AC", country=country, location=name)
        n.add(
            "Bus",
            name + " low voltage",
            carrier="low voltage",
            country=country,
            location=name,
        )
        n.add("Load", name, bus=name + " low voltage", p_set=[100, 200])
        n.add("Load", name + " industry", bus=name + " low voltage", p_set=10)
        n.add(
            "Bus", name + " heat", carrier="rural heat", country=country, location=name
        )
        n.add("Load", name + " heat", bus=name + " heat", p_set=[30, 40])
        for j in (0, 1):
            n.add(
                "Generator",
                name + f" wind {j}",
                bus=name,
                carrier="onwind",
                p_nom_extendable=True,
                p_nom_min=5,
                capital_cost=1,
            )
        n.add(
            "Generator", name + " wind existing", bus=name, carrier="onwind", p_nom=10
        )
        n.add(
            "Generator",
            name + " PV",
            bus=name,
            carrier="solar",
            p_nom_extendable=True,
            capital_cost=2,
        )
        n.add(
            "Generator",
            name + " tracking",
            bus=name,
            carrier="solar-hsat",
            p_nom_extendable=True,
            capital_cost=2,
        )
        n.add(
            "Generator",
            name + " roof",
            bus=name + " low voltage",
            carrier="solar rooftop",
            p_nom_extendable=True,
            capital_cost=2,
        )
    return n


def test_campaign_has_exact_requested_demand_and_adjacent_pairs(config):
    cases = sensitivity.scenarios(config)
    assert len(cases) == 12
    assert {s["demand_factor"] for s in cases} == {1.0, 1.1, 1.2}
    for a, b in zip(cases[::2], cases[1::2]):
        assert a["pair"] == b["pair"]
        assert (a["project"], b["project"]) == ("baseline", "hybrid")
        assert a["index"] + 1 == b["index"]
        assert sensitivity.capacity_limits(config, a) == sensitivity.capacity_limits(
            config, b
        )


def test_af26_and_published_ceiling_derivations(config):
    assert config["deployment"]["central"]["onshore_DK1_MW"] + config["deployment"][
        "central"
    ]["onshore_DK2_MW"] == pytest.approx(7072.7)
    assert config["deployment"]["slow_onshore"]["onshore_DK1_MW"] + config[
        "deployment"
    ]["slow_onshore"]["onshore_DK2_MW"] == pytest.approx(6572.7)
    assert config["deployment"]["slow_onshore"]["onshore_DE_MW"] == pytest.approx(
        (115 + 5 / 7 * (126.6 - 115)) * 1000
    )
    assert config["solar"]["DE_total_MW"] == pytest.approx(
        (215 + 5 / 7 * (345.4 - 215)) * 1000
    )


@pytest.mark.parametrize("factor", [1.0, 1.1, 1.2])
def test_fixed_low_voltage_and_static_industry_scale_but_not_other_end_uses(
    config, network, factor
):
    case = next(
        s for s in sensitivity.scenarios(config) if s["demand_factor"] == factor
    )
    before = network.loads_t.p_set.copy()
    audit = sensitivity.apply(network, config, case)
    for region in ("DKW1", "DKE1", "DE00"):
        assert network.loads_t.p_set[region].equals(before[region] * factor)
        assert network.loads.at[region + " industry", "p_set"] == 10 * factor
        assert network.loads_t.p_set[region + " heat"].equals(before[region + " heat"])
    assert network.loads_t.p_set.SE04.equals(before.SE04)
    assert audit["demand"]["DK"]["after_MWh"] == pytest.approx(
        audit["demand"]["DK"]["before_MWh"] * factor
    )
    sensitivity.validate(network, config, case)
    with pytest.raises(ValueError, match="already applied"):
        sensitivity.apply(network, config, case)


def test_prepared_demand_tampering_is_rejected(config, network):
    case = sensitivity.scenarios(config)[0]
    sensitivity.apply(network, config, case)
    network.loads_t.p_set.loc[network.snapshots[0], "DKE1"] += 1
    with pytest.raises(ValueError, match="Demand changed"):
        sensitivity.validate(network, config, case)


@pytest.mark.parametrize("unit", ["ns", "us", "s"])
def test_prepared_demand_audit_survives_netcdf_roundtrip(
    config, network, tmp_path, unit
):
    network.set_snapshots(network.snapshots.as_unit(unit))
    case = sensitivity.scenarios(config)[10]
    sensitivity.apply(network, config, case)
    path = tmp_path / "composed.nc"
    network.export_to_netcdf(path)
    restored = pypsa.Network(path)
    sensitivity.validate(restored, config, case)


def test_solar_selection_includes_all_bins_and_respects_roof_boundary(network):
    ground = sensitivity.selected_generators(
        network, dict(region="DK2", carriers=["solar", "solar-hsat"])
    )
    assert set(ground) == {"DKE1 PV", "DKE1 tracking"}
    roof = sensitivity.selected_generators(
        network, dict(region="DK2", carriers=["solar rooftop"])
    )
    assert list(roof) == ["DKE1 roof"]
    assert "DKW1 PV" not in ground


def test_capacity_below_committed_background_fails_before_solve(network):
    limit = dict(id="wind_DK1", region="DK1", carriers=["onwind"], maximum_MW=19.0)
    with pytest.raises(ValueError, match="existing/minimum capacity"):
        sensitivity.capacity_report(network, [limit])


def test_country_cap_is_enforced_as_one_sum_including_existing():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035-01-01", periods=2, freq="h"))
    n.add("Bus", "DE00", carrier="AC", country="DE", location="DE00")
    n.add("Load", "load", bus="DE00", p_set=100)
    n.add("Generator", "existing", bus="DE00", carrier="onwind", p_nom=10)
    for j in range(3):
        n.add(
            "Generator",
            f"bin{j}",
            bus="DE00",
            carrier="onwind",
            p_nom_extendable=True,
            capital_cost=1,
            marginal_cost=0,
        )
    n.add("Generator", "gas", bus="DE00", carrier="gas", p_nom=100, marginal_cost=100)
    settings = dict(
        capacity_limits=[
            dict(id="DE_onshore", region="DE", carriers=["onwind"], maximum_MW=30)
        ]
    )
    status, condition = n.optimize(
        solver_name="highs",
        threads=1,
        extra_functionality=lambda n, sns: sensitivity.add_capacity_constraints(
            n, settings
        ),
        include_objective_constant=True,
    )
    assert (status, condition) == ("ok", "optimal")
    report = sensitivity.capacity_report(n, settings["capacity_limits"], solved=True)[0]
    assert report["solved_MW"] == pytest.approx(30)
    assert n.generators_t.p.gas.iloc[0] == pytest.approx(70)


def test_pairs_compare_every_final_demand_even_negative_carbon_load(config, network):
    a = network
    a.add("Bus", "co2", carrier="co2", country="")
    a.add("Load", "emissions", bus="co2", p_set=-3)
    b = a.copy()
    runner.check_pair(a, b)
    b.loads.at["DKE1 industry", "p_set"] += 1
    with pytest.raises(ValueError, match="final-demand"):
        runner.check_pair(a, b)


def test_cost_overlay_changes_only_annuitised_capex(config):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035-01-01", periods=2, freq="h"))
    n.add("Bus", "BEI_HUB", carrier="AC")
    n.add(
        "Generator",
        "BEI_WIND",
        bus="BEI_HUB",
        carrier="offwind-ac",
        p_nom=3000,
        capital_cost=300000,
        marginal_cost=0.015,
        p_max_pu=[0.5, 1],
    )
    n.generators.loc["BEI_WIND", "p_nom_opt"] = 3000
    n.generators_t.p = pd.DataFrame({"BEI_WIND": [1500, 2900]}, index=n.snapshots)
    n.buses_t.marginal_price = pd.DataFrame({"BEI_HUB": [20, 80]}, index=n.snapshots)
    dispatch = n.generators_t.p.copy()
    rows = runner.cost_cases(n, config)
    assert len(rows) == 3
    assert n.generators_t.p.equals(dispatch)
    assert rows[0]["annual_wind_collection_cost_EUR"] == pytest.approx(
        900000000 + 0.015 * 4400
    )
    annuity = 0.07 / (1 - 1.07**-30)
    assert rows[2]["annual_wind_collection_cost_EUR"] - rows[0][
        "annual_wind_collection_cost_EUR"
    ] == pytest.approx((3480 - 2039.704) * 1000 * 3000 * annuity)
    assert not rows[2]["export_cost_included"]


@pytest.mark.parametrize("prepare_fails", [False, True])
def test_controller_submits_twelve_and_collects_afterany_only_after_preparation(
    tmp_path, prepare_fails
):
    log = tmp_path / "calls"
    wrapper = tmp_path / "pixi-sophia"
    wrapper.write_text(
        '#!/bin/bash\nprintf \'%s\\n\' "$*" >> "$CALLS"\n'
        + ("exit 1\n" if prepare_fails else "exit 0\n")
    )
    wrapper.chmod(0o755)
    for command, body in (
        ("sbatch", 'printf \'%s\\n\' "$*" >> "$CALLS"; echo 12345'),
        ("scancel", "exit 0"),
    ):
        p = tmp_path / command
        p.write_text("#!/bin/bash\n" + body + "\n")
        p.chmod(0o755)
    import os

    env = dict(
        os.environ,
        SLURM_JOB_PARTITION="workq",
        SLURM_JOB_ID="9000",
        SLURM_SUBMIT_DIR=str(tmp_path),
        PIXI_SOPHIA=str(wrapper),
        CALLS=str(log),
        PATH=str(tmp_path) + ":" + os.environ["PATH"],
    )
    script = Path("scripts/bornholm_welfare_workq.sbatch").resolve()
    result = subprocess.run(
        ["bash", str(script)], env=env, capture_output=True, text=True
    )
    calls = log.read_text()
    if prepare_fails:
        assert result.returncode != 0
        assert "--array" not in calls
    else:
        assert result.returncode == 0, result.stderr
        assert "--array=0-11%6" in calls
        assert "--dependency=afterany:12345" in calls


def test_solve_array_accepts_task_eleven_before_invoking_srun(tmp_path):
    wrapper = tmp_path / "pixi-sophia"
    wrapper.write_text("#!/bin/bash\nexit 0\n")
    wrapper.chmod(0o755)
    srun = tmp_path / "srun"
    srun.write_text("#!/bin/bash\nprintf '%s\\n' \"$*\"\n")
    srun.chmod(0o755)
    import os

    env = dict(
        os.environ,
        SLURM_JOB_PARTITION="workq",
        SLURM_JOB_ID="9000",
        SLURM_ARRAY_TASK_ID="11",
        SLURM_CPUS_PER_TASK="16",
        SLURM_SUBMIT_DIR=str(tmp_path),
        PIXI_SOPHIA=str(wrapper),
        PATH=str(tmp_path) + ":" + os.environ["PATH"],
    )
    script = Path("scripts/bornholm_welfare_solve.sbatch").resolve()
    result = subprocess.run(
        ["bash", str(script)], env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert "--index 11" in result.stdout


def test_collector_updates_candidate_link_on_retry_without_promoting_flags(
    config, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    case = sensitivity.scenarios(config)[0]
    _, directory = runner.paths(config, case)
    directory.mkdir(parents=True)
    n = pypsa.Network()
    monkeypatch.setattr(pypsa, "Network", lambda path: n)
    for attempt_id in ("first", "retry"):
        attempt = directory / "attempts" / attempt_id
        attempt.mkdir(parents=True)
        candidate = attempt / "candidate_2035.nc"
        candidate.write_text(attempt_id)
        runner.write_json(
            directory / "summary.json",
            dict(
                status="optimal_needs_review",
                scenario=case,
                attempt=str(attempt),
                candidate_network=str(candidate),
                solve_result=dict(status="ok", termination_condition="optimal"),
                screening_issues=["DualVio exceeds threshold"],
            ),
        )
        runner.collect(config)
        link = (
            Path("results")
            / config["prefix"]
            / "transfer"
            / case["name"]
            / "candidate_2035.nc"
        )
        assert link.resolve() == candidate.resolve()
        copied = json.loads((link.parent / "summary.json").read_text())
        assert copied["status"] == "optimal_needs_review"
        assert "accepted_network" not in copied
