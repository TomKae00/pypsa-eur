"""Small independent checks of audit signs, durations and accounting conventions."""

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import pytest
import yaml

from scripts import audit_bornholm_candidate as audit
from scripts.audit_bornholm_candidate import Audit, safe_product, storage_residual


def test_zero_cost_infinite_capacity_is_zero():
    actual = safe_product(np.array([0, 2, -1]), np.array([np.inf, 5, np.inf]))
    np.testing.assert_array_equal(actual, np.array([0, 10, -np.inf]))


def test_three_hour_cyclic_storage_and_shift():
    energy = np.array([[5.0], [8.0]])
    power = np.array([[1.0], [-1.0]])
    hours = np.array([3.0, 3.0])
    loss = np.zeros_like(energy)
    cyclic, initial = np.array([True]), np.array([0.0])
    np.testing.assert_allclose(
        storage_residual(energy, power, hours, loss, cyclic, initial), 0
    )
    np.testing.assert_allclose(
        storage_residual(energy - 5, power, hours, loss, cyclic, initial), 0
    )
    wrong_hours = storage_residual(energy, power, np.ones(2), loss, cyclic, initial)
    np.testing.assert_allclose(wrong_hours[:, 0], [-2, 2])


@pytest.mark.parametrize(
    "restriction", [None, "e_set", "e_nom_set", "marginal_cost_storage"]
)
def test_coal_adjustment_requires_free_inventory_and_capacity(restriction):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035-01-01", periods=2, freq="3h"))
    n.add("Bus", "coal")
    n.add(
        "Store",
        "EU coal Store",
        bus="coal",
        carrier="coal",
        e_nom_extendable=True,
        e_cyclic=True,
        capital_cost=0.1,
    )
    n.stores.loc["EU coal Store", "e_nom_opt"] = 20
    n.stores_t.e = pd.DataFrame({"EU coal Store": [10.0, 13.0]}, index=n.snapshots)
    if restriction is None:
        coal = audit.coal_audit(n)
        assert coal["same_dispatch_shift_avoidable_cost_eur"] == pytest.approx(1.7)
    else:
        n.stores.loc["EU coal Store", restriction] = 0.1
        with pytest.raises(ValueError, match="Coal shift conditions not established"):
            audit.coal_audit(n)


def test_noncyclic_initial_energy_is_not_multiplied_by_retention():
    # Native PyPSA uses initial energy directly as the first-snapshot RHS.
    r = storage_residual(
        np.array([[7.0], [3.5]]),
        np.array([[1.0], [0.0]]),
        np.ones(2),
        np.full((2, 1), 0.5),
        np.array([False]),
        np.array([8.0]),
    )
    np.testing.assert_allclose(r, 0)


def test_objective_subtracts_nominal_not_minimum_and_excludes_fixed_capex():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=2, freq="3h"))
    n.snapshot_weightings.loc[:, :] = 3
    n.add("Bus", "a")
    n.add(
        "Generator",
        "new",
        bus="a",
        p_nom=2,
        p_nom_min=7,
        p_nom_extendable=True,
        capital_cost=10,
        marginal_cost=3,
    )
    n.add("Generator", "fixed", bus="a", p_nom=5, capital_cost=100, marginal_cost=2)
    n.add(
        "Store",
        "minimum_is_not_an_offset",
        bus="a",
        e_nom=0,
        e_nom_min=100,
        e_nom_extendable=True,
        capital_cost=0.1,
    )
    n.generators.loc[:, "p_nom_opt"] = [8, 5]
    n.stores.loc[:, "e_nom_opt"] = 100
    n.generators_t.p = pd.DataFrame(
        {"new": [1.0, 2.0], "fixed": [4.0, 4.0]}, index=n.snapshots
    )
    n._objective_constant, n._objective = 20.0, 145.0
    _, result = Audit(n).costs()
    assert result["gross_investment_eur"] == pytest.approx(90)
    assert result["existing_capacity_constant_eur"] == pytest.approx(20)
    assert result["operating_eur"] == pytest.approx(75)
    assert result["reconstructed_objective_eur"] == pytest.approx(145)


def test_nodal_signs_and_lossy_conversion_detect_wrong_delivery():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=2, freq="3h"))
    n.add("Bus", "a")
    n.add("Bus", "b")
    n.add("Generator", "supply", bus="a", p_nom=10)
    n.add("Load", "demand", bus="b", p_set=5)
    n.add("Link", "converter", bus0="a", bus1="b", p_nom=10, efficiency=0.5)
    n.generators_t.p = pd.DataFrame(10.0, index=n.snapshots, columns=["supply"])
    n.loads_t.p = pd.DataFrame(5.0, index=n.snapshots, columns=["demand"])
    n.links_t.p0 = pd.DataFrame(10.0, index=n.snapshots, columns=["converter"])
    n.links_t.p1 = pd.DataFrame(-5.0, index=n.snapshots, columns=["converter"])
    a = Audit(n)
    buses, _ = a.balance()
    assert buses.max_balance_error.max() == 0
    n.links_t.p1.loc[n.snapshots[1], "converter"] = -4.75
    buses, _ = Audit(n).balance()
    assert buses.set_index("bus").loc["b", "max_balance_error"] == pytest.approx(0.25)


def test_ac_cycle_reports_native_scaled_equation_units():
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035-01-01", periods=2, freq="3h"))
    for name in ["a", "b", "c"]:
        n.add("Bus", name, v_nom=1, carrier="AC")
    for name, source, target in [("ab", "a", "b"), ("bc", "b", "c"), ("ca", "c", "a")]:
        n.add("Line", name, bus0=source, bus1=target, x=1, r=0, s_nom=10)
    n.lines_t.p0 = pd.DataFrame(1.0, index=n.snapshots, columns=n.lines.index)
    a = Audit(n)
    a.kvl()
    assert a.metrics[0]["max_violation"] == pytest.approx(3e5)
    assert "1e5 scaling" in a.metrics[0]["unit"]
    n.lines_t.p0.loc[:, ["ab", "bc", "ca"]] = [1.0, 0.0, -1.0]
    a = Audit(n)
    a.kvl()
    assert a.metrics[0]["max_violation"] == pytest.approx(0)


@pytest.mark.parametrize("layout", ["candidate", "accepted", "original"])
def test_full_saved_audit_preserves_candidate_and_finds_deliberate_imbalance(
    tmp_path, monkeypatch, layout
):
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2035-01-01", periods=2, freq="3h"))
    n.snapshot_weightings.loc[:, :] = 3
    n.add("Carrier", "AC")
    n.add("Carrier", "coal")
    n.add("Bus", "a", carrier="AC")
    n.add("Bus", "fuel", carrier="coal")
    n.add(
        "Generator",
        "supply",
        bus="a",
        p_nom_extendable=True,
        capital_cost=10,
        marginal_cost=4,
    )
    n.generators.loc["supply", "p_nom_opt"] = 5
    n.add("Load", "demand", bus="a", p_set=5)
    n.add(
        "Store",
        "EU coal Store",
        bus="fuel",
        carrier="coal",
        e_nom_extendable=True,
        e_cyclic=True,
        capital_cost=0.1,
    )
    n.stores.loc["EU coal Store", "e_nom_opt"] = 20
    n.add("StorageUnit", "empty", bus="a", p_nom=0)
    n.storage_units.loc["empty", "p_nom_opt"] = 0
    n.add("Link", "unused", bus0="a", bus1="fuel", p_nom=0)
    n.links.loc["unused", "p_nom_opt"] = 0
    n.links["reversed"] = False
    n.links["energy to power ratio"] = 0.0
    for component, attr, name, values in [
        ("Generator", "p", "supply", [5.0, 5.0]),
        ("Load", "p", "demand", [5.0, 5.0]),
        ("Store", "p", "EU coal Store", [0.0, 0.0]),
        ("Store", "e", "EU coal Store", [10.0, 10.0]),
        *[
            ("StorageUnit", attr, "empty", [0.0, 0.0])
            for attr in ["p", "p_dispatch", "p_store", "state_of_charge", "spill"]
        ],
        *[("Link", attr, "unused", [0.0, 0.0]) for attr in ["p", "p0", "p1"]],
    ]:
        n.components[component].dynamic[attr] = pd.DataFrame(
            {name: values}, index=n.snapshots
        )
    # A known 0.25 MW imbalance must be located in the saved physical check.
    n.loads_t.p.loc[n.snapshots[1], "demand"] = 4.75
    n.buses_t.marginal_price = pd.DataFrame(
        {"a": [4.0, 6.0], "fuel": [0.0, 0.0]}, index=n.snapshots
    )
    n._objective, n._objective_constant = 172.0, 0.0
    root = tmp_path / "candidate"
    root.mkdir()
    metadata = {"status": "failed", "scenario": "test", "screening_issues": ["keep"]}
    (root / "summary.json").write_text(json.dumps(metadata))
    config = {
        "solving": {"options": {"transmission_losses": 0}},
        "renewable": {
            "solar": {"capacity_per_sqkm": 1},
            "solar-hsat": {"capacity_per_sqkm": 1},
        },
    }
    if layout == "original":
        network = root / "networks" / "solved_2035.nc"
        n.meta = config
    else:
        network = root / (
            "candidate_2035.nc" if layout == "candidate" else "solved_2035.nc"
        )
        (root / "effective_config.yaml").write_text(yaml.safe_dump(config))
    network.parent.mkdir(parents=True, exist_ok=True)
    n.export_to_netcdf(network)
    before = {
        str(p.relative_to(root)): audit.sha256(p)
        for p in root.rglob("*")
        if p.is_file()
    }

    def forbidden(*args, **kwargs):
        raise AssertionError("The saved audit must not build or solve an LP")

    from pypsa.optimization.optimize import OptimizationAccessor

    monkeypatch.setattr(OptimizationAccessor, "__call__", forbidden)
    monkeypatch.setattr(OptimizationAccessor, "create_model", forbidden)
    output = tmp_path / "audit"
    report = audit.run(root, output)
    assert report["optimization_called"] is False
    assert report["input_network"] == str(network)
    assert (report["effective_config_source"] == "saved_network_meta") == (
        layout == "original"
    )
    assert report["original_run_status"] == "failed"
    assert report["original_screening_issues"] == ["keep"]
    assert report["objective"]["reconstruction_error_eur"] == pytest.approx(0)
    assert report["coal_shifted_candidate_objective_eur"] == pytest.approx(170)
    checks = pd.read_csv(output / "constraint_checks.csv")
    balance = checks.loc[checks["check"].eq("nodal_balance")].iloc[0]
    assert balance.max_violation == pytest.approx(0.25)
    prices = pd.read_csv(output / "electricity_prices.csv.gz", index_col=0)
    assert prices["a"].tolist() == [4.0, 6.0]
    assert before == {
        str(p.relative_to(root)): audit.sha256(p)
        for p in root.rglob("*")
        if p.is_file()
    }
    with pytest.raises(FileExistsError):
        audit.run(root, output)


@pytest.mark.parametrize("task", range(8))
@pytest.mark.parametrize("campaign", ["4930796", "original"])
def test_saved_audit_array_routes_every_scenario_without_a_solve(
    tmp_path, task, campaign
):
    cases = list(yaml.safe_load(Path("config/scenarios_bornholm.yaml").read_text()))
    wrapper = tmp_path / "pixi-sophia"
    wrapper.write_text(
        "#!/bin/bash\npython3 -c 'import sys,json; print(json.dumps(sys.argv[1:]))' \"$@\"\n"
    )
    wrapper.chmod(0o755)
    srun = tmp_path / "srun"
    srun.write_text(
        '#!/bin/bash\n[[ "$1" == --ntasks=1 && "$2" == --cpus-per-task=1 && "$3" == --cpu-bind=cores ]] || exit 3\nshift 3\nexec "$@"\n'
    )
    srun.chmod(0o755)
    script = Path("scripts/bornholm_saved_audit_workq.sbatch").resolve()
    result = subprocess.run(
        ["bash", str(script), campaign],
        env=os.environ
        | {
            "SLURM_SUBMIT_DIR": str(tmp_path),
            "SLURM_JOB_PARTITION": "workq",
            "SLURM_CPUS_PER_TASK": "1",
            "SLURM_ARRAY_TASK_ID": str(task),
            "SLURM_ARRAY_JOB_ID": "555",
            "PIXI_SOPHIA": str(wrapper),
            "PATH": str(tmp_path) + os.pathsep + os.environ["PATH"],
        },
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    args = json.loads(result.stdout)
    assert args[:5] == [
        "run",
        "--as-is",
        "python",
        "-m",
        "scripts.audit_bornholm_candidate",
    ]
    input_root = (
        "results/bornholm_3h_s1_w2_o3"
        if campaign == "original"
        else "results/bornholm_accuracy/4930796"
    )
    assert args[5] == input_root + "/" + cases[task]
    assert args[args.index("--output") + 1] == (
        f"results/bornholm_saved_audit/{campaign}/555/" + cases[task]
    )
