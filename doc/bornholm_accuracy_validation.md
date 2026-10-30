# Paired annual accuracy validation on workq

Replay the **same completed three-hour baseline** twice, without rebuilding,
changing technologies, or overwriting the original eight study solutions.

| Array task | BarConvTol | NumericFocus | Threads |
|---|---:|---:|---:|
| 0 | 1e-8 | 0 | 8 |
| 1 | 1e-8 | 2 | 8 |

Both jobs use barrier (Method=2), crossover off, the existing 20-hour solver limit,
a 24-hour walltime, and the existing memory ceiling capped at 75% of available
physical RAM. Each gets its own exclusive **workq** node. They can run concurrently
when the scheduler has resources; neither depends on the other.

The completed campaign's Slurm peaks were approximately 37–40 GiB on nodes with
about 125 GiB physical RAM. Eight solver threads is a moderate next trial with
useful headroom. Barrier memory can increase with threads, and doubling threads
does not guarantee a twofold speedup. This pair isolates NumericFocus at eight
threads; comparison with the older four-thread run also changes tolerance, so it
is not a clean thread-scaling benchmark. Node hardware can affect runtime as well.

## Submit

From the cluster repository root after pulling this change:

```bash
pair_job=$(sbatch --parsable scripts/bornholm_accuracy_workq.sbatch)
pair_job=${pair_job%%;*}
echo "Accuracy comparison array: $pair_job"
watch -n 10 "squeue -j $pair_job -o '%.18i %.10P %.12T %.12M %.25R'"
```

This submits exactly **two solves**, not the eight-scenario build/solve campaign.
The script refuses other partitions and allocations that do not supply eight CPUs.

## Follow the solves

```bash
tail -n 25 -F "results/bornholm_accuracy/$pair_job/tol8_nf0_t8/solver.log"
```

Use `tol8_nf2_t8` for NumericFocus=2. Before a solver log exists, follow startup:

```bash
tail -n 40 -F "slurm-bornholm-accuracy-${pair_job}_0.out"
tail -n 40 -F "slurm-bornholm-accuracy-${pair_job}_1.out"
```

Memory and final accounting:

```bash
sstat -j "${pair_job}_0.batch,${pair_job}_1.batch" --format=JobID,MaxRSS,AveCPU -P
sacct -j "$pair_job" --format=JobID,State,Elapsed,MaxRSS,ExitCode -P
```

## Results and screening

Each run records the input SHA256, code provenance, effective solver options,
elapsed time, process peak RSS and Gurobi peak memory. Gurobi quality attributes
include original and scaled primal/dual violations, complementarity and iteration
count. Unavailable attributes are recorded explicitly.

The source hash must match the original completed baseline summary. Critical solve
code and configuration hashes must match the completed campaign; the separately
reviewed no-op custom hook is pinned explicitly. This prevents an unannounced
rebuild or a different physical formulation from entering the comparison.

Diagnostic screens go beyond the solver's `optimal` status:

* Avoidable coal-storage cost from shifting inventory and reducing capacity must be
  at most EUR 1 million/year. No inventory or capacity is actually modified.
* Unscaled BoundVio, ConstrVio and ConstrResidual must be at most 1e-3 in their
  respective model units; DualVio and DualResidual at most 1e-5. These broad
  screens do not substitute for carrier-specific physical checks.
* The final displayed primal/dual separation must be at most EUR 5 million.
  The logged dual iterate is **not a certified lower bound**, so this is a
  diagnostic screen, not a rigorous error bound on system cost or Bornholm value.

An accepted run reports `status: checks_passed` and renames its candidate to
`solved_2035.nc`. A rejected candidate remains `candidate_2035.nc` for investigation,
with `status: failed`, screening reasons and a nonzero job exit. Failure of one
array task does not cancel the other. Even a passing run still requires economic
review; no independent optimality certificate is claimed.

Compare the two summaries after both jobs finish:

```bash
~/bin/pixi-sophia run python -m scripts.run_bornholm_accuracy \
  --compare "results/bornholm_accuracy/$pair_job"
```

The comparison checks matching input hashes and common solver settings, then
prints costs, the coal artifact, runtime, memory, iterations and numerical quality
side by side as JSON. Reusing an existing validation output directory is refused.


## Four-profile NumericFocus=1 follow-up

After the NF0 replay returned suboptimally and the NF2 run proved slow, submit a
separate four-task array. This uses the **same original three-hour baseline** as
before, not four different Bornholm scenarios. The old two-task launcher remains
available. No original solve/configuration, network or existing result is changed.

| Task | Profile | ScaleFlag | Crossover | BarHomogeneous | Result folder |
|---|---|---:|---:|---:|---|
| 0 | `nf1` | -1 (automatic) | 0 (off) | -1 (automatic) | `tol8_nf1_t8` |
| 1 | `nf1_scale2` | 2 | 0 | -1 | `tol8_nf1_scale2_t8` |
| 2 | `nf1_crossover` | -1 | -1 (automatic) | -1 | `tol8_nf1_crossover_t8` |
| 3 | `nf1_homogeneous` | -1 | 0 | 1 (on) | `tol8_nf1_homogeneous_t8` |

All four use **NumericFocus=1, BarConvTol=1e-8 and eight threads**. Task 3 now
uses homogeneous barrier instead of the previously proposed 1e-7 tolerance test.
This tests convergence behaviour without relaxing the stopping tolerance. Existing
results from the earlier trial remain untouched and can still be compared.

Each alternative changes exactly one setting relative to task 0. All retain
Method=2, eight threads, 20 hours of solver time, 24 hours of walltime, the inherited
memory ceiling capped at 75% of available RAM, and one exclusive workq node.
The scheduler can run all four concurrently if resources are available; otherwise
they queue independently. There are no dependencies on the previous NF0/NF2 jobs.
Pulling this change does not cancel or restart those jobs.

These are experiments, not an assurance of convergence or faster execution.
Crossover, homogeneous barrier and scaling can alter runtime/memory and must be
evaluated using final unscaled quality metrics. Homogeneous barrier does not imply
that the physical model is infeasible or unbounded.

All trials require native Gurobi status OPTIMAL (2), PyPSA's optimal termination,
finite coal-cost diagnostics with at most EUR 1m/year of avoidable cost, and the
original unscaled bound/constraint/dual violation and residual screens. Missing
or non-finite required metrics fail screening.

When crossover is **off**, the existing EUR 5m displayed-barrier-separation check
remains an additional diagnostic screen. The last row need not describe the
returned candidate, and the logged dual iterate is not a certified lower bound.
When crossover is **enabled**, that row can precede simplex cleanup and is kept
for information only. Acceptance instead uses the returned solution's native
status and original-unit quality metrics, including finite `ComplVioSum` at most
EUR 5m in this model's objective units. This complementarity screen is also a
diagnostic, not a certified objective-error bound. Poor final residuals, a
suboptimal return, or the coal artifact still cause rejection.

Summaries record `screening_basis`; comparisons show final `ComplVio`,
`ComplVioSum` and the earlier displayed separation separately. Neither passing
these screens nor similar system costs alone establishes reliable welfare or
price differences.

Submit **four new solves**, without rebuilding:

```bash
cd ~/projects/pypsa-eur-bornholm
git pull --ff-only origin feature/bornholm-energy-island
sweep_job=$(sbatch --parsable scripts/bornholm_solver_sweep_workq.sbatch)
sweep_job=${sweep_job%%;*}
echo "Four solver trials: $sweep_job"
watch -n 10 "squeue -j $sweep_job -o '%.18i %.10P %.12T %.12M %.25R'"
```

Follow task 0 (substitute another result-folder label for the others):

```bash
tail -n 30 -F "results/bornholm_accuracy/$sweep_job/tol8_nf1_t8/solver.log"
```

Startup, memory and final state:

```bash
tail -n 40 -F "slurm-bornholm-solver-sweep-${sweep_job}_0.out"
sstat -j "${sweep_job}_0.batch,${sweep_job}_1.batch,${sweep_job}_2.batch,${sweep_job}_3.batch" \
  --format=JobID,MaxRSS,AveCPU -P
sacct -j "$sweep_job" --format=JobID,State,Elapsed,MaxRSS,ExitCode -P
```

After the runs finish:

```bash
~/bin/pixi-sophia run python -m scripts.run_bornholm_accuracy \
  --compare "results/bornholm_accuracy/$sweep_job"
```

`--compare` now reads all immediate child `summary.json` files, so it works for
both the old pair and the four-profile array. Its JSON contains a `runs` mapping
keyed by result-folder label, all effective solver options, failure reasons,
quality and timing. The source hashes must agree; solver options outside the
experiment's settings and the adaptive memory ceiling must agree. Missing run
folders are not invented or treated as completed trials; check Slurm for all
submitted tasks before making a final comparison.

If the native solve script raises on a non-optimal return before exporting a
candidate, the wrapper now preserves accessible Gurobi status/quality and the last
barrier row in the failed summary. A NetCDF is not guaranteed in that case. Such a
run remains failed, even if Gurobi exposes a candidate objective. These diagnostics
do not bypass the native rejection or change the physical formulation.

For one trial only, e.g. the NF1 control:

```bash
sbatch --array=0 scripts/bornholm_solver_sweep_workq.sbatch
```

Do not edit the physics or rebuild between trials. No objective-constant removal,
coal-capacity adjustment, balancing slack, technology removal or model rescaling
is introduced by this change.

## Four additional trials: threads and homogeneous barrier

Keep the existing array **4929796** running. Submit the new follow-up launcher
once to add four independent solves of the **same composed three-hour baseline**:

| New task | Profile | Threads | NumericFocus | BarHomogeneous | Crossover | Existing comparison |
|---|---|---:|---:|---:|---:|---|
| 0 | `nf1_t16` | 16 | 1 | -1 | 0 | NF1 control, eight threads |
| 1 | `nf1_homogeneous_t16` | 16 | 1 | 1 | 0 | Homogeneous NF1, eight threads |
| 2 | `nf1_homogeneous_crossover` | 8 | 1 | 1 | -1 | Homogeneous NF1, crossover off |
| 3 | `nf0_homogeneous` | 8 | 0 | 1 | 0 | Homogeneous NF1, crossover off |

Each changes **one solver parameter** from its existing comparison. All keep
BarConvTol=1e-8, Method=2, the same scaling setting (-1), the inherited solver
options, 20-hour solver limit and adaptive memory ceiling. This probes thread
sensitivity, crossover cleanup and whether homogeneous barrier permits a faster
NumericFocus=0 trajectory. Early progress is not evidence of final accuracy.

Each task gets its own exclusive **workq** node with a 24-hour walltime. The new
array reserves **16 CPUs per task** because resource requests are shared by an
array; tasks 2 and 3 deliberately use **eight Gurobi threads** to match their
reference. No other BLAS/OpenMP parallelism is enabled. Node hardware and available
RAM may differ, so timings are not a perfectly isolated hardware benchmark.
The 16-thread runs keep the memory ceiling; it is not automatically increased.

The original four-task launcher and NF0/NF2 pair launcher are unchanged. Pulling
this change does not cancel, restart or modify a running solve. The same input
hash and physical-code checks, output-isolation rules and numerical screens apply
to all new trials. Existing outputs are not overwritten and no rebuild is needed.

Submit only the **new** array:

```bash
cd ~/projects/pypsa-eur-bornholm
git pull --ff-only origin feature/bornholm-energy-island
followup_job=$(sbatch --parsable scripts/bornholm_solver_followup_workq.sbatch)
followup_job=${followup_job%%;*}
echo "Four additional solver trials: $followup_job"
watch -n 10 "squeue -j 4929796,$followup_job -o '%.18i %.10P %.12T %.12M %.25R'"
```

Result folders under `results/bornholm_accuracy/$followup_job/` are:

* `tol8_nf1_t16`
* `tol8_nf1_homogeneous_t16`
* `tol8_nf1_homogeneous_crossover_t8`
* `tol8_nf0_homogeneous_t8`

Inspect the last five lines of all available solver logs:

```bash
tail -n 5 results/bornholm_accuracy/4929796/*/solver.log \
  results/bornholm_accuracy/"$followup_job"/*/solver.log
```

Before those logs exist, inspect startup output:

```bash
tail -n 10 slurm-bornholm-solver-followup-"${followup_job}"_*.out
```

Compare the original four and the new four after they finish:

```bash
~/bin/pixi-sophia run python -m scripts.run_bornholm_accuracy \
  --compare results/bornholm_accuracy/4929796 \
            "results/bornholm_accuracy/$followup_job"
```

Multiple-array comparisons prefix run labels with the array directory name, show
actual solver threads and node names, and permit the planned thread-count
variation. All source hashes must agree and solver options outside the explicit
experiment settings/adaptive memory ceiling must still match. The comparison
remains read-only and does not imply that a rejected candidate is valid.
