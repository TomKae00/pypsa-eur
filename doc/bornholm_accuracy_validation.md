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
