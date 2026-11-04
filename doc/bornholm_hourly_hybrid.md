# Reviewed hourly fixed Hybrid

The hourly option uses the reviewed common base, three onshore resource classes,
2035 assumptions and 8,760 hourly snapshots of 2013 weather. It solves only the
fixed DK2–DE Hybrid project. The main eight-case campaign remains three-hourly.

```bash
sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch 16 standard
```

The launcher reserves one exclusive fatq node with 32 CPUs and uses 16 solver
threads by default. The input prefix is `bornholm_1h_s1_w3_o3_reviewed`. Its
processed-input cache is separate from the main campaign; both the profile
outputs and their logs follow the same scenario wildcards.

The driver measures available RAM after model construction, applies its memory
reserve, limits solver runtime and retains failed candidates and solver logs.
`retry` and `aggfill0` require a previous job ID and a matching source hash;
old two-class hourly attempts cannot serve as a reference for this new input.
The same physical and numerical acceptance checks as the reviewed three-hour
campaign apply, including transmission, coal and chronological storage checks.

See [the common-base instructions](bornholm_calibrated_base.md) before building.
Earlier performance investigations are preserved in the
[hourly experiment history](../archive/bornholm_solver_experiments/doc/bornholm_hourly_hybrid.md).
A successful full-year hourly solve is not claimed by this implementation update.
