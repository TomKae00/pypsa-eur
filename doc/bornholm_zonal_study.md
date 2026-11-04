# Bornholm study

The active workflow is described in [the reviewed common base](bornholm_calibrated_base.md).
It covers the existing eight scenarios with 2035 assumptions, 2013 weather,
three-hour main runs and an hourly fixed Hybrid option.

Use `sbatch scripts/bornholm_calibrated_workq.sbatch` for the reviewed campaign.
`scripts/bornholm_study_workq.sbatch` is a compatibility alias to that controller.

Earlier construction notes, solver trials and result namespaces are preserved
[in the campaign archive](../archive/bornholm_solver_experiments/doc/bornholm_zonal_study.md).
