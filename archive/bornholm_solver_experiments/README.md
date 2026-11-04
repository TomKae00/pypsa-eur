# Superseded Bornholm experiments

This folder preserves the earlier solver sweeps, simplified configurations and
campaign history. Their commands remain relative to the **repository root**.
They are retained for reproducing and interpreting old runs, not for starting
the reviewed campaign. Existing result directories are not deleted.

The active three-hour campaign is `scripts/bornholm_calibrated_workq.sbatch`.
The active hourly Hybrid launcher is `scripts/bornholm_hourly_hybrid_fatq.sbatch`.
The familiar `scripts/bornholm_study_workq.sbatch` now forwards to the reviewed
campaign. `scripts/run_bornholm_study.py` retains legacy helper defaults because
the archived accuracy tools import them; use `scripts.run_bornholm_calibrated`
for current targets and solves. The study module CLI also forwards to that
reviewed campaign; `--legacy` explicitly selects its earlier command behavior.

Reusable accuracy, memory and provenance Python helpers remain in `scripts/`.
Their earlier experiment entry points and configurations are archived here,
with references updated in documentation and regression tests. No result or
input network is removed, rewritten or presented as a reviewed result.
