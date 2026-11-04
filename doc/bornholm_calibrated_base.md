# Bornholm: reviewed common base

The current campaign uses the existing `feature/bornholm-energy-island` branch,
three onshore resource classes, 2035 system assumptions and 2013 weather at
three-hour resolution. All eight scenarios retain the same 34-country footprint,
background grid, resource assumptions and project-area reservation. Eleven
countries excluded from manual onshore CF adjustments remain in the model.

The previous resource/coal correction is integrated into this branch. Denmark's
onshore multiplier is still a **provisional 2035 scenario assumption of 0.60**;
it is not the separately estimated historical DK1/DK2 correction. Switzerland's
4.00 factor remains a provisional spatial-resource proxy requiring terrain,
site representation and existing-fleet checks. Other country factors are in
`data/model_assumptions/onshore_calibration.csv`. All adjustments act on hourly
weather profiles before averaging and clipping; they do not create constant
output or prescribe an optimised technology mix. Native offshore losses are not
applied a second time. The initial assumptions and source discussion are retained
[in the archive](../archive/bornholm_solver_experiments/doc/bornholm_calibrated_base_initial.md).

The fixed 1,000 MWh coal Store is a numerical fuel buffer with an independent
fuel supplier. It does not cap annual coal consumption. It applies only to the
sector-coupled Bornholm case, never to the electricity-only hindcast.

The [reviewed transmission specification](transmission_assumptions.md) supersedes
the earlier four-corridor implementation. In particular, SE2–SE3 uses the cited
10,350 MW 2035 step in both directions. Core fallback NTCs are explicit scenario
approximations; the flow-based market remains incomplete.

## Run on Sophia

Pull this branch, then build new networks. A re-solve of old composed networks
will not install the corrections. Fuel, weather and geographic downloads can
be reused when their workflow dependencies still match.

```bash
git pull --ff-only
sbatch scripts/bornholm_calibrated_workq.sbatch dry-run
sbatch scripts/bornholm_calibrated_workq.sbatch build
```

After reviewing the updated 2023 hindcast, start the baseline and Hybrid first:

```bash
sbatch --array=0,2 scripts/bornholm_calibrated_solve.sbatch
```

The full campaign can be submitted with:

```bash
sbatch scripts/bornholm_calibrated_workq.sbatch all
```

It builds and preflights all eight cases, then submits the solve array and a
comparison check after successful completion. The old study workq command is
an alias for this reviewed controller. The study module CLI also selects the
reviewed campaign; `--legacy` is reserved for archived reproductions. The output prefix is
`bornholm_3h_s1_w3_o3_reviewed`; old results cannot satisfy these new targets.
Do not use a new campaign's comparison checker until all eight results exist.

The hourly fixed Hybrid is retained with the same reviewed assumptions:

```bash
sbatch scripts/bornholm_hourly_hybrid_fatq.sbatch 16 standard
```

It uses `bornholm_1h_s1_w3_o3_reviewed` and its own temporal cache. An older
hourly retry summary cannot be used for the reviewed input, because source
hashes must match. Hourly feasibility and performance still require a full run.

## Acceptance and interpretation

Accepted networks require an optimal solver status and the existing barrier
quality screen, plus physical bus balances, generator availability, chronological
storage balances, the fixed coal buffer and all installed directional/joint
transfer limits. Failed attempts retain diagnostic files and cannot overwrite
an accepted output. Metadata records the input hashes and assumptions.

Compare project effects against the matching baseline within each grid regime.
An absolute `vfixed` versus `vopt` objective difference still includes accounting
for the existing transmission-capacity floor and is not by itself a welfare
estimate. The source-level limitations in the transmission specification and
the provisional CF factors remain relevant when interpreting results.

No full annual HPC solve is claimed by this implementation update.


## Rechecking completed inputs after a timestamp-unit validation failure

The transmission snapshot fingerprint is now normalized to nanoseconds. Earlier
version-1 metadata hashed raw datetime integers; a microsecond index could be
reloaded from NetCDF as nanoseconds and spuriously trigger “Transmission
snapshots changed; recompose”. Existing metadata is accepted only if a lossless
representation of the same ordered timestamps reproduces its digest. Real
changes to dates, ordering or count still fail. Physical limits and network
inputs are not altered by this compatibility fix.

If composition already reached 119/119 and failed only during validation, run
checks on a workq node without rebuilding the inputs, then release just the
baseline and Hybrid solves after success:

```bash
cd ~/projects/pypsa-eur-bornholm
git pull --ff-only &&
check_job=$(sbatch --parsable scripts/bornholm_calibrated_workq.sbatch check) &&
sbatch --array=0,2 --dependency="afterok:${check_job%%;*}" \
  --kill-on-invalid-dep=yes scripts/bornholm_calibrated_solve.sbatch
```

The `check` stage reads all eight existing composed networks and runs their
comparison and solve preflights. It neither launches composition workers nor
submits solves itself. A genuinely changed network will still fail validation;
do not edit its fingerprint manually to bypass that failure.
