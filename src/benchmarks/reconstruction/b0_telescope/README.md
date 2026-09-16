# B0 telescope comparison and prior validation

This benchmark consumes the opt-in `b0_telescope` plugin from PR #6, not the legacy origin-perigee outputs. It is executable analysis and regression code; it does not contain invented detector performance numbers.

## Run

Inside a matching eic-shell with the built plugin installed:

```sh
python src/benchmarks/reconstruction/b0_telescope/run_b0_telescope.py \
  --sim /absolute/sim.edm4hep.root \
  --compact /absolute/epic_ip6_extended.xml \
  --calibration /absolute/calibrations.xml \
  --material /absolute/matching-material-map.cbor \
  --out results/b0-scale1 --weak-prior-scale 1
```

The output directory must not exist. Inputs are explicit; `--dry-run` records the exact argument vector and manifest without claiming reconstruction happened. `--chain` selects one or more of the four chains; omitted means all four. For a cut/field/timing scan, pass e.g. `--set B0TelescopeSeeds:useTiming=true` or `--set B0TelescopeCKF:chi2Cut=25`. `--pdg 211 --pdg -211` selects pions in the analysis, not in reconstruction. Configure fit and timing mass hypotheses separately for controlled comparisons.

The driver writes a reconstruction log, the ROOT output, station manifest, local truth reference, diagnostics, and content-addressed provenance. It hashes explicit inputs, recursively referenced compact XML, the executable and extra `--asset` files before/after reconstruction. Pass loaded detector/plugin libraries and auxiliary calibration files as additional assets. `--source` records actual checkout SHA/dirty state but does not claim those sources produced the installed binary. Analysis verifies the ROOT/truth hashes when supplied a manifest.

No shell command is built from an unquoted concatenated string; subprocess receives an argument list. Failed runs remain marked failed and are not accepted as completed provenance. Absolute include references are supported, but unresolved environment references and remote XML includes fail rather than silently produce incomplete provenance.

## Definitions

A target particle is a supported charged species with actual B0 crossings in at least three physical stations, optionally restricted by signed PDG and production momentum. Generator-status and production-z restrictions from central truth seeding are not inherited. Unsupported species and missing/invalid local truth states are reported separately. Event IDs are local to a run; duplicate IDs fail rather than silently merge unrelated events.

Track/seed matching uses full `(collectionID,index)` identities and strict majority purity. Unknown contributions stay in denominators. A confidently matched real particle outside the selected target sample is counted as **matched non-target**, not fake. Ambiguous and unassociated tracks are separate fake categories. Duplicates are extra matched tracks beyond one per target particle. Representative tracks are chosen by reconstructed measurement count and fit quality, never by distance to truth.

Seed efficiency, reconstruction given a matched seed, and total efficiency are independent counters. CKF can recover a particle from a seed matched to another particle; this is explicitly counted, so the efficiency product is not imposed as an identity. Rates include 95% Wilson intervals and undefined denominators remain null.

Residuals use truth and fitted parameters on the **same actual sensor surface**. Phi is wrapped, full spatial covariance is checked for symmetry and positive definiteness, and D2 is computed with a Cholesky solve. An invertible indefinite covariance is not accepted. Exact-zero curvature does not become a fabricated zero-momentum residual. Missing reference states never remove a reconstructed particle from the efficiency numerator.

Reports: `analysis/report.md`, `summary.json`, and `matched_tracks.jsonl`. The JSON contains momentum/vertex/station/occupancy bins, seed/fit failures, charge errors, residual/pull quantiles and five-dimensional D2 diagnostics. Per-chain elapsed time includes the tracking/refit algorithm but excludes seeding, digitization, ambiguity resolution and writing; use reconstruction wall time and component diagnostics together. Running all four chains is not a single-chain timing benchmark.

## Independent-prior scan

Run the same sample at fresh-prior scales 1, 10 and 100 into separate directories, then:

```sh
python src/benchmarks/reconstruction/b0_telescope/compare_priors.py \
  results/b0-scale1/analysis results/b0-scale10/analysis results/b0-scale100/analysis \
  --out results/prior-comparison.json
```

The comparator verifies common sample/geometry/material hashes and target selections. It compares means and normalized covariance changes only for identical selected measurement IDs on identical reference surfaces; changed associations and missing tracks are reported separately. Thresholds are explicit diagnostic tolerances, not universal detector requirements. Passing this scan does not establish correct material, truth binding or model calibration.

## Required physics matrix

Use prompt protons with vertex spread; both pion charges; displaced daughters; mixed species; three-station recovery; close pairs; front/back single-sided coverage; module/aperture edges; and background overlays. Compare truth-associated fitting, truth-seeded CKF, realistic CKF/refit and realistic direct fitting on identical events. Repeat for supported field/geometry configurations. Beam-spot-constrained fitting, calibrated hardware response and generator/sample production are not supplied by this analysis script.

## Tests

```sh
cd src/benchmarks/reconstruction/b0_telescope
python -m unittest -v
```

The pure tests require NumPy, not ROOT/PODIO. They cover matching/denominator correctness, ties and collection IDs, wrapped residuals, indefinite covariance, missing references, duplicate events, provenance tampering, recursive includes, safe command construction and same-hit prior comparisons. The actual PODIO adapter/driver require detector-backed smoke testing in eic-shell.
