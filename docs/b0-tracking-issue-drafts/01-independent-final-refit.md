# B0: add independent final weak-prior refit after CKF finding

## Motivation

The current B0 chain uses `B0TrackerStubSeeder` to fit B0 measurements and propagate those fitted uncertainties into the seed covariance. `CKFTracking` then uses that seed covariance as its initial prior while processing the B0 measurements again.

That is useful for **track finding**, but it is not a clean statistical definition of the **final fitted covariance** because information from the same measurements can enter twice: once through the hit-derived seed prior and again through the Kalman measurement updates.

Relevant code:

- `src/algorithms/tracking/B0TrackerStubSeeder.cc`
- `src/algorithms/tracking/B0TrackerStubSeederConfig.h`
- `src/algorithms/tracking/CKFTracking.cc`
- `src/global/tracking/tracking.cc`
- `src/benchmarks/reconstruction/b0_tracking/README.md`

The benchmark README already anticipates a follow-up named `doFinalWeakPriorRefit`; this issue should turn that idea into a well-defined implementation and validation task.

Related upstream discussion: https://github.com/eic/EICrecon/issues/2930

## Goal

Use the current hit-derived B0 seed for **candidate finding**, but compute the final track parameters and covariance from an **independent refit of the selected measurements**, so each measurement contributes only once to the reported fit uncertainty.

Central tracking must remain unchanged by default.

## Proposed design

Add an explicit opt-in capability to `CKFTrackingConfig`, for example:

```cpp
bool doFinalWeakPriorRefit = false;
```

Do **not** infer this behavior from `numB0StationsMin`; the switch should be explicit and default off so central tracking is unaffected.

For the B0 CKF factory instances only:

```text
stub seed
  -> CKF track finding / hit selection
  -> selected measurement set
  -> independent final refit with deliberately weak initialization
  -> smoothing
  -> extrapolation to requested output reference surface
```

The first implementation can use the established ACTS Kalman fitter. A global chi2 fitter can be investigated later, but should not block this issue.

### Initialization requirements

The final refit should use the seed parameters only as a numerical starting point. Its prior covariance must be weak enough that the final result is measurement-dominated.

A useful validation is to scale the refit prior covariance over several orders of magnitude and demonstrate that final parameters/covariances are stable once the prior is sufficiently weak.

A real external beam-spot constraint, if intentionally used, should be applied once and documented separately from the weak numerical initialization.

## Acceptance criteria

- [ ] Add a shared, explicit, default-off final-refit configuration option.
- [ ] Enable it only for the B0 reconstruction chain initially.
- [ ] Refitted tracks preserve the CKF-selected measurement set.
- [ ] Final output parameters/covariance come from the independent refit, not from the hit-derived seed prior.
- [ ] Central tracking is unchanged when the option is off.
- [ ] B0 reconstruction efficiency is not significantly degraded relative to the pre-refit baseline.
- [ ] Final residuals and pull widths are added to the B0 benchmark.
- [ ] Demonstrate stability under substantial weakening of the refit prior.
- [ ] Document and count failed refits; either reject explicitly or use a clearly marked fallback.

## Validation

At minimum run the pinned B0 benchmark on:

1. prompt proton gun sample,
2. 3-station and 4-station reconstructible subsets,
3. both charge signs / a charged pion sample,
4. a realistic-background or occupancy sample when available.

Compare before/after:

- total efficiency,
- `epsilon_track|seed`,
- fake and duplicate rates,
- momentum / angle residuals,
- `chi2/ndf`,
- 1D pulls and multivariate Mahalanobis `D^2`,
- runtime.

The pull and `D^2` plots should become physics-quality diagnostics only after this issue is complete.

## Non-goals

- Replacing B0 pattern recognition.
- Changing central tracking behavior.
- Replacing the existing stub seed estimator.
- Tuning the B0 chi2 window before the final-fit statistics are trustworthy.
