# B0 telescope numerical kernels

These kernels implement the experimental local-reconstruction alternative from the B0 review. They do not alter the legacy stub seeder or central tracking. The ACTS adapter and plugin are supplied by dependent PRs.

## State and field convention

The state is `(x,y,dx/dz,dy/dz,q/p)` at the earliest selected sensor, in the ion frame; positions are mm, field tesla, momentum GeV, and charge magnitude one. A weighted five-parameter fit uses both transverse coordinates, QR factorization with longitudinal scaling, the full 2D measurement covariance, and the projection of tilted-plane uncertainties into transverse residuals.

`fieldIntegral=false` uses a midpoint-field parabolic approximation. The optional integral mode transports a reference trajectory and field-response integrals with RK4 and iterates the frozen-path fit. It is an approximate seed estimator, not a material-aware final fit. It needs convergence/step-size tests against the actual DD4hep field and ACTS transport before promotion.

## Candidate association

Station triplets are extended with at most one measurement per actual sensor surface. Front/back hits remain distinct measurements; they are never replaced by a proximity identity. Single-sided stations and genuine three-station candidates remain available. Identical hit sets are deduplicated; close tracks are not merged solely because their hits are nearby.

`maxTrials` charges endpoint tests, triplet tests and extension tests, not just completed fits. A truncated extension is never silently emitted as complete. The candidate count is also capped. This provides bounded work, not guaranteed high-occupancy efficiency or a spatial-index optimization. Monitor truncation and candidate-drop counters in every performance comparison.

Both transverse roads include an explicit `fieldBoundTesla` envelope. This is a user-supplied bound, not a certified bound on an arbitrary map. Revalidate it when the beamline field changes.

## Charge and time

Exactly zero and insignificant curvature produce finite opposite-charge numerical proposals. Changing q/p also shifts correlated position and slope means using the joint local fit; it does not only flip q/p. This is a proposal construction, not a second independent observation or proof of a charge assignment.

Timing is default off. Its pair gate uses time differences minus signed chord time of flight, including a configured mass and model allowance. An unknown common production time cancels. Missing/invalid time variances follow an explicit permissive policy. Curved path, clock calibration and detector-response validation remain necessary.

## Tests

Run `cmake -S src/tests/b0_standalone -B build-b0 -DEIGEN3_INCLUDE_DIR=/path/to/eigen3`, build, and run CTest. Tests include field/no-field fits, exact-zero curvature, correlated charge proposals, displaced offsets, missing/out-of-time measurements, order invariance, hard budgets, three-station acceptance, and front/back extension. These synthetic tests are not detector performance measurements.
