# B0 seeding: make alternate-charge and unresolved-curvature hypotheses self-consistent

## Motivation

The current B0 stub seeder has useful machinery for ambiguous charge sign:

- infer charge from fitted curvature when significant;
- emit both charge hypotheses below a configurable curvature significance;
- optionally force or test both charges.

There are two edge cases that should be made explicit and regression-tested:

1. a candidate with exactly `q/p == 0` is rejected in the compatibility path before the charge-hypothesis helper can treat the sign as unresolved;
2. when both charge hypotheses are emitted, the code flips the sign of `q/p` and the covariance correlations involving `q/p`, while the rest of the seed state comes from the originally fitted signed trajectory.

For sufficiently ambiguous curvature this may be numerically harmless, but an alternate charge hypothesis should ideally be a **self-consistent state hypothesis**, not only a sign change in one component.

Relevant code:

- `src/algorithms/tracking/B0TrackerStubSeeder.cc`
- `src/algorithms/tracking/B0TrackerStubSeeder.h`
- `src/algorithms/tracking/B0TrackerStubSeederConfig.h`

## Goal

Define one consistent policy for unresolved curvature / charge ambiguity and ensure every emitted seed corresponds to a physically and numerically coherent starting state for ACTS propagation.

## Proposed work

### 1. Define the unresolved-curvature regime

Use a significance criterion based on the fitted `q/p` covariance rather than a special exact-zero branch alone.

Questions to settle explicitly:

- what minimum finite `|q/p|` should be emitted when the fit is compatible with zero curvature?
- should the two hypotheses use symmetric `+/- |q/p|` proposals based on the fit uncertainty?
- when should an effectively infinite-momentum hypothesis be rejected rather than regularized?

The policy should avoid literal infinite momentum and unstable propagation.

### 2. Recompute a consistent alternate-charge seed state

When emitting both charge signs, recompute any charge-dependent back-propagation / field correction needed for the alternate hypothesis rather than blindly copying the original signed-fit direction.

At minimum validate whether the current state difference is negligible over the configured ambiguous-curvature region. If it is negligible by construction, document and test that approximation quantitatively.

### 3. Preserve covariance consistency

For a sign transformation of `q/p`, transform all affected covariance cross terms consistently. If the alternate hypothesis is recomputed nonlinearly, propagate/recompute the covariance around that hypothesis rather than relying only on row/column sign flips.

## Acceptance criteria

- [ ] Exact-zero and near-zero curvature candidates follow a documented common policy.
- [ ] No seed contains non-finite momentum/parameters/covariance.
- [ ] Positive and negative truth samples reconstruct with comparable efficiency after accounting for detector/field asymmetries.
- [ ] Wrong-charge rate is reported in the B0 benchmark.
- [ ] Unit tests cover exact zero, sub-threshold significance, threshold boundary, and clearly resolved positive/negative curvature.
- [ ] A test verifies that both emitted charge hypotheses can propagate to the B0 measurement region without pathological failure.
- [ ] If alternate-charge direction/state reuse is retained as an approximation, document the numerical bound showing it is negligible in the ambiguity region.

## Validation

Use controlled particle-gun samples with both signs over the momentum range where curvature significance transitions from resolved to unresolved. Plot:

- fitted curvature significance;
- emitted charge-hypothesis multiplicity;
- correct/wrong charge after final fit;
- CKF success rate by charge hypothesis;
- parameter residuals near the ambiguity threshold.

## Non-goals

- PID determination.
- Forcing both charge hypotheses for all B0 tracks once curvature is well measured.
