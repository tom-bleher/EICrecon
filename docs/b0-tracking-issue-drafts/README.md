# B0 tracking follow-up issue roadmap

These are issue drafts for the `feat/b0-tracking` branch following a source-level review of the current B0 reconstruction and the realistic eight-layer / four-station geometry.

Upstream discussion: https://github.com/eic/EICrecon/issues/2930

## Recommended order

### Priority 1 — establish statistically trustworthy output

1. **[Independent final weak-prior refit](01-independent-final-refit.md)**
   - Keep hit-derived seeds for finding.
   - Refit selected measurements independently for final parameters/covariance.
   - Central tracking remains unchanged by default.

2. **[Local seed reference + prompt/displaced modes](02-local-reference-and-displaced-mode.md)**
   - Stop forcing the natural B0 seed state through an origin-perigee representation.
   - Make the beamline constraint optional and explicit.

3. **[B0-specific truth reference + validation matrix](03-b0-truth-reference-and-validation-matrix.md)**
   - Define reconstructibility from actual B0 station crossings/hits.
   - Add displaced tracks, both charge signs, mixed species, occupancy and correct fake classification.

### Priority 2 — harden pattern recognition

4. **[Candidate building under occupancy + timing](04-candidate-building-occupancy-timing.md)**
   - Instrument and bound endpoint combinatorics.
   - Use both transverse projections earlier.
   - Add optional time-of-flight compatibility.
   - Replace proximity-only front/back equivalence with trajectory-aware compatibility.

5. **[Charge ambiguity / zero-curvature edge cases](05-charge-ambiguity-zero-curvature.md)**
   - Make both-charge hypotheses physically self-consistent.
   - Define exact/near-zero curvature policy.

6. **[Shared physical-station mapping](06-shared-physical-station-map.md)**
   - One authoritative mapping for seeder, CKF station requirements and benchmark provenance.

### Priority 3 — compare simpler/cleaner alternatives

7. **[Direct candidate-to-fit reconstruction vs CKF](07-direct-fit-vs-ckf-benchmark.md)**
   - Test whether global combinatorial finding is still needed after the dedicated B0 candidate builder has already associated measurements.

8. **[Field-integral local seed estimator](08-field-integral-seed-estimator.md)**
   - Experimental alternative to the current parabola/two-field-sample estimator.
   - Keep lower priority unless the current approximation is shown to limit performance.

## Dependencies

```text
01 independent final refit
        |
        +-------> 03 benchmark/pulls become statistically meaningful

02 local seed surfaces/modes
        |
        +-------> 03 displaced benchmark
        +-------> 07 direct candidate fitting

06 shared station map
        |
        +-------> 03 consistent reconstructibility
        +-------> 04 robust candidate building

04 candidate hardening -------> 07 realistic comparison against CKF

08 field-integral estimator is largely independent and should not block 01-07.
```

## Geometry-side companion work

Two related drafts belong primarily in the `tom-bleher/epic` fork:

- validate the realistic front/back B0 surface/material mapping and reproducible ACTS material-map generation;
- replace the effective 70 um readout proxy with a realistic AC-LGAD response/digitization model, coordinated with EICrecon.

## Scope principle

B0-specific behavior should stay B0-specific where possible. Shared changes to `CKFTracking` should be reusable **capabilities** controlled by explicit default-off options, not implicit behavior keyed to `numB0StationsMin` or detector-name assumptions.
