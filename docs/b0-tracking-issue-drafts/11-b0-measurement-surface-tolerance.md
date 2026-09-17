# B0: replace the ad-hoc measurement surface tolerance with a tested precision contract

## Motivation

`TrackerMeasurementFromHits` currently special-cases B0 by increasing the ACTS `globalToLocal` on-surface tolerance from 0.1 um to 1 um. The source comment calls this an "ugly hack".

A source-level review shows that the larger tolerance may be numerically justified rather than arbitrary: reconstructed tracker-hit global positions are stored as single-precision `Vector3f`. Around the B0 z position (~6–7 m, represented in mm), float32 spacing is roughly 0.49 um, so coordinate quantization alone can exceed a 0.1 um surface tolerance even for a correctly placed hit.

The current implementation has two problems:

1. the tolerance is selected through a hard-coded B0 detector-ID branch inside a shared algorithm;
2. there is no detector-backed round-trip test establishing the minimum tolerance required by the actual serialization and DD4hep/ACTS surface transforms.

The goal is **not** to force B0 back to the ACTS default tolerance. The goal is to replace an undocumented detector-name exception with a measurable, configurable precision contract.

## Required study

For every representative B0 sensitive surface, including module/sensor edges and both front/back orientations, test:

```text
DD4hep cell center / reconstructed hit position
 -> EDM4eic Vector3f storage
 -> ACTS globalToLocal
 -> local coordinates
 -> ACTS localToGlobal
```

Record at least:

- normal-distance residual to the ACTS surface before `globalToLocal`;
- local-coordinate round-trip residual;
- dependence on global z / station;
- dependence on sensor orientation and edge position;
- failures as a function of tolerance.

Separate expected float serialization error from genuine geometry/surface misalignment.

## Implementation direction

Introduce an explicit configuration/capability for measurement conversion rather than branching on `B0Tracker_Station_1_ID` inside the generic algorithm. Possible approaches include:

- a configurable `onSurfaceTolerance` parameter on `TrackerMeasurementFromHits`, with detector factories choosing appropriate values; or
- preserving the generic default and creating a B0-specific factory/configuration wrapper that supplies the validated tolerance.

Whichever approach is used, keep central-detector defaults unchanged unless separately validated.

The tolerance should be interpreted as a numerical conversion tolerance, **not** as sensor thickness or detector spatial resolution.

## Failure accounting

A failed surface lookup or failed `globalToLocal` currently causes the measurement builder to warn and skip the hit. For B0 validation runs, add explicit counters and an optional strict mode so silently lost measurements cannot masquerade as tracking inefficiency.

Recommended counters:

- input tracker hits;
- missing volume/surface mapping;
- `globalToLocal` failures;
- successfully created measurements;
- failures by B0 station/layer where possible.

## Acceptance criteria

- [ ] The required B0 tolerance is established with an end-to-end float-storage/surface round-trip test.
- [ ] B0-specific tolerance is supplied through explicit configuration/capability, not a hidden detector-ID branch in the shared algorithm.
- [ ] The configured value is documented as a numerical tolerance, not a detector resolution.
- [ ] Central detector behavior is unchanged by default.
- [ ] Missing surface and coordinate-conversion failures are counted; strict validation mode can fail a run on unexpected losses.
- [ ] Tests cover all four stations, both faces, representative sensor centers and near-edge cells.
- [ ] Geometry mismatch remains distinguishable from expected float precision.

## Companion geometry work

`tom-bleher/epic` PR #2 provides the DD4hep/ACTS surface audit and source inventory. The geometry-side stable-identifier issue in ePIC PR #1 should ensure the same physical sensor can be identified throughout this round trip.

## Non-goals

- Reducing the tolerance merely to remove a special case.
- Using tolerance inflation to hide a displaced or incorrectly oriented ACTS surface.
- Treating a successful coordinate conversion as proof of calibrated hit covariance.
