# B0: support local seed reference surfaces and separate prompt/displaced reconstruction modes

## Motivation

The B0 stub seeder currently estimates a local trajectory through the B0 spectrometer but ultimately converts it to parameters on an origin-centered perigee because `CKFTracking` reconstructs every seed on a hard-coded `PerigeeSurface(0,0,0)`.

In the default B0 configuration, `constrainToBeamline=true`, so the seed direction is additionally constructed from the origin to the estimated B0 magnet-entrance point. That is a reasonable prompt-proton mode, but it is not the right abstraction for tracks first observed in B0 after a displaced decay.

This matters directly for the physics use case discussed in https://github.com/eic/EICrecon/issues/2930: B0 should be able to reconstruct secondary charged particles even when central detectors see little or nothing.

Relevant code:

- `src/algorithms/tracking/B0TrackerStubSeeder.cc`
- `src/algorithms/tracking/B0TrackerStubSeederConfig.h`
- `src/algorithms/tracking/CKFTracking.cc`
- `src/global/tracking/tracking.cc`

## Important distinction

An origin-centered perigee surface does **not** itself imply that a track originates at the interaction point: displaced tracks can have non-zero impact parameters.

The prompt assumption enters through the **construction of the seed mean** (`constrainToBeamline`) and through the current inability of the generic CKF wrapper to initialize from a natural local B0 surface.

## Goal

Allow B0 reconstruction to initialize from a state near the first useful B0 measurements and make the beamline/beam-spot constraint an explicit optional physics constraint rather than a structural requirement of the tracking chain.

Target conceptual state near a local reference plane:

```text
(x, y, tx = dx/dz, ty = dy/dz, q/p, t)
```

converted to ACTS bound parameters on a real B0 sensor/reference plane.

## Proposed work

### 1. Make seed reference surfaces explicit

Extend the seed-to-ACTS bridge so a seed can specify the surface on which its bound parameters are defined.

Possible implementation directions:

- add/reference a geometry/surface identifier in the EICrecon seed representation where feasible;
- provide a B0-specific conversion layer that constructs `BoundTrackParameters` on the intended ACTS plane;
- keep origin-perigee support for existing central seeds.

Do **not** store local-plane coordinates in fields that `CKFTracking` will interpret as perigee coordinates.

### 2. Split B0 seeding into explicit modes

Suggested modes:

- **Prompt mode**: may use the beamline/beam-spot constraint intentionally.
- **Displaced/local mode**: infer direction and momentum from B0 measurements only, with no origin constraint.

The distinction should be configuration-visible and reflected in benchmark provenance.

### 3. Avoid unnecessary 6 m extrapolation in the initialization

Initialize close to the first relevant B0 surface, then let ACTS propagate through the actual field/material. Extrapolate to an origin perigee only as an output/reference transformation when an analysis requires it.

### 4. Keep beam-spot information separate from the numerical seed prior

A beam-spot constraint is external information and may be useful for prompt tracks. It should not be silently folded into a mode intended for displaced tracks.

## Acceptance criteria

- [ ] `CKFTracking` (or a small reusable helper around it) can initialize from a non-origin ACTS reference surface.
- [ ] Existing central tracking behavior remains unchanged.
- [ ] B0 has explicit prompt and displaced/local configuration modes.
- [ ] Prompt mode reproduces the current prompt-proton performance within statistical uncertainty.
- [ ] Displaced/local mode reconstructs charged tracks whose production vertex is far outside the default prompt vertex window.
- [ ] Output tracks can still be extrapolated to the origin perigee when requested.
- [ ] Seed covariance is defined on, and transformed consistently from, the actual seed surface.
- [ ] Unit tests cover coordinate/covariance transformations between local B0 surfaces and output perigee parameters.

## Validation samples

- prompt protons across the current 8–41 GeV / 4–22 mrad benchmark region;
- positive and negative pions;
- charged daughters from displaced Lambda-like decays with vertices upstream of and near the B0 magnet;
- three-station and four-station hit patterns;
- tracks crossing only one front/back face at a station because of acceptance/module edges.

## Non-goals

- Changing central seeding to local surfaces in the same PR.
- Replacing ACTS propagation/fitting.
- Assuming that simply setting `constrainToBeamline=false` validates displaced reconstruction; the covariance, reference surface, and benchmark must be treated consistently.
