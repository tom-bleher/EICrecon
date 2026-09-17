# B0 benchmark: define a B0-crossing truth reference and add displaced/mixed-species validation

## Motivation

The current B0 benchmark is deliberately narrow: prompt generator-stable protons in the solenoid-off `epic_ip6_extended` configuration. That is appropriate for the first baseline, but it is not sufficient to judge whether B0 tracking solves the broader use cases discussed in eic/EICrecon#2930.

There is also a subtle validation trap in the current truth-seeded chain: `TrackParamTruthInit` applies default production-vertex cuts (`|x|,|y| <= 80 mm`, `|z| <= 200 mm`) before B0 truth seeds are formed. A displaced B0-visible charged particle can therefore disappear before the B0 CKF is tested.

Relevant code:

- `src/algorithms/tracking/TrackParamTruthInit.cc`
- `src/algorithms/tracking/TrackParamTruthInitConfig.h`
- `src/global/tracking/tracking.cc`
- `src/benchmarks/reconstruction/b0_tracking/analyze_b0_tracking.py`
- `src/benchmarks/reconstruction/b0_tracking/README.md`

Related upstream discussion: https://github.com/eic/EICrecon/issues/2930

## Goal

Define a B0-specific truth reference whose denominator is based on **actual B0 reconstructibility**, not on prompt-track assumptions inherited from central tracking, and extend the benchmark to the particle species / vertex topologies that B0 is expected to reconstruct.

## Proposed truth definition

A useful B0 truth reference should start from MC particles and ask whether the particle leaves sufficient usable B0 information, e.g.:

- charged particle;
- crosses / produces associated hits in at least the configured number of **physical B0 stations**;
- optional minimum momentum / acceptance cuts stated explicitly in provenance;
- no implicit `|z_vertex| < 200 mm` prompt requirement unless the sample is intentionally prompt-only.

For each reconstructible MC particle, derive a truth state on the same reference surface used for comparison. For local/displaced modes this should be a propagated truth state at a B0 reference surface rather than a straight-ray approximation at the IP.

## Required benchmark categories

### Baseline prompt sample

- prompt protons;
- preserve the existing 8–41 GeV, ion-frame 4–22 mrad sample for continuity;
- split results by 3 vs 4 physical stations and front/back coverage.

### Charge/species checks

- pi+ and pi-;
- proton and anti-proton if relevant;
- optionally kaons/electrons where B0 acceptance makes sense;
- explicitly report wrong-charge reconstruction rate.

### Displaced samples

Include charged daughters produced:

- upstream of B0;
- near/inside the B0 magnet acceptance region where appropriate;
- across a broad `z_vertex` range motivated by Lambda/spectator studies from issue #2930.

### Occupancy/background samples

Add at least one sample with realistic or stress-test hit multiplicity. The benchmark should distinguish pattern-recognition failures from fit failures.

## Metrics

For every bin/sample report:

- `epsilon_seed`;
- `epsilon_track|seed`;
- `epsilon_total`;
- fake rate;
- duplicate rate;
- wrong-charge rate;
- momentum, phi, theta, and impact/reference-state residuals;
- pull widths and multivariate `D^2` once the independent final refit exists;
- reconstruction/runtime tails;
- fraction of events where combination/fallback budgets are hit.

Useful binning:

- momentum;
- ion-frame angle;
- production `z`;
- physical station multiplicity;
- front/back face multiplicity;
- total B0 hit occupancy;
- distance to module/aperture edges when available.

## Matching definitions

Do not define every reconstructed track without a match to the *target sample* as a fake. Once mixed species/backgrounds are included, distinguish:

1. a genuinely fake/combinatorial track with no valid MC association;
2. a real reconstructed particle that is outside the selected signal denominator;
3. a duplicate reconstruction of an already matched truth particle.

Persist these definitions in `provenance.json`.

## Acceptance criteria

- [ ] Add a B0-specific reconstructibility definition independent of `TrackParamTruthInit` prompt vertex cuts.
- [ ] Preserve the existing prompt-proton benchmark as a named baseline.
- [ ] Add both charge signs and at least one non-proton charged species.
- [ ] Add displaced-track sample(s) with production vertices well beyond `|z|=200 mm`.
- [ ] Add occupancy/background validation or a reproducible synthetic stress test.
- [ ] Compare truth-seeded and data-seeded chains using the same B0 reconstructible denominator.
- [ ] Correct fake/duplicate classification for mixed samples.
- [ ] Record all geometry, material-map, field, code SHA, selection and matching definitions in provenance.

## Non-goals

- Globally relaxing central-tracking truth cuts as a shortcut.
- Claiming a displaced-track result from a truth denominator that already removes displaced tracks.
