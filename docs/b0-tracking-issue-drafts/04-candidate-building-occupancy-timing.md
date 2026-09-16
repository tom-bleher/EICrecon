# B0 seeding: harden candidate building for timing, two-projection compatibility, and high occupancy

## Motivation

`B0TrackerStubSeeder` already does significantly more than brute-force station combinations: it groups hits by physical station, ranks outer-station endpoint pairs, applies a non-bending-plane road, splits the combination budget across station subsets, and retains 3-station fallbacks.

The remaining candidate-building logic is still tuned for low occupancy and prompt proton studies:

- endpoint generation initially scales as `N_first * N_last` and all compatible endpoint pairs are stored before the fit-combination budget is applied;
- early pruning is dominated by the non-bending projection;
- front/back near-duplicates are identified by a transverse proximity heuristic (`sharedHitDistance`);
- B0 hit reconstruction already carries ~30 ps timing, but the seeder currently does not use hit-time compatibility.

Relevant code:

- `src/algorithms/tracking/B0TrackerStubSeeder.cc`
- `src/algorithms/tracking/B0TrackerStubSeederConfig.h`
- `src/detectors/B0TRK/B0TRK.cc`
- `src/algorithms/tracking/B0ReconstructionCounters.*`

## Goal

Make B0 candidate building robust under realistic occupancies without sacrificing efficiency for missing/single-sided measurements, and reduce reliance on purely geometric proximity for front/back duplicate handling.

## Proposed work

### 1. Add explicit occupancy instrumentation before changing algorithms

Record per event:

- hits per physical station;
- front/back hit counts where resolvable;
- number of outer endpoint pairs considered;
- number surviving each compatibility stage;
- number of middle-station candidates per endpoint pair;
- number of fit combinations attempted;
- whether any budget/truncation was hit;
- candidate-family multiplicities before/after deduplication.

Expose these in the benchmark JSON/report.

### 2. Bound endpoint-pair work

`maxCombinations` currently bounds fitted combinations, not the initial endpoint-pair cross product. Add an explicit strategy so pathological events cannot allocate/process arbitrarily many endpoint pairs.

Possible approaches to benchmark:

- spatial binning / KD-tree lookup in one or both transverse coordinates;
- sorted-window scans;
- top-K compatible endpoints per first/last hit;
- configurable global endpoint budget with deterministic ranking.

Any cap must be observable through counters and tested for efficiency bias.

### 3. Use both projections earlier

Add a cheap bending-plane compatibility estimate before the full polynomial/field fit. It does not need to be a precise momentum fit; it should reject combinations that imply impossible curvature/slope given broad B0 momentum and field bounds.

This should complement, not replace, the current conservative non-bend curvature allowance.

### 4. Add time-of-flight compatibility as an optional candidate gate

B0 digitization currently configures 30 ps timing. Use hit-time **differences** so an unknown common production time cancels:

```text
r_t(i,j) = (t_j - t_i) - (s_j - s_i)/(beta*c)
```

Start with a loose, configurable gate that includes measurement resolution, path-length uncertainty and the mass/momentum hypothesis uncertainty.

Timing should initially be optional and benchmarked independently; it may be most useful for out-of-time/background rejection rather than separating two prompt tracks from the same collision.

### 5. Replace proximity-only front/back equivalence with trajectory compatibility

`sharedHitDistance` is useful for suppressing front/back clones, but two genuinely close trajectories can also be within the same distance threshold.

Prefer a compatibility decision based on the current candidate prediction and uncertainty at the partner face. Preserve exact measurement identity in all cases.

Requirements:

- a missing front or back hit must not invalidate a station;
- aperture/module-edge inefficiency must remain supported;
- two close physical tracks must be allowed to form separate seed families.

## Acceptance criteria

- [ ] Add occupancy and candidate-stage counters to the benchmark.
- [ ] Bound outer endpoint work independently of `maxCombinations`.
- [ ] Add a cheap bend-plane preselection and quantify its efficiency/fake/runtime impact.
- [ ] Add optional hit-time compatibility with a default that is either disabled or demonstrably safe.
- [ ] Replace or augment `sharedHitDistance` with uncertainty-aware trajectory compatibility.
- [ ] Add close-track regression tests showing that two nearby tracks are not merged solely because their front/back measurements are close.
- [ ] Add tests for single-sided station coverage and module/aperture edges.
- [ ] Show runtime scaling versus hit occupancy and record truncation rates.

## Validation

Compare candidate-building variants on:

1. single prompt proton events;
2. two nearby tracks with controlled separation;
3. signal + random in-time hits;
4. signal + out-of-time hits;
5. realistic background overlays when available.

Measure efficiency, fake/duplicate rates, candidate multiplicity, endpoint-pair counts, fit-combination counts and wall time.

## Non-goals

- Replacing the fitter.
- Introducing ML before deterministic candidate-building limits are understood.
- Requiring a front/back pair at every physical station.
