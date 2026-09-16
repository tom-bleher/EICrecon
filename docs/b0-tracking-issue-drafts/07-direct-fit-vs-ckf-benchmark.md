# B0: benchmark direct candidate-to-Kalman fitting against the CKF chain

## Motivation

The current B0 chain is:

```text
B0 measurements
  -> dedicated station-aware stub seeder
  -> ACTS combinatorial Kalman filter
  -> ambiguity resolution
```

The dedicated seeder already performs substantial hit association: it selects one hit per physical station, applies compatibility cuts, ranks complete candidates, forms primary/fallback seed families, and stores the selected hits in each `TrackSeed`.

This raises a useful architecture question: once a B0 candidate already has a plausible measurement set, how much value is the **combinatorial** part of CKF still providing relative to fitting that candidate directly with ACTS?

A direct fit could be simpler, cheaper, and statistically clearer for low-ambiguity B0 events, while CKF could remain available when candidate extension/branching is genuinely required.

Related upstream discussion: https://github.com/eic/EICrecon/issues/2930

## Goal

Implement a controlled comparison between:

1. current `stub seed -> CKF` reconstruction;
2. `station-aware candidate -> direct ACTS Kalman fit` using the candidate-associated measurements;
3. optionally a hybrid mode that performs limited candidate extension before the direct fit.

Do not replace CKF until this comparison demonstrates an advantage on realistic samples.

## Proposed implementation

### A. Build an explicit measurement-set representation

For each accepted B0 seed family/candidate, preserve the exact associated B0 measurements/surfaces needed for fitting.

The direct-fit path should not need to rediscover the candidate from the full event measurement collection.

### B. Fit the candidate with ACTS

Use the established ACTS Kalman fitter with full field/material propagation. Initialization should follow the local-reference work in the corresponding B0 issue and final statistics should follow the independent-refit requirements.

### C. Optional limited extension

If realistic data show that candidate seeds often miss one front/back measurement or one station hit, test a bounded extension step:

- propagate candidate to expected B0 surfaces;
- collect compatible measurements within uncertainty-aware windows;
- branch only when there is genuine ambiguity;
- fit resulting measurement sets directly.

This is intentionally smaller in scope than a global CKF search.

## Comparison metrics

For the same input events and reconstructible denominator compare:

- `epsilon_total`;
- `epsilon_track|seed`;
- fake rate;
- duplicate rate;
- wrong-charge rate;
- 3-station recovery efficiency;
- final momentum/angle resolution;
- fit `chi2/ndf` and pulls after independent refit semantics are in place;
- number of propagated surfaces / branches;
- CPU time and high-percentile event latency;
- memory/candidate multiplicity.

## Required samples

- prompt proton baseline;
- positive/negative pion samples;
- displaced charged tracks;
- module-edge and single-sided station cases;
- two nearby tracks;
- realistic or stress-test occupancy/background overlays.

## Acceptance criteria

- [ ] Add an opt-in B0 direct-fit path without changing central tracking.
- [ ] Preserve the current CKF path as a reference.
- [ ] Run both paths from the same candidate-building output where possible.
- [ ] Benchmark accuracy, efficiency, fake/duplicate rates, runtime and memory.
- [ ] Document specific failure classes where CKF succeeds and direct fitting fails, and vice versa.
- [ ] If a hybrid candidate-extension mode is added, bound its branching explicitly and instrument it.
- [ ] Make no default switch until the benchmark shows the new path is at least as robust on the intended B0 physics samples.

## Interpretation

The purpose of this issue is not to prove that CKF is inappropriate for B0. Kalman fitting is well suited to a short telescope. The question is narrower: **does B0 still need global combinatorial track finding after the dedicated seeder has already associated measurements?**

## Non-goals

- Replacing ACTS fitting.
- Introducing a GNN/transformer before deterministic baselines are understood.
- Removing CKF support from EICrecon globally.
