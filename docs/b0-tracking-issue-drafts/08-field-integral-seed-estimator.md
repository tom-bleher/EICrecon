# B0 seeding: prototype a local field-integral seed estimator as an alternative to the parabola approximation

## Motivation

The current B0 stub estimator is physically motivated and already performs well for the prompt-proton sample: it rotates to the ion frame, fits a parabola in the bend plane, samples the ACTS/DD4hep field, corrects non-bend quadrupole curvature, and analytically steps the result back toward the magnet entrance.

This is a reasonable production baseline and should **not** be replaced merely for elegance.

However, the implementation has accumulated detector-specific corrections for:

- midpoint field sampling;
- the entrance-to-first-station field gap;
- the tilted hard field entrance;
- non-bend quadrupole correction;
- back-extrapolation to the seed reference.

A compact alternative worth benchmarking is a **linearized field-integral local-state fit** that uses the same field provider but expresses both transverse projections in one least-squares model.

Relevant code:

- `src/algorithms/tracking/B0TrackerStubSeeder.cc`
- `src/algorithms/tracking/B0TrackerStubSeeder.h`
- `src/algorithms/tracking/B0TrackerStubSeederConfig.h`

## Goal

Implement an experimental seed estimator that returns a local state near B0 from field integrals and compare it with the current parabola/two-sample estimator. Keep the current estimator as the default unless the new model demonstrably improves robustness, simplicity, or portability.

## Model sketch

For small transverse slopes and a local reference `z0`, define

```text
a0 = (x0, y0, tx0, ty0, eta)^T,   eta = q/p
```

with approximately

```text
x''(z) ~= -kappa * eta * By(z)
y''(z) ~= +kappa * eta * Bx(z)
```

For a fixed reference path, precompute

```text
Fx_i = -kappa * integral[z0->zi] (zi-u) By(r_ref(u)) du
Fy_i = +kappa * integral[z0->zi] (zi-u) Bx(r_ref(u)) du
```

and solve the weighted linear system

```text
x_i ~= x0 + tx0 * dz_i + eta * Fx_i
y_i ~= y0 + ty0 * dz_i + eta * Fy_i
```

for the five local parameters. One or two iterations can update the reference path/field integrals if necessary.

This is an **initialization model**, not a replacement for full ACTS transport/material fitting.

## Implementation requirements

- use the same ACTS/DD4hep magnetic-field provider as CKF/fitting;
- use measured hit covariance in the weighted fit;
- return a full local covariance including correlations;
- support 3- and 4-station candidates;
- support the realistic front/back geometry without assuming exactly one plane per station;
- expose diagnostics comparing estimated state against the current seed estimator and truth.

## Comparison to current estimator

Measure on identical candidates:

- q/p bias and resolution;
- direction bias/resolution;
- charge-sign correctness/significance;
- covariance pull widths;
- CKF/direct-fit convergence rate;
- sensitivity to field configuration;
- CPU time;
- implementation/configuration complexity.

Test specifically whether the current estimator's documented midpoint/two-sample approximation is already comfortably below measurement/material uncertainties. If so, retain it and close this issue as an evaluated alternative.

## Acceptance criteria

- [ ] Add the field-integral estimator behind an explicit experimental configuration switch or benchmark-only helper.
- [ ] Unit-test the linear model in uniform dipole and dipole+quadrupole toy fields.
- [ ] Compare against numerical ACTS propagation/truth in the real B0 field.
- [ ] Compare 3- and 4-station seed resolution with the current estimator.
- [ ] Compare downstream track-finding efficiency and runtime.
- [ ] Document whether the extra generality is worth the complexity.
- [ ] Do not change the default estimator without benchmark evidence.

## Non-goals

- Replacing full ACTS propagation with the linearized equations.
- Treating the small-angle model as exact for arbitrary fields/slopes.
- Blocking higher-priority work on the independent final refit, local seed surfaces, or displaced-track validation.
