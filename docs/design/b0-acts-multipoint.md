# Experimental ACTS multipoint seeding for B0

This branch starts at upstream main `cef494fd67c4c523e1c99a396c5f9d7c9e51108c`
(fetched 2026-09-06). It replaces B0's solenoidal seed estimation with the
previously tested B0 stub candidate finder followed by ACTS multipoint parameter
estimation. Central tracking and the B0 truth-seeded chain retain their settings.

## Algorithm and interfaces

The [ACTS multipoint estimator](https://github.com/acts-project/acts/blob/2790f1b05c0c94cb3b569cbb18262b24867c7855/Core/include/Acts/Seeding/EstimateTrackParamsFromSeed.hpp)
fits a transverse Taubin circle and a longitudinal line in the magnetic-field
frame. It accepts three or more ordered points and an arbitrary homogeneous
field direction. It is a parameter estimator: the B0-specific candidate finder
still supplies one hit per station, including three-station fallback candidates.

`B0TrackerStubSeeds` feeds `B0TrackerSeeds`, together with
`B0TrackerMeasurements`. The latter factory uses all candidate measurements,
ordered along the outgoing ion direction, with uniform fit weights and no
geometric refinement. It samples the actual field at the fitted path midpoint.
The field approximation and the existing candidate selection remain limitations
in a combined dipole/quadrupole magnet.

ACTS 44 lacks this API. A private MPL-2.0 extraction from ACTS 47.7.0 commit
`2790f1b05c0c94cb3b569cbb18262b24867c7855` supplies it, preserving upstream
copyright and source hashes. CMake selects the native implementation when its
compile/link probe succeeds. A whole-stack ACTS 47 build has not been validated;
native dispatch and the extracted arithmetic were checked separately.

**Changed seed convention:** `B0TrackerSeedParameters.surface` identifies the
first registered sensor. `loc`, angles and covariance refer to that sensor,
not the origin. CKF now resolves nonzero surface IDs; zero still means the
legacy origin perigee. An invalid surface is skipped without shifting subsequent
seed relations. The seed's `perigee` metadata and the final fitted track target
remain at the origin. Consumers must honor the parameter surface ID.

The full five-dimensional measurement covariance is propagated from both local
coordinates of every hit, including local `xy` covariance, through a central
numerical Jacobian on fixed sensor surfaces. Phi differences are wrapped. Local
coordinates are in mm and q/p in GeV^-1. Provisional angular and relative-q/p
model terms are added; time uses approximate relativistic flight time with an
independent 100 ns² variance. Charge reversal also reverses the q/p covariance
cross terms.

CKF reuses the seed's hits. Multiplying the seed covariance by 100 weakens this
correlated prior; it does **not** make the prior independent or establish
calibrated final uncertainties. B0 CKF and ambiguity resolution accept three
measurements. The upstream CKF chi-square cutoff of 15 and geometry padding of
1 mm are retained for both comparison modes.

## Running the comparison

Build and source this checkout inside eic-shell, then load the local
`epic_ip6_extended_5x41` geometry. The unsuffixed `epic_ip6_extended` selects a
different beam-field configuration in the tested geometry installation.

```bash
source /path/to/epic/install/bin/thisepic.sh epic_ip6_extended_5x41
source install/bin/eicrecon-this.sh
eicrecon -Pjana:nevents=100 -Pnthreads=1 \
  -Pacts:MaterialMap="$DETECTOR_PATH/calibrations/materials-map.cbor" \
  -Ppodio:output_file=multipoint.root input.edm4hep.root
```

Use the original simulation as input. Replaying a reconstruction file containing
`B0TrackerSeeds` can supply that existing collection instead of running the new
factory. Run from a directory where the geometry's relative field-map and
calibration paths resolve.

For the controlled stub baseline, add
`-Ptracking:B0TrackerSeeds:useMultipoint=false`. This copies the same candidate
finder's original origin-bound parameters into the usual output collections.
It is not upstream main's generic solenoidal B0 seeder. Parameters of the
candidate finder use the prefix `tracking:B0TrackerStubSeeds`.

Useful multipoint options, under `tracking:B0TrackerSeeds`, are:

| Parameter | Default | Purpose |
|---|---:|---|
| `covarianceInflation` | 100 | Weaken the reused measurement prior |
| `diagonalCovariance` | false | Diagnostic removal of covariance cross terms |
| `scatteringScale` | 0.0003 | Angular model sigma coefficient, GeV rad |
| `qOverPRelativeUncertainty` | 0.025 | Relative momentum model uncertainty |
| `minCurvatureSignificance` | 2 | Threshold for emitting both charge hypotheses |
| `charge` | 0 | Infer charge; ±1 forces a diagnostic hypothesis |
| `testBothCharges` | false | Emit both hypotheses for every candidate |

Persist `EventHeader`, the normal B0 seed/measurement/track collections and their
relation closure, plus `MCParticles`, `B0TrackerHits` and raw-hit associations,
for validation. Intermediate `B0TrackerStubSeeds` and
`B0TrackerStubSeedParameters` are optional diagnostic outputs, not default file
contents. `src/scripts/check_b0_multipoint.py` checks these outputs and reports
truth matching; `compare_b0_multipoint.py` compares its JSON reports on an
identical truth denominator. Run the relevant algorithm tests with:

```bash
build/src/tests/algorithms_test/algorithms_test \
  '[B0TrackerStubSeeder],[b0acts],[multipoint]'
```

## Validation provenance

Validated in `eic_xl-26.04.0-stable.sif`, ACTS 44.4.0, with the local epic source
at `2200b0aa76ed398e3e3d752e9d5346b43e171405` and its installed 5×41 geometry.
That epic checkout also contains local changes: the commit alone does not fully
identify its installation. Entry XML SHA256:
`611f6fe1fbef6f52d4c1ab05ab2268d1cb7327cde08e40552c3c7385194eb94f`;
material map SHA256:
`62cbf103e6c41b77397ddba591c53664359172d486da886702c6eca0da2c77b7`.

Inputs are the first 5,000 DVCS events and first 1,000 proton-gun events under
`/nas/nas1/tomble/presentation16/simulation_outputs/20260825_030157/combined/`:

- `dvcs_devepic/DVCS.1.ab.hiAcc.5x41/sim.edm4hep.root`
- `pg_devepic_p38/pgun.theta4to22/sim.edm4hep.root` (38 GeV gun, 4–22 mrad).

Input SHA256, in the same order:
`b8709613b58528b8a8da844c7843ef73b3b072be30c6ea77cb7cf71947a1d499`,
`44538ce81d78f377bbbd210471411953adf71c9cc71a1c648a4175205aec612c`.

Logs, input manifests, parameter scans, ROOT outputs and sensor-truth audits are
under `/nas/nas1/tomble/b0_validation/20260906_multipoint_integration/`.
Eight-worker runs were compared with one worker using event identities, rather
than file entry order.

The truth denominator is primary protons with at least three separated B0
stations, retaining more than half their original momentum and moving forward
at each selected station. Matching requires association purity ≥0.5. When
several final tracks match one proton, select by measurement count, then chi2,
then track index; no truth-nearest selection. Resolution is half the 16th–84th
percentile span of `100*(reconstructed q/p / truth q/p - 1)`, without tail or
wrong-charge clipping. Percentiles use NumPy linear interpolation; displayed
percentages are rounded to two decimal places.

The full install build and CTest suite pass. Dedicated tests cover both charges,
arbitrary field direction, tilted sensors, three/four points, invalid inputs and
6,000 independent local-noise throws. Their entire whitened five-dimensional
empirical covariance agrees with the predicted matrix within 0.10 per element.
The compatibility helper reproduces all 3,253 earlier trial candidates for five
estimator variants exactly. Both native dispatch and fallback pass the same
analytic tests.

Runtime checks exercise an invalid seed before a valid seed: reconstructed
parameters and covariance are bit-identical, with the surviving seed association
correctly retaining its original index. Sensor truth is projected from its
associated simulation hit onto the actual first sensor; origin truth is not
used to assess a sensor-bound state.

## Reconstruction results

The multipoint integration increases matched reconstruction efficiency but worsens final momentum
resolution, including on the same successfully reconstructed protons. It is an
experimental alternative, not a demonstrated improvement over the stub seeder.

| Sample | Accepted truth protons | Stub matched | Multipoint matched | Stub width | Multipoint width |
|---|---:|---:|---:|---:|---:|
| DVCS, 5,000 events | 2,375 | 1,976 (83.20%) | 2,151 (90.57%) | 3.12% | 3.98% |
| Proton gun, 1,000 events | 464 | 370 (79.74%) | 422 (90.95%) | 2.81% | 3.76% |

For the 1,964 DVCS protons found by both methods, the widths are 3.13% and
3.87%; the difference is +0.73 percentage points with paired bootstrap 95%
interval [0.55, 0.93]. For 368 common gun protons, they are 2.80% and 3.73%,
difference +0.93 [0.56, 1.48]. Bootstrap resampling uses 2,000 draws, seed
271828, with the matched primary proton as the resampling unit. These intervals are
conditional on the common-success sample, not systematic-error bounds. Multipoint gains
187 and loses 12 accepted DVCS protons; the gun counts are 54 and 2.

Among quality-selected accepted protons, multipoint has four wrong charges and
32 relative-q/p errors above 20% in magnitude in DVCS, and two/nine in the gun.
The stub baseline has none in these samples. No final tracks have association
purity below 0.5 in either mode. Multipoint produces three duplicates across
all truth particles in DVCS (one for an accepted primary proton); the gun has
none. These observations do not establish performance in denser backgrounds.

| DVCS covariance multiplier | Matched / 2,375 | Inclusive width | Common-track width |
|---|---:|---:|---:|
| 1 | 2,152 | 3.99% | 3.91% |
| 10 | 2,157 | 4.00% | 3.89% |
| 100 (default) | 2,151 | 3.98% | 3.87% |

Changing inflation does not remove the resolution regression. The default of
100 deliberately weakens the prior; it is not a tuned optimum. The stub
baseline uses an origin/beam-spot constraint, whereas the multipoint state is
estimated from telescope measurements alone. That difference and the local
homogeneous-field approximation are possible explanations, not isolated causes.
Testing a consistent vertex constraint in the downstream fit is a useful next
experiment. Adding the origin as another point in a single constant-field helix
would incorrectly extend the B0 field through the upstream region.

All 2,931 DVCS and 599 gun seed covariances are finite and positive definite
after float persistence. All 2,195 DVCS and 431 gun tracks contain a measurement
on the starting sensor. Two in each sample select a different hit on that same
sensor; these are traced CKF alternatives, not skipped starting surfaces.
The one/eight-worker DVCS comparison passes for all 5,000 events and 277 fields,
including recursive PODIO relations, object ordering and exact stored bits.

The first-sensor truth audit retains 2,883 unambiguous DVCS associations, reports
41 ambiguous cases and seven non-primary-proton cases separately, and evaluates
2,114 clean first-emitted candidates. Their joint five-parameter nominal 95%
coverage is 61.40% at inflation 1, 97.82% at 10 and 99.24% at 100. Thus the
fixed-field measurement Jacobian passes its controlled noise test, but the
measurement-plus-model covariance is not calibrated for full simulation. The
broad initialization must not be interpreted as a validated final-track error.
