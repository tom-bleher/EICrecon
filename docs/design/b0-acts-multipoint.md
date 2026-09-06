# ACTS multipoint seeding and CKF for B0

This branch uses the native ACTS multipoint parameter estimator followed by
ACTS CKF. It requires **ACTS 47.7.0 or newer**; 47.7.0 is the validated version.
The dependency is enforced by the normal CMake version requirement. There is
one estimator implementation, without an extracted compatibility copy or a
runtime switch to the previous seeder.

## Algorithm and output

The B0 stub candidate finder supplies three or more station measurements.
`B0TrackerStubSeeds` and `B0TrackerMeasurements` feed `B0TrackerSeeds`, which
re-estimates every candidate with ACTS using all its measurements, ordered along
the outgoing ion direction. The estimator uses uniform weights without
geometric refinement and samples the magnetic field at the candidate midpoint.
Its homogeneous-field approximation is a seed approximation; CKF propagates
through the detector field and material.

`B0TrackerSeedParameters.surface` identifies the first sensor. Its local
coordinates, angles and covariance refer to that sensor, while the final fitted
track target remains the origin. CKF resolves the sensor ID and preserves the
original seed relation when rejecting invalid input. Consumers of seed
parameters must honor their surface ID. The intermediate stub collections are
optional diagnostic outputs.

The seed covariance propagates both local coordinates of every hit, including
local xy covariance, through a numerical Jacobian on fixed sensor surfaces.
Phi differences are wrapped. Angular and relative-q/p model uncertainties are
added before covariance inflation. CKF reuses those hits: inflation weakens this
finite, hit-derived prior but does not make it independent or calibrate the
final covariance.

## Environment and running

On the validation machine, enter the installed native ACTS environment from
this checkout, then build and source this checkout's installation:

```bash
~/eic/eic-shell-acts47
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --target install --parallel
source install/bin/eicrecon-this.sh
source ~/eic/dev/epic/install/bin/thisepic.sh epic_ip6_extended_5x41
```

The launcher provides ACTS 47.7.0 headers and libraries in the pinned
`eic_xl-26.04.0-stable.sif` environment. On another machine, use an eic-shell
installation with matching ACTS headers, CMake package and runtime libraries.
For a separate ACTS prefix, source its `bin/this_acts.sh` and configure with
`-DActs_DIR=/path/to/acts/install/lib/cmake/Acts`; EICrecon and its loaded plugins
must all be built against that installation. An older container's headers or
plugins must not override it. Use a fresh build directory when changing ACTS.

The tested geometry is **`epic_ip6_extended_5x41`**, electron 5 GeV × proton
41 GeV, crossing angle −0.025 rad. The unsuffixed `epic_ip6_extended` selects a
different beam-field configuration in this local geometry installation.
Run from a directory where its relative field-map and calibration paths resolve:

```bash
eicrecon -Pjana:nevents=100 -Pnthreads=1 \
  -Pacts:MaterialMap="$DETECTOR_PATH/calibrations/materials-map.cbor" \
  -Ppodio:output_file=multipoint.edm4eic.root input.edm4hep.root
```

Use the original simulation input: a reconstruction file containing
`B0TrackerSeeds` can supply its existing collection instead of running the new
factory. When trimming output collections, retain the complete PODIO relation
closure. The truth checker also needs `EventHeader`, `MCParticles`,
`B0TrackerHits`, B0 raw-hit associations and the B0 seed/measurement/track outputs.

The retained defaults are:

| Parameter | Default | Meaning |
|---|---:|---|
| `tracking:B0TrackerSeeds:covarianceInflation` | 100 | Weaken the reused hit prior |
| `tracking:B0TrackerSeeds:scatteringScale` | 0.0003 | Angular model coefficient, GeV rad |
| `tracking:B0TrackerSeeds:qOverPRelativeUncertainty` | 0.025 | Relative q/p model uncertainty |
| `tracking:B0TrackerCKFTrajectories:Chi2CutOff` | 15 | Individual hit compatibility χ² cutoff |
| `tracking:B0TrackerCKFTrajectories:NumMeasurementsCutOff` | 10 | Maximum measurement candidates per surface |
| `tracking:B0TrackerCKFTrajectories:NumMeasurementsMin` | 3 | Minimum measurements per track |

There is no additional final-track χ²/ndf or beam-compatibility cut. B0 ambiguity
resolution also accepts three measurements. Candidate-finder settings use the
prefix `tracking:B0TrackerStubSeeds`.

## Physics validation and limits

The same multipoint + CKF physics configuration was tested with native ACTS
47.7.0 (`2790f1b05c0c94cb3b569cbb18262b24867c7855`) and EICrecon
`0bff1394558a8d41fdd3c4259534b82264d30013`. A fresh sample, disjoint from the
initial 5,000 DVCS and 1,000 gun events, gave:

| Fresh sample | Matched / accepted protons | Efficiency | q/p core width | >20% tails | Joint 95% coverage |
|---|---:|---:|---:|---:|---:|
| DVCS, 2,000 events | 872 / 964 | 90.46% | 4.06% | 19 | 73.05% |
| 38 GeV proton gun, 1,000 events, 4–22 mrad | 434 / 474 | 91.56% | 3.82% | 10 | 71.43% |

Acceptance requires a primary forward proton with at least three clean B0
stations; matching requires association purity ≥0.5. Multiple matched tracks
are selected by measurement count, then χ², then index, never nearest truth.
Core width is half the 16th–84th percentile span of
`100*((q/p)_reco/(q/p)_truth - 1)`, with NumPy linear interpolation and no tail
or wrong-charge clipping. Percentages are rounded to two decimals. Covariance
coverage uses the persisted five-parameter origin covariance, excluding time.

Multipoint improves matching efficiency over the previous beam-informed stub,
but has worse momentum resolution and tails. The final covariance remains
uncalibrated despite finite positive-definite matrices. Restricting CKF
branching and simple quality cuts did not provide a reliable improvement;
constrained refit experiments failed covariance or convergence validation.
The branch retains the tested unconstrained reconstruction. These samples do
not validate denser backgrounds or other beam configurations.

Local reproducibility records are under
`/nas/nas1/tomble/b0_validation/20260906_multipoint_practical/`: its `README.md`,
`ckf_trial/provenance.json`, `quality_audit/heldout_summary.json` and
`covariance_audit/summary.json` record commands, software, inputs and results.
The installed geometry comes from epic
`2200b0aa76ed398e3e3d752e9d5346b43e171405` with local modifications; the entry
XML SHA256 is `611f6fe1fbef6f52d4c1ab05ab2268d1cb7327cde08e40552c3c7385194eb94f`
and material-map SHA256 is
`62cbf103e6c41b77397ddba591c53664359172d486da886702c6eca0da2c77b7`.

For local algorithm and persisted-output checks:

```bash
build/src/tests/algorithms_test/algorithms_test \
  '[B0TrackerStubSeeder],[b0acts]'
python3 src/scripts/check_b0_multipoint.py multipoint.edm4eic.root
```

The native-only cleanup was built and installed with ACTS 47.7.0. Both CTest
suites pass, including 7,852 assertions in 72 algorithm test cases. On 500 DVCS
and 500 gun events, all 224 persisted B0/EventHeader leaf branches match the
previous native build exactly after aligning event identities, including seed
and fitted covariances and hit relations. CMake rejects ACTS 44.4.0 with the
normal minimum-version error. Local validation records are in
`/nas/nas1/tomble/b0_validation/20260906_native_branch/`.
