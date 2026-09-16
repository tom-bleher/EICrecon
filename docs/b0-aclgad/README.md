# Opt-in B0 physical-pitch AC-LGAD response

This implements the detector-response follow-up to `tom-bleher/epic#1` and the
local tracking stack in EICrecon #4–#7. It does **not** change central tracking,
replace the legacy B0 default chain, or claim a calibrated B0 sensor response.

## What is implemented

`b0_aclgad` is a separately loaded JANA plugin. Its `B0ACLGADResponse` factory
converts `B0TrackerHits` into physical-channel ADC/TDC values, channel truth
relations, diagnostic channel hits, and clustered `Measurement2D` values on the
actual sensitive ACTS planes. Loading `b0_telescope` with
`-Pb0_telescope:response=aclgad` feeds those clusters to all four existing local
truth/realistic and direct/CKF comparison chains. The default remains `effective`.

The geometry must have `B0ACLGADReadoutVersion=1`, a matching `CartesianGridXY`
pitch/offset and complete cells inside each rectangular sensor. The companion
`epic/scripts/b0_readout/make_aclgad_profile.py` creates a fresh physical-pitch
profile without modifying the baseline geometry. Re-simulate with that profile:
a previously produced effective-70-um file must not be relabelled as physical
500-um data. The adapter checks geometry and channel round trips; the benchmark
also requires a completed, hash-matched simulation receipt.

## Response model and units

The dependency-light kernel `B0ACLGADResponse.h` uses mm, ns, eV and equivalent
output electrons. Each simulated energy deposit gives `energy / pairEnergy *
gain` nominal signal. Gaussian pad-area integrals distribute it across cells.
This is a **phenomenological spatial response**, not a microscopic drift model,
a claim that AC charge sharing is Gaussian, or a waveform simulation. The model
has no Fano fluctuations, per-deposit avalanche fluctuations, bipolar pulse
shape, dead-channel map, or multi-pulse channel reconstruction.

Deposits inside one configured time gate are summed on each physical pad before
fixed per-channel gain variation, additive electronic noise, ADC quantization,
and thresholding. There is no per-deposit threshold that would discard two
subthreshold contributions whose sum is detectable. The model does not
renormalize charge lost beyond sensor edges. Noise-only channels are supported
on every supplied sensor. Sensor/clock jitter and independent channel jitter
are separate; reconstruction does not divide common timing uncertainty by the
number of pads. A configurable amplitude time-walk term is corrected using the
measured amplitude. The TDC step must be an integer number of picoseconds for
the current RawTrackerHit representation.

Four-neighbor connected pads within one sensor form clusters with a bounded
**whole-cluster** time span. One channel can contribute to only one cluster.
Distinct particles in the same pad/time gate are unresolved pileup, not split
using Monte Carlo identity. Multiple clusters may exist on a sensor. Saturated
channels remain in the raw output; saturated clusters are rejected by default.
Disabling that rejection is a diagnostic mode, not a validated saturation model.
Input, channel and work budgets fail the event explicitly rather than quietly
returning a physically biased truncated event.

Cluster positions are amplitude-weighted centroids. Covariance includes the
centroid/noise/ADC Jacobian and explicit calibration terms for model error,
single-pad response and edges, with a spatial correlation parameter. Position-
time cross terms are retained. Both spatial and full three-dimensional
covariance positivity are checked. This is **not** proof of calibrated errors:
threshold selection, gain variation and nonlinear edge bias require measured
calibration and pull/coverage studies. The spatial CKF still fits two coordinates;
cluster time is available to optional candidate timing, not silently promoted
to a six-dimensional timing measurement in the CKF.

## Calibration is required, not fabricated

`calibrationFile` is required. Every field is explicit in a strict schema-versioned
key/value file; missing, duplicate, unknown, non-finite and overflowing values
fail. `experimental.cfg` is an **uncalibrated example** and is rejected unless
`allowExperimental=true` is supplied. Its 500-um pitch, gain/noise parameters,
50-um model-error term and timing components are examples, **not** a measurement
of the B0 design. A `status=calibrated` label is supplied by the calibration
producer; it does not independently validate the model. Record a real source
in `reference` and retain the calibration file hash.

The spatial response width and covariance floors need sensor/electronics data
covering subpixel positions, thresholds, incidence angles, fluence/operating
conditions, gains, timing and sensor edges. Replacing the model should preserve
this input/output contract. Do not tune tracking cuts to conceal miscalibrated
hit uncertainty. Relevant primary context: https://arxiv.org/abs/2211.13809 ;
these studies motivate charge-sharing models but do not validate these example
parameters or this Gaussian surrogate.

## Collections and truth accounting

| Output | Semantics |
|---|---|
| `B0ACLGADRawHits` | Physical channel ID; `charge` is ADC counts, timestamp is ps. Do not use the legacy raw-energy conversion. |
| `B0ACLGADRawHitLinks/Associations` | Every contributing SimTrackerHit and its signal fraction; positive noise remains unassociated. |
| `B0ACLGADChannelHits` | Pad-centre diagnostics linked one-to-one to raw channels; pitch-based errors are **not** cluster tracking errors. |
| `B0ACLGADMeasurements` | Cluster position/time/full Cov3f on a sensitive ACTS surface, all channel-hit relations, normalized amplitude weights. |

`B0MeasurementTruth.h` weights channels inside each measurement, then the track
exporter weights physical measurements equally. Missing/noise fractions remain
in the denominator. Legacy unnormalized raw associations are normalized only
downward; a small signal tail is never promoted to unit purity. The local truth
seeder uses the same rules. The AC-LGAD benchmark adapter also matches seeds at
the measurement level, so a large cluster does not count as several independent
tracking measurements or dominate the truth vote merely because it has more pads.
The legacy B0 seeder/CKF and central tracker are unchanged.

## Running (requires a built, matched eic-shell stack)

First build the companion physical profile and simulate through its `simulate`
subcommand. Then, from an EICrecon checkout containing this stack:

```sh
python3 src/benchmarks/reconstruction/b0_telescope/run_aclgad_comparison.py \
  --profile /path/to/physical-profile \
  --simulation-receipt /path/to/simulation-run/simulation.json \
  --geometry-calibration /path/to/calibration.xml \
  --material /path/to/matching-material-map.cbor \
  --response docs/b0-aclgad/experimental.cfg --allow-experimental \
  --out /path/to/new-response-run
```

Use `--asset` for dynamic libraries and plugin-opened calibration/field assets,
`--source` to record source state, `--set FactoryTag:parameter=value` for supported
tracking settings, `--prior-scale` for fresh-prior scans, and `--dry-run` to inspect
an invocation without executing reconstruction. No shell command interpolation
is used. Profile, simulation, calibration and output hashes are checked. Source
identity is not proof of installed-library identity. A completed reconstruction
and a completed analysis are recorded separately.

Outputs include the existing same-surface track/seed report and
`cluster_residuals.jsonl`, with per-cluster multiplicity, majority matching,
missing-reference status, spatial/time residuals, pulls and D2 when truth exists.
Compare against the existing effective-response benchmark using paired physics
samples generated with the corresponding geometry; do not relabel one simulated
ROOT file for both profiles. No measured efficiency gain is asserted here.

## Tests and promotion gates

```sh
cmake -S src/tests/b0_aclgad -B build-response -DCMAKE_BUILD_TYPE=Release
cmake --build build-response --parallel 2
ctest --test-dir build-response --output-on-failure
python3 -m unittest discover -s src/benchmarks/reconstruction/b0_telescope -p test_aclgad.py -v
```

The standalone tests exercise threshold-after-summation, loss at edges, channel
and sensor identity, deterministic noise/gains, timing, saturation, budgets,
calibration rejection and cluster-aware truth weights. The Python tests use
synthetic records; mocked simulator outputs are not detector events.

Before promotion from draft: compile the plugin against supported ACTS/DD4hep/
EDM/JANA versions; run the actual geometry/profile audit; demonstrate channel-ID
round trips, rotated local covariance and PODIO relations; test both signs,
multiple species, missing layers and background; compare effective and physical
response on reconstructed samples; calibrate residuals/coverage versus subpixel
position, edge distance, angle, threshold and charge; scan the independent final
prior and verify material mapping. Actual detector response and covariance
coverage have not been measured by these dependency-light tests.
