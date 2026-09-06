# B0 seeding: implementation and comparison options

Research checked 6 September 2026 against `feat/b0-stub-seeder` at
`1cab33ff`, upstream EICrecon `cef494fd`, and the installed ACTS 44.4.0 in
eic-shell `26.04.0-stable`. The target is the realistic B0 geometry at epic
`2200b0aa`, using `epic_ip6_extended_5x41` and 5 × 41 GeV beams.

## Assessment

Recommendation (inferred from the source and comparisons below): keep the
dedicated B0 candidate finder and the existing ACTS CKF. The most
useful next estimator experiment is ACTS's multipoint estimator, following
the three-point comparison below. No turnkey B0 replacement was
verified. Distinguish a candidate finder, a seed parameter estimator, and a
track fitter: replacing one does not provide the others.

Confirmed from source: the [current algorithm](../../src/algorithms/tracking/B0TrackerStubSeeder.cc)
groups hits by physical station, prunes endpoint pairs and intermediate hits
in the non-bending projection, fits a parabola and line, and samples the field
at the candidate midpoint and entrance gap. It corrects the non-bending
quadrupole contribution, ranks and deduplicates candidates, and constructs
origin-perigee parameters with propagated measurement covariance, beam-spot
covariance, and calibrated scattering/model terms. The CKF then uses the
full field and material model. The August comparative DOCX describes an older
sampled-field-integral implementation; its implementation-specific advice
should not be applied to this revision without rechecking the source.

This division of work has professional precedent. LHCb Hybrid Seeding uses
separate bending and non-bending projections, a detector-specific polynomial
model, search windows, and subsequent track extension. Its station layout and
stereo measurements differ from B0, so its implementation is not a direct
substitute. [Hybrid Seeding paper, sections 2–3](https://arxiv.org/html/2007.02591v2).
FASER likewise constructs candidates before fitting with its field map and
material model. [FASER tracker description](https://faser.web.cern.ch/tracker).

Upstream EICrecon still registers the generic `TrackSeeding` for B0 at the
checked revision. Its defaults and origin assumptions are not a validated
baseline for the B0 telescope; the reported production failure is documented
in [issue 2746](https://github.com/eic/EICrecon/issues/2746). That issue also
concerns displaced secondaries: successful prompt-proton reconstruction does
not establish a solution for those decays.

## Existing alternatives, in trial order

| Alternative | Availability | Useful experiment and limitation |
|---|---|---|
| Truth-seeded B0 CKF | Already wired into this branch | Diagnostic reference for losses after seeding. Use identical CKF cuts, material and truth acceptance. It is not usable on real data or an absolute efficiency ceiling. |
| ACTS `estimateTrackParamsFromSeed` | Installed in ACTS 44.4.0 | Use three ordered B0 hits and the actual field vector. Compare local direction and signed q/p first. Requires an adapter, covariance treatment, and transport to the origin reference expected by EICrecon's CKF. |
| ACTS multipoint estimator | Confirmed in release 47.7.0; absent from installed 44.4.0 | `estimateTrackParamsFromSpacePoints` uses all three or four ordered points, optional weights and geometric refinement. It assumes a homogeneous field and supplies parameters, not fitted covariance. Test in an isolated newer environment or review a backport. |
| ACTS fixed-hit Kalman fitter | Installed in 44.4.0 | Fit each candidate's assigned measurements using the field and material map, reusing EICrecon geometry/source links. Requires an adapter and initial parameters; isolates fitting from CKF hit selection. |
| ACTS global chi-square fitter | Library header installed; no EICrecon factory wiring found | Lower-priority assigned-hit cross-check: this version explicitly leaves energy-loss handling unfinished. Requires integration and does not find candidates. |
| FASER circle / three-station seed tools | Public source in Calypso | Athena/Gaudi tools require substantial porting. The active circle path consumes station tracklets and derived pseudo-points, which B0's individual layer hits do not supply. Prefer the ACTS multipoint kernel for a circle-fit comparison. |
| GENFIT Kalman/DAF/GBL | External toolkit, not wired into this chain | Possible independent fitter, but needs measurement, geometry and field adapters. Lower priority than using installed ACTS components. Its repository reports that active development has ended. |

The installed ACTS three-point estimator explicitly accepts a field vector
along any direction and estimates the state at the first point. Its separate
covariance initializer uses configured diagonal uncertainties and inflation;
it does not automatically supply a full measurement-derived B0 covariance.
[ACTS 44.4.0 estimator API](https://raw.githubusercontent.com/acts-project/acts/v44.4.0/Core/include/Acts/Seeding/EstimateTrackParamsFromSeed.hpp).

This must not be confused with stock ACTS triplet *finding*, whose documented
assumptions include a nearly homogeneous field, field-aligned coordinates,
radial ordering and central origins. Rotating the coordinates alone does not
validate those assumptions across the B0 field boundary. EICrecon's own
`TrackSeeding` parameter calculation also assumes Bz.
[ACTS seeding documentation](https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html).

Other option references:
[47.7.0 multipoint API](https://raw.githubusercontent.com/acts-project/acts/v47.7.0/Core/include/Acts/Seeding/EstimateTrackParamsFromSeed.hpp),
[44.4.0 Kalman fitter](https://raw.githubusercontent.com/acts-project/acts/v44.4.0/Core/include/Acts/TrackFitting/KalmanFitter.hpp),
[44.4.0 global chi-square implementation](https://raw.githubusercontent.com/acts-project/acts/v44.4.0/Core/include/Acts/TrackFitting/GlobalChiSquareFitter.hpp),
[FASER circle tool](https://gitlab.cern.ch/faser/calypso/-/blob/2fab73eb3ff7e08f52bbeebe2c44c1ce16056abe/Tracking/Acts/FaserActsKalmanFilter/src/CircleFitTrackSeedTool.cxx),
[ACTS fitters](https://acts.readthedocs.io/en/latest/core/reconstruction/track_fitting.html),
[GENFIT repository](https://github.com/GenFit/GenFit).

LHCb Rec and FASER Calypso are publicly readable source repositories. Rec's
`COPYING` and Hybrid Seeding source specify GPLv3; Calypso's top-level license
states Apache 2.0 except where other licenses are indicated. Public source
does not make either framework a drop-in B0 package.
[LHCb license](https://gitlab.cern.ch/lhcb/Rec/-/blob/f095bfb07716fccd4375609725ece065db334faa/COPYING),
[FASER license](https://gitlab.cern.ch/faser/calypso/-/blob/2fab73eb3ff7e08f52bbeebe2c44c1ce16056abe/LICENSE).

## Validation priorities

1. Compare the current estimator and installed ACTS estimator on identical
   candidates and field samples. For four-station candidates, choose a
   deterministic wide-lever-arm triplet and reserve the fourth hit as a
   prediction check. Compare states on the same surface before comparing CKF
   output. A field-integral or numerical-transport estimator is a later
   controlled comparison if field-dependent bias is observed.
2. Use truth hit associations, not the seed with the closest truth momentum,
   to measure seed efficiency and pulls. Define the denominator as primary
   protons with hits in at least three distinct stations, and report any
   additional momentum, angle or cleanliness cuts. Report fakes and
   duplicates with their own explicit denominators.
3. Check residual bias, per-object pulls, covariance correlations and tails
   versus momentum, angle, station count and vertex position. Separate gun
   vertices from the afterburned beam spot. The default beam sizes and
   scattering terms are sample/geometry calibrations, not universal values.
4. Keep geometry, material, digitization, CKF cuts and ambiguity selection
   fixed when comparing seeders. The branch's chi-square cutoff of 50 remains
   provisional without background validation. Measure candidate-budget
   saturation, runtime and fake rate under representative occupancy before
   replacing enumeration with a more elaborate finder.
5. Validate displaced decays separately. The origin constraint and current
   CKF seed-reference convention need explicit treatment; switching off the
   constraint alone is not evidence of a production-ready displaced mode.

## Three-point trial

The standalone experiment is published on the fork's
[`study/b0-acts-three-point` branch](https://github.com/tom-bleher/EICrecon/tree/study/b0-acts-three-point/src/benchmarks/b0_acts_three_point).
It calls the installed ACTS estimator on the first, nearest-midpoint and last
hits of each persisted candidate. All four hits remain available to the
current estimator. Both use a fitted-path midpoint field sample; this trial
uses equal-weight fits, valid for the uniform hit variances in this sample.
Direction comparisons transport the current origin seed numerically to the
first hit's global-z plane through the field, without material interactions.
This is a seed-parameter comparison, not a replacement factory or CKF trial.

On 5,000 afterburned DVCS events with the geometry above, 3,253 unambiguous
primary-proton-associated candidates yielded finite ACTS estimates. Select
the first emitted candidate per truth particle, never the closest to truth.
For the 2,288 particles whose selected hits all retain over half the origin
momentum and have positive lab pz, q/p residual central widths were **3.79%**
(current) and **4.16%** (ACTS), with median biases **+0.22%** and **−0.45%**.
Both residuals use truth momentum at the first simulated hit. Width is
half the 84th-minus-16th percentile span, rounded here to two decimals.
A paired bootstrap put the ACTS-minus-current width difference at
**+0.22 to +0.55 percentage points** (95% interval).

The ACTS local projected-slope residual widths were smaller: x/y **0.166/0.092**
versus **0.347/0.268** in units of 10^-3. These are slope differences, not
angular distances. The clean three-hit subset had only 167 particles
(q/p widths 7.00% versus 6.91%); it does not establish an improvement.
These results motivate an all-four-hit stock estimator test, while retaining
the existing candidate finder and treating covariance/CKF integration as
separate work.

## Exact regression comparison

Use [compare_tracking_outputs.py](../../src/scripts/compare_tracking_outputs.py)
inside eic-shell:

```bash
python3 src/scripts/compare_tracking_outputs.py before.root after.root
python3 -m unittest discover -s src/tests/scripts_test -v
```

The default scope is all persisted B0 collections and their recursive relation
closure. Persist `EventHeader` and all referenced collections. The comparison
requires unique run/event identities and identical event sets, checks every
event and primitive field without rounding, and preserves object ordering and
relation indices. It permits event-entry reordering. A pass concerns selected
event data, not ROOT file bytes or unselected collections.

The earlier 5,000-event validation files omitted `EventHeader` and referenced
collections. They cannot establish this exact comparison. The previous
rounded-q/p, 2,000-entry check does not substantiate a bit-identical-output
claim; fresh output with complete identity and persistence is required.

The fallback fix was checked first on 100 events, then on 5,000 events of the
existing afterburned `DVCS.1.ab.hiAcc.5x41` simulation with the configuration
above. All 257 fields in the selected B0 relation closure match the preceding
revision exactly; no event or analysis cuts were applied to the comparison.
The same 5,000 events also match exactly between one and eight threads, and
the relation checker reports no missing referenced collections.
The missing-error behavior is covered separately by the B0 regression test.
Build/install, both CTest suites, 21 B0 test cases (469 assertions), 10 Python
comparator tests, pinned clang-format 23.1.0 and targeted clang-tidy passed.

## Joint covariance correction

The non-bending correction subtracts a term proportional to fitted q/p from
each y coordinate. Treating its subsequent line fit as independent of the
bending fit omitted momentum-induced y uncertainty and cross correlations.
The corrected implementation applies the weighted line estimator to
`dy_i/d(q/p)`, then propagates the joint entrance covariance before converting
to perigee parameters. Field samples remain fixed in this covariance model.
Seed values and candidate selection are unchanged; the covariance supplied
to CKF changes. The old detached `prehard` comparison worktree also has the
omission; it remains a historical baseline rather than being edited.

An independent original-hit finite-difference regression checks all covariance
elements for 24 combinations of 3/4 hits, both charge signs, positive/zero/
negative Bx and constrained/unconstrained seeds, using unequal hit variances.
Full build/install and both CTest suites passed; B0 now has 22 test cases and
1,117 assertions. Formatting, targeted clang-tidy and all ten comparator
tests passed. Fresh 100-event checks preceded 5,000-event DVCS and 1,000-event
38 GeV proton-gun checks. The DVCS output matches exactly across one/eight
threads over all 257 fields; all 29 relation branches pass validation.

The truth-associated covariance audit on the same DVCS sample gives:

| Selection | Particles | Joint 95% Gaussian-region coverage |
|---|---:|---:|
| First emitted seed per matched primary proton | 2,389 | 90.00% |
| Additionally clean forward hits | 2,288 | 93.53% |
| Clean, three hits | 167 | 82.04% |
| Clean, four hits | 2,121 | 94.44% |

Covariance pulls use origin truth and the full five-dimensional perigee
covariance; the clean selection is defined in the trial above. All 3,253
associated seed covariances are finite and positive definite. The first-per-
truth phi pull width moves from 1.095 to 1.020, and clean joint coverage from
92.66% to 93.53%. Percentages are rounded to two decimals and pull widths to
three. Remaining non-Gaussian tails, especially in three-hit DVCS seeds,
are **not resolved** by this correction. Do not globally inflate covariance
to fit this one sample; validate phase-space dependence and background first.

The independent gun input is `pg_devepic_p38/pgun.theta4to22` from the same
20260825 simulation production and geometry. Of 1,000 events, 676 candidates
have unambiguous primary-proton associations; 447 first clean matches give
94.85% joint coverage. Its truth vertex is fixed, so the default beam-spot
contribution overstates parts of its uncertainty (theta pull width 0.062).
That result cannot establish calibration for the afterburned DVCS beam spot.
The reusable extraction and tail audit are included with the study branch.

The default CKF output changes from 2,285 to 2,283 tracks: events 2338, 3578
and 4691 lose a track, while 4392 gains one. Seed states and selected hit
triplets for all 3,253 associated candidates are unchanged exactly. A
diagnostic rerun with `tracking:B0TrackerCKFTrajectories:Chi2CutOff=100`
recovers all three lost tracks, consistent with covariance-dependent CKF
measurement gating; their fitted chi-square/ndf values are 116.64/12,
108.21/8 and 102.92/6. The production cutoff remains 50. This small count
change is not an efficiency improvement claim, and the looser diagnostic
cut is not a recommended default.
