# B0 seeding: implementation and comparison options

Research checked 6 September 2026 against `feat/b0-stub-seeder` at
`1cab33ff`, upstream EICrecon `cef494fd`, and the installed ACTS 44.4.0 in
eic-shell `26.04.0-stable`. The target is the realistic B0 geometry at epic
`2200b0aa`, using `epic_ip6_extended_5x41` and 5 × 41 GeV beams.

## Assessment

Recommendation (inferred from the source and comparisons below): keep the
dedicated B0 candidate finder and the existing ACTS CKF. The most
useful additional experiment is the installed ACTS three-point parameter
estimator, applied to the same candidate hits. No turnkey B0 replacement was
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
| ACTS multipoint estimator | Present in current ACTS main; absent from this installation | `estimateTrackParamsFromSpacePoints` offers weighted circle/line estimation with a homogeneous-field assumption. Try in an isolated newer environment after the three-point trial; no release availability is asserted here. |
| ACTS global chi-square fitter | Library header installed; no EICrecon factory wiring found | Independent fit/covariance cross-check on assigned hits. Requires integration and does not find the initial candidates. |
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
[current multipoint API](https://raw.githubusercontent.com/acts-project/acts/main/Core/include/Acts/Seeding/EstimateTrackParamsFromSeed.hpp),
[ACTS fitters](https://acts.readthedocs.io/en/latest/core/reconstruction/track_fitting.html),
[GENFIT repository](https://github.com/GenFit/GenFit).

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

Performance gains from the alternatives are unverified. These are recommended
experiments, not measured improvements from a stock
replacement. The implementation changes accompanying this note fix covariance
fallback handling and exact output validation only.

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
