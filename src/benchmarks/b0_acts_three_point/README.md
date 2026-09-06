# ACTS three-point B0 seed study

A standalone, offline comparison of the installed ACTS
`estimateTrackParamsFromSeed` against persisted B0 stub seeds. It calls the
actual C++ ACTS implementation and the same `DD4hepFieldAdapter` class used by
EICrecon. It does not change reconstruction factories or run a replacement CKF.

## Run

Run in eic-shell, with a local geometry matching the simulation and
reconstruction. ACTS 44.4.0, DD4hep, Python uproot/awkward/numpy and scipy (for
covariance auditing) were used. This benchmark has its own CMake project;
normal EICrecon builds are unaffected.

```bash
# From the repository root, inside eic-shell:
source /path/to/epic/install/bin/thisepic.sh epic_ip6_extended_5x41
STUDY=$PWD/src/benchmarks/b0_acts_three_point
cmake -S "$STUDY" -B build-b0-study
cmake --build build-b0-study -j$(nproc)
mkdir -p results/b0-study
python3 "$STUDY/extract.py" reconstruction.root results/b0-study/seeds.txt \
  --events 100 --crossing-angle=-0.025
build-b0-study/trial "$DETECTOR_PATH/$DETECTOR_CONFIG.xml" \
  results/b0-study/seeds.txt results/b0-study/estimates.txt
python3 "$STUDY/analyze.py" results/b0-study/seeds.txt \
  results/b0-study/estimates.txt > results/b0-study/comparison.json
python3 "$STUDY/seed_covariance_audit.py" results/b0-study/seeds.txt.json \
  > results/b0-study/covariance.json
python3 "$STUDY/test_audit.py"
```

Use tmux or another detachable session for builds and larger samples. Check
the actual exit code and output. Geometry field-map paths must resolve from
the run directory, as in reconstruction. Preserve the input identity,
EICrecon/geometry revisions, container release, detector configuration and
beam/generator settings alongside results. Do not mix geometries between
reconstruction and this comparison.

The ROOT input must persist `EventHeader`, `MCParticles`, `B0TrackerHits`,
`B0TrackerRawHits`, `B0TrackerRawHitAssociations`, `B0TrackerRecHits`,
`B0TrackerSeeds` and `B0TrackerSeedParameters`. Extraction checks relation
collection IDs, index bounds and unique run/event identities. It requires
ordinary owning collections with this schema. Missing input fails explicitly.

## Comparison definition

- Follow seed → reconstructed hit → raw hit → simulation hit → MCParticle.
  Each hit must have exactly one raw-to-simulation association; every hit in
  the seed must belong to the same primary proton (`PDG == 2212`,
  `generatorStatus == 1`). Report exclusions in the counts JSON.
- Order the persisted candidate hits in the ion frame. Use first and last
  hits and the interior hit nearest their midpoint; ties retain seed order.
  This compares parameter estimation on already accepted candidates, not
  candidate-finding efficiency.
- Query the field on an **unweighted** parabola/line fitted midpoint. This
  approximates the production weighted-fit query; the fits coincide for equal
  per-hit variances. It is not a claim of identical field queries for arbitrary
  sensor resolutions. ACTS uses the full sampled field vector and assumes
  constant field across its triplet.
- ACTS returns a free state at the first hit. Propagate the stored current
  seed from origin perigee to that hit's global-z plane with RK4 and the full
  field, without energy loss or scattering. Both directions are compared to
  the associated simulation-hit momentum there. The 1 mm step is configurable
  through a final C++ argument; compare with 0.5 mm to assess convergence.
  The two estimates can have different transverse positions on that plane;
  this is not a comparison on the actual tilted sensor surface.
- Compare both signed q/p values with truth at that same first sensor. Current
  q/p is transported unchanged because material is omitted. This reference
  choice has a small energy-loss bias relative to the production seed's
  original origin reference. A separate covariance audit uses origin truth.
- Report all associated candidates and the first emitted seed per truth
  particle, without selecting closest q/p or correct charge. The additional
  clean-forward subset requires every associated simulation hit to have
  positive lab pz and retain more than half its particle's origin momentum.
  This removes strongly degraded/re-entering hits using truth and is a
  diagnostic subset, not a reconstruction cut or efficiency denominator.

Widths are half the 16th–84th percentile interval, using NumPy's default
linear quantile interpolation. q/p residuals are percentages; the direction
separation is an angle in mrad. The projected-slope residual fields contain
1,000 times the difference in dx/dz or dy/dz and are small-angle proxies, not
angular distances.
The paired bootstrap uses 1,000 resamples and a fixed random seed; its interval
is conditional on the clean selected particles, not a systematic uncertainty.
Raw numerical outputs retain double precision; no truth-momentum rounding or
truth-based seed ranking is performed.

The executable checks exact circles for both charge signs in axial and rotated
fields before loading geometry. Covariance auditing reports marginal pulls,
positive definiteness and joint five-dimensional coverage. Its origin-truth
reference assumes the upstream path is field-free and is inappropriate for
samples with substantial upstream bending.

## Interpretation limits

The study does not provide an ACTS covariance adapter, origin transport with
material, end-to-end CKF efficiency, fake-rate comparison, or background
validation. Better first-hit direction residuals do not by themselves establish
better reconstruction. A successful run establishes that this stock parameter
estimator can consume B0 triplets, not that it is a turnkey seeder replacement.
No large input/output ROOT files or machine-specific paths are stored here.

## Recorded baseline study (6 September 2026)

These results predate the later covariance correction and identify the exact
baseline, rather than claiming to describe every future branch revision:

- EICrecon `f97b79c14d9a0566d1914ac2b5364b89b599d812` (seed calculation from
  `cf40fc36`); eic-shell `26.04.0-stable`; ACTS `44.4.0`.
- Geometry `2200b0aa76ed398e3e3d752e9d5346b43e171405`,
  `epic_ip6_extended_5x41`, electron 5 GeV × proton 41 GeV.
- Afterburned `DVCS.1.ab.hiAcc.5x41`, first 5,000 events from the local
  `20260825_030157/combined/dvcs_devepic` simulation sample. This is a local
  sample identifier, not a public campaign download URL.
- Entry XML SHA256:
  `611f6fe1fbef6f52d4c1ab05ab2268d1cb7327cde08e40552c3c7385194eb94f`.
  Reconstruction material-map SHA256:
  `62cbf103e6c41b77397ddba591c53664359172d486da886702c6eca0da2c77b7`.
- Default branch seed/CKF settings, beam widths (0.17, 0.02, 37) mm; eight
  reconstruction threads. The offline propagation includes field only.

Of 3,310 emitted seeds, extraction excludes 48 with nonunique hit truth
associations and nine belonging to other particles. The remaining 3,253 seeds
represent 2,389 distinct primary protons. The clean first-per-particle subset
contains 2,288 protons (167 three-hit and 2,121 four-hit seeds).

| Clean first-per-particle metric | Current stub | ACTS three-point |
|---|---:|---:|
| Signed q/p residual median | +0.22% | −0.45% |
| Signed q/p residual half-68% width | 3.79% | 4.16% |
| First-hit direction error median | 0.350 mrad | 0.158 mrad |
| First-hit direction error 95th percentile | 1.005 mrad | 0.405 mrad |

Table percentages are rounded to two decimal places and angles to three.
The paired-bootstrap 95% interval for ACTS minus current q/p width is
[+0.22, +0.55] percentage points. This sample supports better local direction
estimates and slightly worse momentum resolution for ACTS. For the 167 clean
three-hit seeds the q/p widths are 7.00% and 6.91%, respectively; the overall
momentum difference is dominated by four-hit candidates, where ACTS discards
one measurement.

The unfiltered first-per-particle sample has large tails: q/p RMS residuals
are 48.34% (current) and 46.37% (ACTS), with direction-error RMS near 96 mrad
for both. The clean-subset improvements do not remove these failures. These
are seed-level conditional measurements, with no end-to-end CKF claim.

The published harness was rebuilt and rerun on the first 100 events: its 90
associated seed rows reproduce the original fitted-midpoint trial's numerical
input/output and comparison. Four synthetic circle sign/rotation checks also
pass. Large validation artifacts remain external to the repository.

Covariance audit coverage uses **all selected seeds** as denominator, including
invalid covariance/residual cases as failures. Marginal quantiles and joint
chi-square medians use only valid entries, with valid/invalid counts reported.
Nonfinite and nonpositive diagonals are excluded before eigendecomposition;
invalid covariance categories may overlap. Undefined quantiles are JSON null.
