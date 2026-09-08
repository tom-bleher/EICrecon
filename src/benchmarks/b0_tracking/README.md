# B0 tracking validation

Run inside the EIC environment used to build reconstruction. `analyze.py` requires
uproot, awkward and NumPy. It checks direct-primary proton associations through
measurement → reconstructed hit → raw hit → simulated hit, including collection
IDs, bounds, cell IDs and digitized charge closure. It never matches by proximity.

```bash
python -m unittest discover -s src/benchmarks/b0_tracking
python src/benchmarks/b0_tracking/analyze.py baseline.root \
  --variant refit.root --threshold-keV 0.54 --provenance provenance.json \
  --output comparison.json
# Also compare candidates before ambiguity resolution:
python src/benchmarks/b0_tracking/analyze.py baseline.root \
  --variant refit.root --suffix Unfiltered --threshold-keV 0.54 \
  --provenance provenance.json --output unfiltered.json
```

Persist `EventHeader`, `MCParticles`, `B0TrackerHits`, `B0TrackerRawHits`,
`B0TrackerRawHitAssociations`, `B0TrackerRecHits`, `B0TrackerMeasurements`,
`B0TrackerSeeds`, `B0TrackerSeedParameters`, and the chosen
`B0TrackerCKF{Tracks,Trajectories,TrackParameters}` collections, plus their
recursive relation closure. For unfiltered analysis append `Unfiltered` to each
of those three CKF collection names. For the truth-seeded control use
`--prefix B0TrackerCKFTruthSeeded` with its corresponding collections and closure.
Before a refit comparison, verify identical seed priors and measurements with
`src/scripts/compare_tracking_outputs.py --collections B0TrackerSeeds,B0TrackerMeasurements`;
matching seed ObjectIDs alone does not prove identical prior contents.
The threshold must match digitization; a mismatch is an error. Current geometry
has two CAD layers per physical station; `--layers-per-station 1` supports the
older single-layer layout explicitly.

Efficiency is matched generator-status-1 protons divided by such protons with
**direct simulated hits in at least three physical stations**. Every measurement
has one vote shared between active simulated deposits; secondary contributions
remain in the denominator but cannot establish a direct-primary match. Purity
must be at least 0.5. Ties resolve by purity, measurement count, chi-square, then
track index. The reported unmatched-primary fraction is not a general charged
particle fake rate. Covariance coverage uses the full five-dimensional
origin/lab-z perigee and a chi-square(5) 95% boundary. The truth ray assumes a
field-free vertex region. No clean-hit or Gaussian-core cuts remove tails.
Invalid covariance matrices remain counted; storage-roundoff-limited matrices
are identified without clamping. JSON retains full precision; percentiles use
NumPy's linear interpolation.

Each summary and phase-space bin also reports matched-proton momentum bias,
central 68% half-width, wrong-charge counts, and absolute residual tails above
20% and 100%. These are diagnostics, not cuts. Momentum residuals compare
`1/abs(q/p)` with generated proton momentum, so material transport can contribute.
Zero/nonfinite q/p remains counted separately; finite outliers are never clipped.

## Assigned-hit refit experiment

`tracking:B0TrackerCKFTrajectories:RefitSeedCovarianceScale` defaults to **0**
(disabled). Values at least 1 enable a Kalman refit of exactly the CKF-assigned
measurements, with the original seed mean and scaled seed covariance. Failed,
incomplete or invalid refits are rejected. Finding cuts are unchanged, but
refit chi-square can change downstream ambiguity winners. Finite scaling still
reuses a data-dependent prior: it is a sensitivity experiment, not a proof of
statistical independence or calibrated uncertainties.

```bash
# Source the matching local geometry and reconstruction installation first.
# Repeat the same input/parameters for scale 0, 1, 100, 1000 and 10000.
eicrecon -Pjana:nevents=1000 -Pnthreads=8 \
  -Pacts:MaterialMap="$B0_MATERIAL_MAP" \
  -Ptracking:B0TrackerCKFTrajectories:ParticleHypothesisPdg=2212 \
  -Ptracking:B0TrackerCKFTruthSeededTrajectories:ParticleHypothesisPdg=2212 \
  -Ptracking:B0TrackerCKFTrajectories:RefitSeedCovarianceScale=100 \
  -Ppodio:output_collections="$B0_VALIDATION_COLLECTIONS" \
  -Ppodio:output_file=refit.root "$B0_SIM_INPUT"
```

Use zero seed beam-spot widths only for fixed-origin guns, with
`-Ptracking:B0TrackerSeeds:beamSpotSizeX=0` and the analogous Y/Z parameters.
Do not apply that override to afterburned DVCS. The proton hypothesis above is
specific to this proton study; the general reconstruction default remains pion.

The refit explicitly uses the CKF material-aware extrapolator. ACTS 47.7's
[built-in fitter reference extrapolation](https://github.com/acts-project/acts/blob/v47.7.0/Core/include/Acts/TrackFitting/KalmanFitter.hpp#L965-L974)
constructs empty-actor propagation options, so using it directly would omit
material effects between the fitted state and the reference surface.

## ACTS 47.7 validation, 2026-09-08

Pinned provenance and numerical summaries are in `validation-20260908.json`.
Each sample contains the first 1,000 events of the specified simulation chunk;
geometry is `epic_ip6_extended_5x41`, beams 5×41 GeV, crossing angle −25 mrad,
threshold 0.54 keV, CKF chi-square cutoff 50, proton fit hypothesis. The source
baseline is `e2d0391a97b42c46a546f19a2cf19d63e214da1a`. Reconstruction uses the
current local geometry; its physical XML matches the simulation snapshot apart
from material-map references and version metadata.

| Sample | Eligible | Matched, scale 0 | Joint 95% coverage, scale 0 / 100 / 1000 / 10000 |
|---|---:|---:|---:|
| Afterburned DVCS | 472 | 458 | 71.18 / 74.45 / 73.36 / 73.58% |
| 15 GeV proton gun, 4–22 mrad | 438 | 423 | 63.12 / 65.17 / 65.17 / 64.61% |
| 41 GeV proton gun, 4–22 mrad | 443 | 433 | 69.98 / 72.98 / 71.82 / 69.28% |

All matched stored bound covariances in this table are positive definite.
Nevertheless, coverage is substantially below nominal and does not establish a
common improvement plateau. At scale 10000 the 15 GeV sample loses two matches.
For 543 identical unfiltered DVCS seed/hit candidates, scaling 1000→10000 changes
q/p by a median absolute 0.303%, 84th percentile 1.006%, and maximum 55.42%.
Thus ambiguity selection alone cannot explain the instability. **Keep the refit
disabled by default. Final covariance calibration remains unresolved; do not
replace that problem with a global inflation factor.**

The Cartesian covariance fix propagates the full bound covariance using the
surface Jacobian, including displaced-perigee direction dependence. This fixes
missing output information; it cannot calibrate the underlying bound fit.
The B0 material guard requires actual mapped material on both disc approaches
of every B0 layer and rejects prototypes. A runtime test removing one mapped
approach exits with status 1 instead of silently omitting tracking collections.

## Endpoint and occupancy regressions

`make_gun.py` writes deterministic HepMC3 protons in the ion frame rotated by
−25 mrad into lab coordinates. Single-gun defaults are 8 GeV, tx=−0.002,
ty=0.012; use `--momentum 15 --ty 0.020` for the second endpoint case.
`--multiplicity 8` generates eight independent origin protons per event,
uniform in momentum 8–41 GeV, azimuth, and cos(theta) over 4–22 mrad.
This is a combinatoric stress test, **not a beam-background model**.

Simulate 100 events with npsim and random seed 170008, using a minimal steering
file pointing to the recorded geometry. This study used the recorded opt-in
`Geant4TrackerWeightedPlacementAction` for B0 with `HitPositionCombination=2`
and `CollectSingleDeposits=False`; reproduce that action as well as the XML.
Set reconstruction seed beam-spot widths to zero for these origin guns.

| Sample | Eligible | Before endpoint fix | After endpoint fix |
|---|---:|---:|---:|
| 8 GeV, ty=0.012 | 91 | 2 | 89 |
| 15 GeV, ty=0.020 | 77 | 0 | 74 |
| Eight protons/event | 336 | 304 | 325 |

Counts are matched direct-primary protons, not total tracks. On the occupancy
sample, increasing `(maxCombinations,maxSeeds)` from `(512,20)` to `(8192,100)`
also gives 325 matches. Tightening CKF chi-square cutoff from 50 to 15 reduces
that to 314. Keep the existing defaults. The new early road uses a sampled-field
curvature estimate appropriate to the smooth B0 multipole, then reapplies the
original beamline tolerance after signed quadrupole correction. It is not a
rigorous envelope for arbitrary unsampled field structure. Real background
occupancy and phase-space-dependent efficiency remain necessary production
validation.
