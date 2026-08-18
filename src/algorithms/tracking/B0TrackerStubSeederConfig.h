// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

namespace eicrecon {

/// Configuration for the B0 stub seeder.
///
/// The B0 tracker sits inside the B0pf dipole (field along y), where the
/// solenoid-model helix estimation of the orthogonal seeder returns invalid
/// parameters (random charge sign). This seeder instead fits the hits of the
/// four B0 stations directly: a straight line in the non-bend plane and a
/// sampled field-integral model in the bend plane, giving direction and signed q/p, then
/// back-extrapolates analytically to the origin perigee that CKFTracking
/// expects.
struct B0TrackerStubSeederConfig {

  // --- geometry / field ---
  /// The dipole field is sampled from the ACTS/DD4hep field provider along
  /// each fitted candidate. This is only a guard against a zero-field query,
  /// not an analytic B0 field substitute.
  float minAbsFieldY = 0.05;
  /// Ion-frame z of the B0pf field entrance [mm]. Used for the
  /// origin-constrained upstream chord and as the field-integral reference.
  /// <= 0 (the default) derives it from B0PF_CenterPosition, B0PF_Length and
  /// B0PF_XPosition. A positive value overrides that for diagnostics.
  float zFieldEntrance = 0.0;
  /// Minimum number of ACTS field mesh points from the B0pf entrance through
  /// the outermost selected station; all selected hit positions are added too.
  unsigned int fieldSamples = 5;
  /// Fit once on the parabola reference path, then repeat on the fitted path.
  unsigned int fieldFitIterations = 2;

  // --- seeding logic ---
  /// Hits whose ion-frame z differs by more than this (mm) are different
  /// stations. Official B0 disks are ~270 mm apart; realistic front/back
  /// sensors on one disk are ~7 mm. 50 mm splits the first and merges the
  /// second, so the same seeder works for 4- and 8-layer geometries.
  float stationZGap = 50.0;
  /// Minimum number of distinct stations in a seed candidate
  unsigned int minStations = 3;
  /// Half-width (mm) of the non-bend-plane road used to preselect hits before
  /// the fit. The dipole does not bend in y, so the hits of one track are
  /// collinear in y to well within this width, while hits of different tracks
  /// generally are not. This keeps the enumeration below maxCombinations in
  /// busy events instead of truncating it.
  float yRoadWidth = 5.0;
  /// Maximum non-bending residual [mm] of an endpoint pair relative to a ray
  /// from the beamline. Applied only when constrainToBeamline is true; <=0
  /// disables it. Compatible pairs are sorted by this residual before the
  /// combination budget is spent.
  float maxYBeamlineResidual = 5.0;
  /// Cap on hit combinations tried per event
  unsigned int maxCombinations = 512;
  /// Maximum seeds emitted per event (after overlap deduplication)
  unsigned int maxSeeds = 20;
  /// Two accepted seeds may share at most this many hits. Keeping this below
  /// three prevents leave-one-out subsets from duplicating their parent seed.
  unsigned int maxSharedHits = 2;
  /// Hits at the same station closer than this (mm, ion-frame transverse
  /// distance) count as the same hit for the sharing test above. Front/back
  /// sensors of one disk see the same track ~0.1 mm apart, so without this a
  /// single proton with one doubled station emits four seeds, with two doubled
  /// stations up to sixteen. 0 falls back to exact hit identity.
  float sharedHitDistance = 0.5;
  /// Charge policy: 0 infers the sign from fitted curvature and sampled By;
  /// +1 or -1 forces a diagnostic hypothesis.
  int charge = 0;
  /// Emit both charge hypotheses for every compatible candidate. Intended for
  /// validation of opposite-sign tracks; normal reconstruction infers one.
  bool testBothCharges = false;
  /// Minimum accepted momentum [GeV]. This rejects clearly incompatible
  /// low-momentum curvature; unresolved high momentum is retained as q/p~0.
  float pMin = 3.0;
  /// Below this |q/p| significance, emit both charge hypotheses because the
  /// curvature sign is unresolved by the fitted hit precision.
  float minCurvatureSignificance = 3.0;

  // --- compatibility and ranking ---
  /// Maximum RMS residual of the y(z) straight-line fit [mm]. This rejects
  /// cross-station combinations that cannot be one telescope trajectory.
  float maxYResidual = 0.5;
  /// Maximum RMS residual of the x(z) parabola fit [mm].
  float maxXResidual = 2.0;
  /// Maximum absolute transverse slope at the first selected station.
  float maxAbsTransverseSlope = 0.10;

  /// Constrain the seed to the beamline: direction from the chord
  /// origin -> field-entrance point, loc = (0,0) with beam-spot covariance.
  /// Kills the 6 m lever-arm amplification of fit errors into the perigee.
  /// NOTE: false is a diagnostic mode only -- the loc covariances below are
  /// tuned for the constrained mode and grossly underestimate the
  /// extrapolated perigee uncertainty when this is disabled. Reconstruction
  /// of displaced decays needs a seed anchored away from the origin, which
  /// CKFTracking does not yet support.
  bool constrainToBeamline = true;

  // --- covariance fallbacks (diagonal variances, edm4eic units: mm^2,
  // rad^2, (1/GeV)^2, ns^2) ---
  // The fitted hit covariance is propagated to loc/angles/q/p, including
  // correlations. These values are used only for parameters that the fit does
  // not constrain (notably loc0/loc1 in beamline-constrained mode), or when
  // input hit variances are unavailable.
  float locaVariance   = 4.0;
  float locbVariance   = 1600.0;
  float phiVariance    = 0.01;
  float thetaVariance  = 1.0e-5;
  float qOverPVariance = 2.0e-4;
  float timeVariance   = 100.0;

  /// Seed time prior [ns], paired with the deliberately loose timeVariance
  /// above. B0 has no timing digitization yet -- the generic silicon front end
  /// only smears the Geant4 time -- so this is an explicit prior rather than an
  /// estimate. Replace with a hit-derived, time-of-flight-corrected estimator
  /// once a real AC-LGAD timing response exists.
  float seedTime = 10.0;

  // --- calibrated field/material/model additions ---
  // Derived from 5k-event proton-gun samples at 15, 25, 35, and 41 GeV and
  // validated at 20, 30, and 38 GeV. These diagonal additions supplement,
  // rather than replace, the propagated hit-measurement covariance.
  //
  // Provenance of the values below -- regenerate with
  // presentation16/calibrate_b0_seed_covariance.py before changing them:
  //   calibration artifact : presentation16/assets/b0_seed_covariance/calibration.json
  //                          sha256 e039bbc0d9378776cb7c757873b459fcaad1ae275a471dd
  //                                 23038104e834595e8
  //   reconstruction run   : presentation16/b0_covariance_validation_20260817
  //   simulation input     : presentation14/run_p14dev_pcurve_20260814_134439
  //   geometry             : presentation14/p14dev_detector_20260814_124519
  //                          (epic @ b0-tracker-realistic-geometry, epic_ip6_extended)
  //   material map         : pg_p14dev_p41_6b9ec99108/material-map.cbor
  //                          sha256 74a4b052e1154b5b2babe708944ec04fe838faafd83b79a
  //                                 8d7c53716c6e654a3
  // The constants are only valid for that geometry/field/material combination;
  // a geometry change invalidates them and requires a re-calibration.
  /// Momentum-independent residual phi variance [rad^2].
  float phiModelVariance = 1.934e-5;
  /// Coefficient of the material/model term sigma(phi) = scale * |q/p| [GeV rad].
  float phiQOverPScale = 0.1360;
  /// Momentum-independent residual theta variance [rad^2].
  float thetaModelVariance = 2.201e-9;
  /// Coefficient of the material/model term sigma(theta) = scale * |q/p| [GeV rad].
  float thetaQOverPScale = 1.043e-3;
  /// Momentum-independent residual q/p variance [(1/GeV)^2].
  float qOverPModelVariance = 6.290e-8;
  /// Relative residual q/p uncertainty from the trajectory model.
  float qOverPRelativeUncertainty = 0.02253;
  /// External relative B-field integral uncertainty. The calibration samples
  /// use the same map in simulation and reconstruction, so this is not
  /// identifiable there and deliberately defaults to zero.
  float fieldRelativeUncertainty = 0.0;
};

} // namespace eicrecon
