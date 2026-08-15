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
/// parabola in the bend plane, giving direction, momentum and charge, then
/// back-extrapolates analytically to the origin perigee that CKFTracking
/// expects.
struct B0TrackerStubSeederConfig {

  // --- geometry / field ---
  /// Ion-beam rotation about y (rad); B0 hits are fitted in this rotated frame
  float crossingAngle = -0.025;
  /// The dipole field is sampled from the ACTS/DD4hep field provider along
  /// each fitted candidate. This is only a guard against a zero-field query,
  /// not an analytic B0 field substitute.
  float minAbsFieldY = 0.05;
  /// Ion-frame z where the B0pf field region starts (mm), used only for the
  /// origin-constrained upstream chord. Momentum always uses the sampled ACTS
  /// field at the fitted candidate.
  float zFieldEntrance = 5800.0;
  /// Number of ACTS field samples between the innermost and outermost selected
  /// B0 stations. A single sample is allowed for a uniform-field test.
  unsigned int fieldSamples = 5;

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
  /// Cap on hit combinations tried per event
  unsigned int maxCombinations = 512;
  /// Maximum seeds emitted per event (after overlap deduplication)
  unsigned int maxSeeds = 20;
  /// Two accepted seeds may share at most this many hits. Keeping this below
  /// three prevents leave-one-out subsets from duplicating their parent seed.
  unsigned int maxSharedHits = 2;
  /// Charge policy: 0 infers the sign from fitted curvature and sampled By;
  /// +1 or -1 forces a diagnostic hypothesis.
  int charge = 0;
  /// Emit both charge hypotheses for every compatible candidate. Intended for
  /// validation of opposite-sign tracks; normal reconstruction infers one.
  bool testBothCharges = false;
  /// Accepted range for the field-aware fitted momentum [GeV]. Candidates
  /// outside this range are rejected rather than assigned a nominal momentum.
  float pMin = 3.0;
  float pMax = 400.0;

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

  // --- seed covariance (diagonal variances, edm4eic units: mm^2, rad^2,
  // (1/GeV)^2, ns^2), same storage convention as OrthogonalTrackSeedingConfig ---
  /// d0 variance: beam-spot scale (constrained) rather than fit-extrapolation
  float locaError = 4.0;
  /// z0 variance: bunch-length scale (~35 mm RMS)
  float locbError = 1600.0;
  /// phi variance: transverse-slope error / small |t_xy| -> sigma ~ 0.1 rad
  float phiError = 0.01;
  /// theta variance: chord + field-model systematics -> sigma ~ 3 mrad
  float thetaError = 1.0e-5;
  /// q/p variance: quad-gradient bias on the sagitta -> sigma ~ 50% of q/p
  float qOverPError = 2.0e-4;
  float timeError   = 100.0;
};

} // namespace eicrecon
