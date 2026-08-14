// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <string>

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
  /// DD4hep readout used to decode the B0 layer field.
  std::string readout = "B0TrackerHits";

  // --- geometry / field (defaults match the analytic B0pf description) ---
  /// Ion-beam rotation about y (rad); B0 hits are fitted in this rotated frame
  float crossingAngle = -0.025;
  /// B0pf dipole field along y (T), analytic MultipoleMagnet value for the
  /// nominal 5x41 beamline. NOTE (pre-upstream TODO): the applied field is
  /// B0PF_Bmax * FieldScaleFactor, and FieldScaleFactor is ~0.30 for the
  /// light-ion beamline configs -- this should be read from the DD4hep field
  /// provider (as CKFTracking does) instead of configured.
  float bFieldY = 1.184;
  /// Ion-frame z where the B0pf field region starts (mm); trajectory is
  /// treated as straight upstream of this plane
  float zFieldEntrance = 5800.0;

  // --- seeding logic ---
  /// Hits whose ion-frame z differs by more than this (mm) are different
  /// stations. Official B0 disks are ~270 mm apart; realistic front/back
  /// sensors on one disk are ~7 mm. 50 mm splits the first and merges the
  /// second, so the same seeder works for 4- and 8-layer geometries.
  float stationZGap = 50.0;
  /// Minimum number of distinct stations in a seed candidate
  unsigned int minStations = 3;
  /// Cap on hit combinations tried per event
  unsigned int maxCombinations = 512;
  /// Maximum seeds emitted per event (after overlap deduplication)
  unsigned int maxSeeds = 20;
  /// Two accepted seeds may share at most this many hits. 3 lets a
  /// leave-one-out subset seed coexist with its parent full-station seed.
  unsigned int maxSharedHits = 3;
  /// Assumed charge (outgoing hadron side). The fitted curvature sign is
  /// currently not used to determine charge; tracks bending opposite to this
  /// assumption will not be reconstructed.
  int charge = 1;
  /// Fallback |p| (GeV) when the sagitta fit is unusable
  float momentumPrior = 41.0;
  /// Accept window for the fitted momentum (GeV); outside -> momentumPrior
  float pMin = 3.0;
  float pMax = 400.0;

  /// Constrain the seed to the beamline: direction from the chord
  /// origin -> field-entrance point, loc = (0,0) with beam-spot covariance.
  /// Kills the 6 m lever-arm amplification of fit errors into the perigee.
  /// NOTE: false is a diagnostic mode only -- the loc covariances below are
  /// tuned for the constrained mode and grossly underestimate the
  /// extrapolated perigee uncertainty when this is disabled.
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
