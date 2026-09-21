// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

namespace eicrecon {

/// Configuration for the B0 stub seeder.
///
/// The B0 tracker sits inside the B0pf combined-function magnet (dipole along
/// y plus a quadrupole gradient), where the solenoid-model helix estimation of
/// the orthogonal seeder returns invalid parameters (random charge sign). This
/// seeder instead fits the hits of the four B0 stations directly: a parabola
/// in the bend plane gives direction and signed q/p from the local dipole
/// field, the quadrupole term Bx = G*y is removed from the non-bend plane, and
/// the result is back-extrapolated analytically to the origin perigee that
/// CKFTracking expects.
struct B0TrackerStubSeederConfig {

  // --- geometry / field ---
  /// Minimum |By| [T] at the candidate's middle hit, sampled from the
  /// ACTS/DD4hep field provider. A guard against a zero-field query only.
  float minAbsFieldY = 0.05;
  /// Ion-frame z of the B0pf field entrance [mm]. Used for the
  /// origin-constrained upstream chord and as the bend-plane reference.
  /// <= 0 (the default) derives it from B0PF_CenterPosition, B0PF_Length and
  /// B0PF_XPosition. A positive value overrides that for diagnostics.
  float zFieldEntrance = 0.0;

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
  /// Maximum RMS residual of the y(z) straight-line fit [mm] after the
  /// quadrupole (Bx) bending has been removed. This rejects cross-station
  /// combinations that cannot be one telescope trajectory.
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

  // --- interaction-vertex spread ---
  // Gaussian sigma of the interaction vertex in the lab frame [mm]. The
  // beamline-constrained seed reports the parameters of a track leaving the
  // origin, so the spread of the true vertex is an uncertainty on those
  // parameters, and the longitudinal term dominates the polar angle:
  // sigma(theta) = theta * sigmaZ / 6 m, about 180 urad for a 30 mrad track
  // and a 37 mm bunch, against 13 urad from the hit resolution.
  // The defaults are the vertex distribution of the afterburned 5x41 DVCS
  // sample (RMS 0.168, 0.019, 37.2 mm; the crossing angle leaves it diagonal
  // in the lab frame to better than 0.07 in correlation). Beam settings
  // change them. Set all three to zero for a gun thrown from a fixed point.
  float beamSpotSizeX = 0.17;
  float beamSpotSizeY = 0.02;
  float beamSpotSizeZ = 37.0;

  // --- covariance fallbacks (diagonal variances, edm4eic units: mm^2,
  // rad^2, (1/GeV)^2, ns^2) ---
  // The fitted hit covariance is propagated to loc/angles/q/p, including
  // correlations, and the vertex spread above is added. These values are used
  // only for parameters that neither constrains, or when input hit variances
  // are unavailable; with a non-zero beam spot that no longer includes
  // loc0/loc1 in beamline-constrained mode.
  float locaVariance   = 4.0;
  float locbVariance   = 1600.0;
  float phiVariance    = 0.01;
  float thetaVariance  = 1.0e-5;
  float qOverPVariance = 2.0e-4;
  /// The seed time is the time at the origin perigee, i.e. the vertex time
  /// (0 ns for prompt tracks); this variance keeps it unconstrained until B0
  /// has a timing digitization.
  float timeVariance = 100.0;

  // --- scattering terms added to the propagated hit covariance ---
  // The propagated measurement covariance describes the hit resolution only.
  // Multiple scattering in the beam pipe and the B0 stations, and the residual
  // trajectory-model error, scale with 1/p and are added here. Measured on a
  // proton-gun sample (3000 protons, 8-41 GeV, ion-frame 4-22 mrad,
  // epic_ip6_extended with the realistic B0 modules and a material map
  // generated for that geometry): the seed direction is good to 13-24 urad
  // and q/p to 2.4-4.3 %, with pull widths of 0.85-1.0 after these terms.
  // CKFTracking then reconstructs 698 of 765 seeded four-station events
  // (723 with truth seeds); a wider window gains under 2 % of tracks and
  // degrades the fitted momentum. NOTE: with a material map that does not
  // match the geometry the CKF becomes strongly window-limited; regenerate
  // the map (epic scripts/refresh_local_material_map.sh) rather than
  // inflating these.
  /// Scattering angle coefficient [GeV rad]: sigma = scale * |q/p| is added
  /// to theta and, divided by sin(theta), to phi.
  float scatteringScale = 3.0e-4;
  /// Relative q/p uncertainty from scattering and the trajectory model.
  float qOverPRelativeUncertainty = 0.025;
};

} // namespace eicrecon
