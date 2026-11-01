// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <Acts/Definitions/Units.hpp>

#include <vector>

namespace eicrecon {
struct CKFTrackingConfig {
  std::vector<double> etaBins                    = {};
  std::vector<double> chi2CutOff                 = {15.};
  std::vector<std::size_t> numMeasurementsCutOff = {10};

  std::size_t numMeasurementsMin = 4;

  // --- Per-branch stopping for the CKF ---
  // A branch whose fit has diverged keeps inflating its covariance until it
  // overflows double precision and produces a NaN chi2. ACTS's default branch
  // stopper never stops, so nothing contains this. A branch is dropped once
  // the filtered variance of q/p exceeds this limit, or is no longer positive.
  // The default sits midway between the variance a seed asserts (0.1) and what
  // a diverged fit reaches (1e10 and up); it corresponds to a sigma(q/p) that
  // is already many times |q/p| for any track worth keeping.
  // See the implementation in CKFTracking.cc.
  double maxQOverPVariance = 1e5 / (Acts::UnitConstants::GeV * Acts::UnitConstants::GeV);

  /// Include the measurement time (in addition to loc0 and loc1) in the calibrated measurements
  bool useTime = false;

  /// Extrapolate tracks to the perigee at the origin in two stages (see TrackExtrapolation.h):
  /// from the first measurement backward to a plane across the beam at referencePlaneZ, with the
  /// material of the tracking geometry, then from there to the perigee without material. Tracks
  /// that cannot reach the perigee from the plane are kept on the plane. For a spectrometer far
  /// from the interaction point, such as the B0 tracker. Off by default, which extrapolates with
  /// Acts::extrapolateTrackToReferenceSurface.
  bool twoStageExtrapolation = false;
  double referencePlaneZ     = 0.;
  /// Tolerance on reaching the perigee in stage 2. The Acts default (1e-4) makes a perigee line
  /// unreachable for tracks within about 10 mrad of it.
  double perigeeSurfaceTolerance = 1e-6 * Acts::UnitConstants::mm;
  /// Path limit of stage 2
  double maxExtrapolationPath = 20. * Acts::UnitConstants::m;
};
} // namespace eicrecon
