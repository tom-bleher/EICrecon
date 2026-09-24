// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <cstddef>

namespace eicrecon {

struct B0TripletSeedingConfig {
  /// Hits further apart in z than this belong to different stations
  double stationGap = 50. * Acts::UnitConstants::mm;
  /// Maximum distance of the middle hit from the chord of the outer two,
  /// along the local magnetic field, in which a helix does not bend
  double maxResidual = 1. * Acts::UnitConstants::mm;
  /// Minimum momentum of the three-hit estimate
  double minMomentum = 1. * Acts::UnitConstants::GeV;
  /// Distance upstream of the first hit at which the seed is expressed
  double anchorDistance = 10. * Acts::UnitConstants::mm;
  /// Maximum number of seeds per event, smallest residual first
  std::size_t maxSeeds = 1000;

  /// Seed uncertainties at the first hit
  double positionError       = 0.1 * Acts::UnitConstants::mm;
  double angleError          = 2. * Acts::UnitConstants::mrad;
  double qOverPRelativeError = 0.25;
  double timeError           = 10. * Acts::UnitConstants::ns;
};

} // namespace eicrecon
