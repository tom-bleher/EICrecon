// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <cstddef>

#include "algorithms/tracking/HelixSeed.h"

namespace eicrecon {

/// Configuration of TruthTrackSeeding.
struct TruthTrackSeedingConfig {
  /// Minimum distance from the last hit kept for a particle to the next one, in time order, so
  /// that the hits a particle is seeded from lie on different layers or stations
  double minHitDistance = 10. * Acts::UnitConstants::mm;

  /// Maximum number of hits per particle, the first in time order; every triplet of them is a seed
  std::size_t maxHitsPerParticle = 6;

  /// Where the seeds start and their prior
  HelixSeedPrior prior{};
};

} // namespace eicrecon
