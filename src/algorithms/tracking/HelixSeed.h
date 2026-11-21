// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/unit_system.h>
#include <array>
#include <cstdint>
#include <span>

namespace eicrecon {

/// Where a seed starts and its prior, a loose diagonal envelope for the track finding
struct HelixSeedPrior {
  /// Distance of the seed perigee upstream of its first hit, where the track finding starts
  double anchorDistance = 10. * Acts::UnitConstants::mm;

  // In EDM4eic units
  double positionError       = 0.1 * edm4eic::unit::mm;
  double angleError          = 2.e-3; // rad
  double qOverPRelativeError = 0.25;
  double timeError           = 10. * edm4eic::unit::ns;
};

/// Adds the seed of three hits on a helix in the field at the middle one, unless the helix has no
/// finite non-zero q/p or runs along the beam. Hit i of the collection is at positions[i], and
/// positions and field are in Acts units. Returns whether a seed was added.
bool addHelixSeed(const HelixSeedPrior& prior, const std::array<std::uint32_t, 3>& triplet,
                  float quality, const edm4eic::TrackerHitCollection& hits,
                  std::span<const Acts::Vector3> positions, const Acts::Vector3& field,
                  edm4eic::TrackSeedCollection& seeds,
                  edm4eic::TrackParametersCollection& parameters);

} // namespace eicrecon
