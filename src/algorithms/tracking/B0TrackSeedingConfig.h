// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <Acts/Seeding2/BroadTripletSeedFilter.hpp>
#include <Acts/Seeding2/DoubletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeedFinder.hpp>
#include <edm4eic/unit_system.h>

namespace eicrecon {

/// Configuration of B0TrackSeeding.
struct B0TrackSeedingConfig {
  // The Acts seeder, in Acts units, with only what differs from its defaults. Every hit is tried as
  // a middle hit; in its frame, r' is about the distance along the beam and z' lies along the local
  // field.
  static constexpr float stationGap = 150. * Acts::UnitConstants::mm;
  static constexpr float noLimit    = 1.e5 * Acts::UnitConstants::mm;

  /// Doublets across stations. Tracks may come from decays upstream, so the origin is not cut on.
  Acts::DoubletSeedFinder::Config doublets{.deltaRMin          = stationGap,
                                           .deltaRMax          = noLimit,
                                           .collisionRegionMin = -noLimit,
                                           .collisionRegionMax = noLimit};
  Acts::TripletSeedFinder::Config triplets{.minPt     = 1. * Acts::UnitConstants::GeV,
                                           .impactMax = noLimit};
  Acts::BroadTripletSeedFilter::Config filter{.deltaRMin = stationGap, .impactWeightFactor = 0.F};

  /// Scale on the hit errors given to the Acts finders. Their slope cut has no floor on its
  /// scattering term, so at 100 GeV it cuts at about 1 sigma of the hit errors; this makes it at
  /// least hitErrorScale sigma.
  double hitErrorScale = 3.;

  /// Distance of the seed perigee upstream of its first hit, where the track finding starts
  double anchorDistance = 10. * Acts::UnitConstants::mm;

  // Seed prior, a loose diagonal envelope for the track finding, in EDM4eic units
  double positionError       = 0.1 * edm4eic::unit::mm;
  double angleError          = 2.e-3; // rad
  double qOverPRelativeError = 0.25;
  double timeError           = 10. * edm4eic::unit::ns;
};

} // namespace eicrecon
