// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <edm4eic/unit_system.h>
#include <cstddef>

namespace eicrecon {

/// Configuration of B0TrackSeeding, in Acts units unless noted
struct B0TrackSeedingConfig {
  /// Hits are seeded if they lie within this range of distance along the ion beam (r' in the
  /// seeding frame) and of distance from the beam along the field (z'); middle hits of a seed
  /// within the middle range
  double rMin       = 5700. * Acts::UnitConstants::mm;
  double rMax       = 7100. * Acts::UnitConstants::mm;
  double zMax       = 300. * Acts::UnitConstants::mm;
  double rMinMiddle = 6000. * Acts::UnitConstants::mm;
  double rMaxMiddle = 6800. * Acts::UnitConstants::mm;

  /// Doublets: distance along the beam between the middle and the outer hits, and the largest
  /// slope along the field
  double deltaRMin   = 150. * Acts::UnitConstants::mm;
  double deltaRMax   = 1100. * Acts::UnitConstants::mm;
  double cotThetaMax = 0.5;
  /// Minimum momentum across the field of the Acts finders
  double minPt = 0.5 * Acts::UnitConstants::GeV;
  /// Triplets: tolerance of the Acts scattering cut, and the material it assumes per station.
  /// With the frame along the field at each middle hit, 7 keeps 2-15 GeV pions and protons that
  /// the Acts default of 5 loses at low momentum.
  double sigmaScattering  = 7.;
  double radLengthPerSeed = 0.03;
  /// The finders take the field magnitude at construction, so they are built for bins of |B|
  /// of this width, in tesla
  double fieldBinWidth = 0.05;
  /// Maximum number of seeds per middle hit kept by the Acts seed filter
  std::size_t maxSeedsPerMiddle = 5;

  /// Hit-time window. A particle is never faster than light, so the later of two hits may lag
  /// behind a flight at the speed of light by at most that of the slowest seeded particle, of
  /// mass slowParticleMass at minMomentum, and may not be early, both within timeWindowSigmas
  /// times the hit time errors (at least timeErrorFloor each). 0 turns it off. For nuclear
  /// fragments set slowParticleMass to the mass per unit charge, as minMomentum is a cut on
  /// rigidity.
  double timeWindowSigmas = 5.;
  /// in EDM4eic units (ns)
  double timeErrorFloor   = 0.01 * edm4eic::unit::ns;
  double slowParticleMass = 0.938272 * Acts::UnitConstants::GeV;

  /// Minimum momentum of a seed
  double minMomentum = 1. * Acts::UnitConstants::GeV;
  /// Seeds are expressed on a perigee this far upstream of their first hit, where the track
  /// finding starts
  double anchorDistance = 10. * Acts::UnitConstants::mm;
  /// Maximum number of seeds per event; the seeds with the highest Acts weight are kept
  std::size_t maxSeeds = 1000;

  /// Diagonal prior of the seed parameters, a loose envelope that the track finding refits
  double positionError       = 0.1 * Acts::UnitConstants::mm;
  double angleError          = 2. * Acts::UnitConstants::mrad;
  double qOverPRelativeError = 0.25;
  /// in EDM4eic units (ns)
  double timeError = 10. * edm4eic::unit::ns;
};

} // namespace eicrecon
