// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <edm4eic/unit_system.h>

namespace eicrecon {

/// Configuration of B0TrackSeeding, in Acts units unless noted.
/// In the seeding frame of a middle hit, r' is about the distance along the ion beam and z' the
/// distance from the beam along the local field. No station positions are configured: every hit is
/// tried as a middle hit, and the doublet distances below select the stations.
struct B0TrackSeedingConfig {
  /// Minimum distance in r' between the middle and an outer hit
  double deltaRMin = 150. * Acts::UnitConstants::mm;
  /// Maximum distance in r' between the middle and an outer hit
  double deltaRMax = 1100. * Acts::UnitConstants::mm;
  /// Maximum doublet slope along the field
  double cotThetaMax = 0.5;
  /// Minimum momentum across the field in the Acts finders
  double minPt = 0.5 * Acts::UnitConstants::GeV;
  /// Number of sigmas of the Acts scattering cut. With radLengthPerSeed 0.03, 7 keeps the 2-15 GeV
  /// pions and protons that 5 loses.
  double sigmaScattering = 7.;
  /// Radiation lengths per station assumed by the scattering cut
  double radLengthPerSeed = 0.03;
  /// Scale on the hit errors given to the Acts finders. Their slope cut has no floor on its
  /// scattering term, so at 100 GeV it cuts at about 1 sigma of the hit errors; this makes it at
  /// least hitErrorScale sigma.
  double hitErrorScale = 3.;

  /// Minimum seed momentum
  double minMomentum = 1. * Acts::UnitConstants::GeV;
  /// Distance of the seed perigee upstream of its first hit, where the track finding starts
  double anchorDistance = 10. * Acts::UnitConstants::mm;

  /// Seed prior: a loose diagonal envelope for the track finding
  double positionError = 0.1 * Acts::UnitConstants::mm;
  double angleError    = 2. * Acts::UnitConstants::mrad;
  /// Relative q/p error
  double qOverPRelativeError = 0.25;
  /// Seed time error, in EDM4eic units (ns)
  double timeError = 10. * edm4eic::unit::ns;
};

} // namespace eicrecon
