// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <edm4eic/unit_system.h>
#include <cstddef>

namespace eicrecon {

struct B0TripletSeedingConfig {
  /// A z gap between consecutive hits larger than this starts a new station. It
  /// must exceed the z spread within a station (up to 12 mm) and stay below the
  /// station spacing (at least 270 mm).
  double stationGap = 50. * Acts::UnitConstants::mm;
  /// Maximum distance of the middle hit from the chord of the outer two,
  /// measured along the local magnetic field
  double maxResidual = 1. * Acts::UnitConstants::mm;
  /// Minimum momentum of the three-hit estimate
  double minMomentum = 1. * Acts::UnitConstants::GeV;
  /// Distance upstream of the first hit at which the seed is expressed. The seed is moved there
  /// along a straight line, and its position error grows by the sagitta of the track over it.
  double anchorDistance = 10. * Acts::UnitConstants::mm;
  /// Maximum number of seeds per event. Candidates are ranked by residual alone and the rest are
  /// dropped, so in a busy event genuine seeds can be lost. A warning is logged the 1st, 2nd, 4th,
  /// 8th, ... time it happens.
  std::size_t maxSeeds = 1000;

  /// Hit-time window. A particle is never faster than light, so the later of two hits may only
  /// lag by up to the flight time of the slowest seeded particle, distance * (1 / beta_min - 1),
  /// and may not be early. Both bounds are widened by timeWindowSigmas times the combined hit time
  /// error (TrackerHit::getTimeError(), at least timeErrorFloor each) and by timeWindowTolerance.
  /// beta_min is that of a particle of mass slowParticleMass at minMomentum. With hit time errors
  /// of several ns (before the B0 timing was set to 30 ps) the window rejects nothing.
  /// 0 turns the window off.
  double timeWindowSigmas = 5.;
  /// in EDM4eic units (ns)
  double timeWindowTolerance = 0.;
  /// in EDM4eic units (ns)
  double timeErrorFloor = 0.01;
  /// Mass per unit charge of the slowest particles to seed. minMomentum is a cut on rigidity, so
  /// for nuclear fragments set the mass per unit charge: about 1.876 GeV for deuterons, alpha
  /// particles and carbon nuclei.
  double slowParticleMass = 0.938272 * Acts::UnitConstants::GeV;

  /// Prior uncertainties of the seed parameters. The prior is a deliberately loose diagonal
  /// envelope, not an error estimate: the track finder refits from it, so it only has to cover
  /// the true error. The relative q/p error is a floor that dominates for most seeds. Widening
  /// the angle error to 10 or 30 mrad did not change the efficiency for 2-15 GeV pions. Seeds
  /// built from other particles or from stray hits can have the wrong charge.
  double positionError       = 0.1 * Acts::UnitConstants::mm;
  double angleError          = 2. * Acts::UnitConstants::mrad;
  double qOverPRelativeError = 0.25;
  /// in EDM4eic units (ns), unlike the other values: Acts measures time as a length
  double timeError = 10. * edm4eic::unit::ns;
};

} // namespace eicrecon
