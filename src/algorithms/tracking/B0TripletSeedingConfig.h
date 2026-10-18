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
  /// Maximum difference between the time difference of the first hit and each of the other two and
  /// the flight time between them at the speed of light, in EDM4eic units (ns). It rejects hits
  /// from other bunch crossings before the residual is computed. A proton of 1 GeV (beta = 0.73)
  /// lags light by 1 ns across the four stations. 0 (the default) turns the check off.
  double maxTimeResidual = 0.;
  /// Minimum momentum of the three-hit estimate
  double minMomentum = 1. * Acts::UnitConstants::GeV;
  /// Distance upstream of the first hit at which the seed is expressed
  double anchorDistance = 10. * Acts::UnitConstants::mm;
  /// Maximum number of seeds per event. Candidates are ranked by residual alone and the rest are
  /// dropped, so in a busy event genuine seeds can be lost. A warning is logged the 1st, 2nd, 4th,
  /// 8th, ... time it happens.
  std::size_t maxSeeds = 1000;

  /// Skip a candidate whose three hits all belong to better-ranked seeds. A track on four stations
  /// then gives two seeds instead of four. Off by default: the track finding does not succeed from
  /// every seed of a track, and keeping only one seed per track cost efficiency (dd47ba24).
  bool dropSeedsWithUsedHits = false;

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
