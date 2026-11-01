// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <cstddef>
#include <limits>
#include <string>

#include "algorithms/tracking/B0TripletSeedingConfig.h"

namespace eicrecon {

/// Configuration of B0Seeding2: the Acts Seeding2 triplet seeder run on B0 hits in a frame whose
/// z axis is along the magnet field. Lengths in Acts units (mm), momenta in GeV, fields in tesla.
/// All values are finite (JANA turns infinity into 0). The cuts that assume a beam line
/// (collision region, impact, impact and z-origin weights) are opened up.
struct B0Seeding2Config {
  /// Seed writing (anchor, prior), minMomentum, the hit-time window, the effective field and
  /// maxSeeds are shared with B0TripletSeeding so that the two write identical seeds for the same
  /// hit triplet
  B0TripletSeedingConfig triplet;

  /// How the field enters:
  ///  "perMiddle" frame z' along the field at each middle hit, finders built for |B| there (in
  ///              bins of bFieldBinWidth), seed parameters in the effective field along the
  ///              triplet (B0TripletSeeding fieldSamples). The default.
  ///  "effective" frame z' along global y, finders in bFieldInZ, seed parameters as "perMiddle".
  ///              The quadrupole turns the field off the midplane, so the Acts r-z' cut then
  ///              needs sigmaScattering of about 150 not to lose tracks there.
  ///  "scalar"    as "effective", but seed parameters in bFieldInZ as well, which is off by up to
  ///              40 % in q/p across the acceptance (for comparison only)
  std::string fieldModel = "perMiddle";
  /// Field magnitude of the Acts finders in the fixed frame, in tesla
  double bFieldInZ = 1.58;
  /// Bin width of |B| for the per-middle finders, in tesla
  double bFieldBinWidth = 0.05;
  /// Angle of the ion beam about global y; NaN reads the DD4hep constant CrossingAngle
  double frameAngle = std::numeric_limits<double>::quiet_NaN();

  /// Space point window in the primed frame: r' is the distance along the ion beam, z' the
  /// coordinate along the field
  double rMin       = 5700. * Acts::UnitConstants::mm;
  double rMax       = 7100. * Acts::UnitConstants::mm;
  double zMin       = -300. * Acts::UnitConstants::mm;
  double zMax       = 300. * Acts::UnitConstants::mm;
  double rMinMiddle = 6000. * Acts::UnitConstants::mm;
  double rMaxMiddle = 6800. * Acts::UnitConstants::mm;
  /// Added in quadrature to the hit error along z' (a field-model error for the r-z' cuts)
  double extraSigmaZ = 0. * Acts::UnitConstants::mm;

  /// DoubletSeedFinder
  double deltaRMinBottom    = 150. * Acts::UnitConstants::mm;
  double deltaRMaxBottom    = 1100. * Acts::UnitConstants::mm;
  double deltaRMinTop       = 150. * Acts::UnitConstants::mm;
  double deltaRMaxTop       = 1100. * Acts::UnitConstants::mm;
  double deltaZMin          = -1000. * Acts::UnitConstants::mm;
  double deltaZMax          = 1000. * Acts::UnitConstants::mm;
  double collisionRegionMin = -10000. * Acts::UnitConstants::mm;
  double collisionRegionMax = 10000. * Acts::UnitConstants::mm;
  double cotThetaMax        = 0.5;
  /// DoubletSeedFinder and TripletSeedFinder
  double impactMax         = 1.e5 * Acts::UnitConstants::mm;
  double minPt             = 0.5 * Acts::UnitConstants::GeV;
  double helixCutTolerance = 1.;
  /// TripletSeedFinder. In the per-middle frame 7 keeps the efficiency of B0TripletSeeding on
  /// DVCS protons, 2-15 GeV pions and protons and Lambda daughters (the Acts default 5 loses up
  /// to 2 % of low-momentum tracks).
  double sigmaScattering  = 7.;
  double radLengthPerSeed = 0.03;
  bool sortedByCotTheta   = true;

  /// BroadTripletSeedFilter
  double deltaInvHelixDiameter          = 0.00003 * (1. / Acts::UnitConstants::mm);
  double filterDeltaRMin                = 150. * Acts::UnitConstants::mm;
  double compatSeedWeight               = 200.;
  std::size_t compatSeedLimit           = 2;
  double impactWeightFactor             = 0.;
  double zOriginWeightFactor            = 0.;
  std::size_t maxSeedsPerSpM            = 5;
  std::size_t maxSeedsPerSpMConf        = 10;
  std::size_t maxQualitySeedsPerSpMConf = 10;
};

} // namespace eicrecon
