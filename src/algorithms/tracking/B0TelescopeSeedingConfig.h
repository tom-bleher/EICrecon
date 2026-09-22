// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#pragma once

#include <edm4eic/unit_system.h>
#include <cstddef>

namespace eicrecon {

struct B0TelescopeSeedingConfig {
  // Explicit mapping of readout layer IDs: 2 for realistic front/back, 1 for ideal disks.
  unsigned layersPerStation = 2;
  double minDeltaZ = 50.0 * edm4eic::unit::mm;
  double maxSlope = 0.15;
  // Loose commissioning windows, not tuned resolutions or background rejection.
  double maxResidualX = 20.0 * edm4eic::unit::mm;
  double maxResidualY = 5.0 * edm4eic::unit::mm;
  std::size_t maxCombinations = 100000;
  std::size_t maxSeedsPerMiddle = 8;
  std::size_t maxSeeds = 64;

  // Diagonal initial uncertainties at the first-hit plane, propagated by ACTS.
  double positionSigma = 0.1 * edm4eic::unit::mm;
  double angleSigma = 0.01; // radians
  double relativeMomentumSigma = 0.5;
  double qOverPSigmaFloor = 0.002 / edm4eic::unit::GeV;
  double timeSigma = 10.0 * edm4eic::unit::ns;
  double propagationPathLimit = 15000.0 * edm4eic::unit::mm;
};

} // namespace eicrecon
