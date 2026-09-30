// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck, Dmitry Romanov

#pragma once

namespace eicrecon {
struct TrackerHitReconstructionConfig {
  float timeResolution = 10;
  // Raw charge counts per keV, as set in the digitization
  int countsPerKeV = 1;
};
} // namespace eicrecon
