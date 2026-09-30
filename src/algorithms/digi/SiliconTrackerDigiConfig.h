// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten, Dmitry Romanov

#pragma once

#include <DD4hep/DD4hepUnits.h>

namespace eicrecon {

struct SiliconTrackerDigiConfig {
  // sub-systems should overwrite their own
  // NB: be aware of thresholds in npsim! E.g. https://github.com/eic/npsim/pull/9/files
  double threshold      = 0 * dd4hep::keV;
  double timeResolution = 8; /// TODO 8 of what units??? Same TODO in juggler. Probably [ns]
  // Apply the threshold to the energy summed over all SimTrackerHits in a cell, after adding
  // Gaussian noise of width hypot(noise, relativeNoise * energy), instead of to each
  // SimTrackerHit. Needed when charge sharing splits one deposit over several cells.
  bool thresholdOnCellSum = false;
  double noise            = 0 * dd4hep::keV;
  double relativeNoise    = 0;
  // Charge counts per keV of deposited energy in RawTrackerHit; TrackerHitReconstruction must use
  // the same value. 1 gives integer keV
  int countsPerKeV = 1;
};

} // namespace eicrecon
