// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten, Dmitry Romanov

#pragma once

#include <DD4hep/DD4hepUnits.h>

namespace eicrecon {

struct SiliconTrackerDigiConfig {
  // Sub-systems should overwrite their own values.
  // NB: be aware of thresholds in npsim! E.g. https://github.com/eic/npsim/pull/9/files
  double threshold = 0 * dd4hep::keV;

  // Gaussian sigma [ns] of the generic *channel-level* time measurement.
  // Contributions in one cell are first aggregated. The channel's deterministic
  // time is the earliest unsmeared SimTrackerHit time; one Gaussian error is then
  // applied to that channel, independent of contribution multiplicity. The
  // resulting RawTrackerHit timestamp is stored in integer picoseconds.
  // Detector-specific discriminator/time-walk/common-clock models should use a
  // dedicated response implementation rather than changing these semantics.
  double timeResolution = 8.0;
};

} // namespace eicrecon
