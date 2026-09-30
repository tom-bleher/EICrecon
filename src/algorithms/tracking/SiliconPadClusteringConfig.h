// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/unit_system.h>
#include <string>

namespace eicrecon {

struct SiliconPadClusteringConfig {
  std::string readout;
  // Pads of one sensor that touch, including diagonally, and are closer than deltaT in time
  // form one cluster
  double deltaT = 1 * edm4eic::unit::ns;
  // Discriminator threshold: pads below it are read out only next to (including diagonally) a
  // pad above it, as for a readout of fired pads and their neighbours; 0 reads all pads
  double seed_threshold = 0;
  // Metal electrode size (DD4hep length units) used to invert the trapezoid pad response of
  // TrapezoidPadResponse.h; 0 gives the two-pad centroid
  double electrode_x = 0;
  double electrode_y = 0;
  // Position resolution (DD4hep length units) along an axis where the cluster is one pad wide,
  // or two pads wide; 0 means pitch/sqrt(12), which is also used for wider clusters
  double single_pad_resolution = 0;
  double shared_pad_resolution = 0;
  // Clusters wider than this many pads along an axis are irregular and get pitch/sqrt(12); two
  // for a trapezoid response, more when the response leaks to further pads
  int max_regular_width = 2;
  // Optional calibrated two-pad position curve (SharedOffsetCurve.h) replacing the trapezoid
  // inversion, for sensors whose sharing is not trapezoidal
  std::string inversion_file;
};

} // namespace eicrecon
