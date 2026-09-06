// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <vector>

namespace eicrecon {
struct CKFTrackingConfig {
  std::vector<double> etaBins                    = {};
  std::vector<double> chi2CutOff                 = {15.};
  std::vector<std::size_t> numMeasurementsCutOff = {10};

  // Absolute PDG for the fit mass/charge hypothesis, not reconstructed PID.
  // Supported charged species: 11, 13, 211, 321, 2212. Charge sign comes from q/p.
  int particleHypothesisPdg = 211;

  std::size_t numMeasurementsMin = 4;
  // Disabled for central tracking; B0 requires independent physical stations.
  std::size_t numB0StationsMin = 0;
  double b0StationZGap         = 50.0; // mm in the ion frame, matching B0 stub seeding
};
} // namespace eicrecon
