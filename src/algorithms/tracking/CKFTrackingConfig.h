// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <cstddef>
#include <string>
#include <vector>

namespace eicrecon {
struct CKFTrackingConfig {
  std::vector<double> etaBins                    = {};
  std::vector<double> chi2CutOff                 = {15.};
  std::vector<std::size_t> numMeasurementsCutOff = {10};

  std::size_t numMeasurementsMin = 4;

  /// Optional telescope coverage cut, disabled for ordinary CKF instances.
  std::size_t numStationsMin = 0;
  std::string stationReadout;
  double stationZGap = 50. * Acts::UnitConstants::mm;
};
} // namespace eicrecon
