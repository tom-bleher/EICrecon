// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <vector>

namespace eicrecon {
struct CKFTrackingConfig {
  std::vector<double> etaBins                    = {};
  std::vector<double> chi2CutOff                 = {15.};
  std::vector<std::size_t> numMeasurementsCutOff = {10};

  std::size_t numMeasurementsMin = 4;

  /// For forward telescopes, extrapolate upstream from the first measurement.
  bool extrapolateBackwardFromFirst = false;

  /// Optional minimum number of longitudinal groups of accepted measurements.
  /// Zero disables the requirement, as appropriate for central tracking.
  std::size_t numMeasurementGroupsMin = 0;
  /// Gaps larger than this separate physical stations in a forward telescope.
  double measurementGroupZGap = 50. * Acts::UnitConstants::mm;
};
} // namespace eicrecon
