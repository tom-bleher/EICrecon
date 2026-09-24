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

  std::size_t numMeasurementsMin = 4;

  /// Express track parameters on the z = 0 plane instead of the perigee
  /// surface, for tracks nearly parallel to the beam line
  bool transverseReferencePlane = false;

  /// Apply material effects when extrapolating tracks to the reference surface
  bool referenceMaterialEffects = true;
};
} // namespace eicrecon
