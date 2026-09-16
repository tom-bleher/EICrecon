// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <stdexcept>

namespace eicrecon::b0 {
/// Fresh absolute initialization covariance, deliberately independent of the
/// hit-derived seed/posterior covariance. All variances are in EDM bound units.
inline Eigen::Matrix<double, 6, 6> weakPrior(const std::array<double, 6>& variances,
                                           double scale = 1.0) {
  if (!(scale > 0) || !std::isfinite(scale)) throw std::invalid_argument("Invalid weak-prior scale");
  Eigen::Matrix<double, 6, 6> result = Eigen::Matrix<double, 6, 6>::Zero();
  for (int i = 0; i < 6; ++i) {
    const double value = variances[i] * scale;
    if (!(value > 0) || !std::isfinite(value)) throw std::invalid_argument("Invalid weak-prior variance");
    result(i, i) = value;
  }
  return result;
}
} // namespace eicrecon::b0
