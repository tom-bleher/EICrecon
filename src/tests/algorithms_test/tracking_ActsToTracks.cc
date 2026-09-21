// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <catch2/catch_test_macros.hpp>
#include <array>
#include <cmath>

#include "algorithms/tracking/BoundToCartesianCovariance.h"

TEST_CASE("Cartesian track covariance agrees with independent finite differences",
          "[tracking][actstotracks]") {
  using Matrix = Eigen::Matrix<double, 6, 6>;
  using Vector = Eigen::Matrix<double, 6, 1>;
  using namespace Acts::UnitConstants;
#if Acts_VERSION_MAJOR >= 45
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  const Acts::GeometryContext gctx{};
#endif

  // A correlated covariance with different physical scales, including time.
  Matrix factor = Matrix::Identity();
  for (int row = 0; row < 6; ++row) {
    for (int col = 0; col < row; ++col) {
      factor(row, col) = 0.07 * (row + col + 1);
    }
  }
  Vector scales;
  scales << 0.2 * mm, 4.0 * mm, 0.003, 0.0002, 0.001 / GeV, 2.0 * ns;
  factor                  = scales.asDiagonal() * factor;
  const Matrix covariance = factor * factor.transpose();
  const std::array<double, 6> steps{1e-4 * mm, 1e-4 * mm, 1e-6, 1e-6, 1e-7 / GeV, 1e-4 * ns};

  for (const double tilt : {0.0, 0.025}) {
    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.linear()         = Eigen::AngleAxisd(tilt, Acts::Vector3::UnitY()).toRotationMatrix();
    transform.translation()    = Acts::Vector3(7.0 * mm, -11.0 * mm, 37.0 * mm);
    const auto surface         = Acts::Surface::makeShared<Acts::PerigeeSurface>(transform);
    for (const double charge : {-1.0, 1.0}) {
      for (const double theta : {0.03, 1.2, 2.8}) {
        for (const double impact : {-3.0 * mm, 0.0 * mm, 2.0 * mm}) {
          Acts::BoundVector parameters;
          parameters << impact, 43.0 * mm, 0.7, theta, charge / (41.0 * GeV), 5.0 * ns;
          CAPTURE(tilt, charge, theta, impact);

          // Differentiate the actual position/momentum map, not its Jacobian.
          const auto cartesian = [&](const Acts::BoundVector& pars) -> Vector {
            const double phi   = pars[Acts::eBoundPhi];
            const double polar = pars[Acts::eBoundTheta];
            const Acts::Vector3 direction(std::cos(phi) * std::sin(polar),
                                          std::sin(phi) * std::sin(polar), std::cos(polar));
            Vector result;
            result.head<3>() = surface->localToGlobal(gctx, pars.head<2>(), direction) / mm;
            result.tail<3>() = std::abs(1.0 / pars[Acts::eBoundQOverP]) * direction / GeV;
            return result;
          };
          Matrix numericalJacobian;
          for (int col = 0; col < 6; ++col) {
            auto plus  = parameters;
            auto minus = parameters;
            plus[col] += steps[col];
            minus[col] -= steps[col];
            numericalJacobian.col(col) = (cartesian(plus) - cartesian(minus)) / (2.0 * steps[col]);
          }
          const Matrix expected = numericalJacobian * covariance * numericalJacobian.transpose();
          const Matrix actual =
              eicrecon::boundToCartesianCovariance(*surface, gctx, parameters, covariance);
          REQUIRE(actual.allFinite());
          CHECK((actual.block<3, 3>(0, 3).norm() > 0.0));
          for (int row = 0; row < 6; ++row) {
            for (int col = row; col < 6; ++col) {
              const double scale = std::sqrt(expected(row, row) * expected(col, col));
              CHECK(std::abs(actual(row, col) - expected(row, col)) < 2e-6 * scale);
              CHECK(std::abs(actual(row, col) - actual(col, row)) < 1e-12 * scale);
            }
          }
        }
      }
    }
  }
}
