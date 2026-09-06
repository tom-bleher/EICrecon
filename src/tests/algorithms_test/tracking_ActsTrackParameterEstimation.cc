// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Acts/Definitions/Units.hpp>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "algorithms/tracking/ActsTrackParameterEstimation.h"

using eicrecon::acts_compat::estimateTrackParamsFromSpacePoints;

TEST_CASE("Multipoint parameters reproduce helices in arbitrary fields", "[tracking][multipoint]") {
  for (double charge : {-1., 1.}) {
    for (double angle : {0., .7}) {
      for (double pitch : {-.3, 0., .4}) {
        const double radius = 20. / (charge * .000299792458 * 1.2);
        const auto rotation =
            Eigen::AngleAxisd(angle, Acts::Vector3(1., 2., 3.).normalized()).toRotationMatrix();
        const auto point = [&](double length) -> Acts::Vector3 {
          return rotation * Acts::Vector3(radius * (std::cos(length / radius) - 1.), pitch * length,
                                          radius * std::sin(length / radius));
        };
        const std::vector<Acts::Vector3> points{point(0.), point(300.), point(600.), point(1000.)};
        const Acts::Vector3 field = rotation * Acts::Vector3(0., 1.2 * Acts::UnitConstants::T, 0.);
        const Acts::Vector3 direction = rotation * Acts::Vector3(0., pitch, 1.).normalized();
        for (std::size_t refine : {0U, 3U}) {
          const auto result = estimateTrackParamsFromSpacePoints(points, field, 2., refine);
          REQUIRE(result.has_value());
          CHECK((*result)[Acts::eFreeQOverP] ==
                Catch::Approx(charge / (20. * std::sqrt(1. + pitch * pitch))).margin(1e-8));
          CHECK((result->segment<3>(Acts::eFreeDir0) - direction).norm() < 1e-8);
          CHECK((result->segment<3>(Acts::eFreePos0) - points[0]).norm() < 1e-12);
          CHECK((*result)[Acts::eFreeTime] == 2.);
        }
        // A zero-weight fourth point is a valid way to fit only the first three.
        const std::array<double, 4> weights{1., 1., 1., 0.};
        const auto weighted = estimateTrackParamsFromSpacePoints(points, field, 0., 0, weights);
        REQUIRE(weighted.has_value());
        CHECK((*weighted)[Acts::eFreeQOverP] ==
              Catch::Approx(charge / (20. * std::sqrt(1. + pitch * pitch))).margin(1e-8));
      }
    }
  }
}

TEST_CASE("Multipoint parameters reject invalid input and handle zero field",
          "[tracking][multipoint]") {
  const std::vector<Acts::Vector3> line{{0., 0., 0.}, {1., 2., 3.}, {2., 4., 6.}};
  const auto result = estimateTrackParamsFromSpacePoints(line, Acts::Vector3::Zero());
  REQUIRE(result.has_value());
  CHECK((*result)[Acts::eFreeQOverP] == 0.);
  CHECK((result->segment<3>(Acts::eFreeDir0) - Acts::Vector3(1., 2., 3.).normalized()).norm() <
        1e-12);
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(std::span(line).first(2), Acts::Vector3::UnitY()));
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(line, Acts::Vector3::UnitY(), 0., 0, {}, 3));
  const std::vector<Acts::Vector3> coincident(3, Acts::Vector3::Zero());
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(coincident, Acts::Vector3::UnitY()));
  const std::array<double, 2> mismatched{1., 1.};
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(line, Acts::Vector3::UnitY(), 0., 0, mismatched));
  const std::array<double, 3> negative{-1., 1., 1.};
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(line, Acts::Vector3::UnitY(), 0., 0, negative));
  const std::array<double, 3> zero{0., 0., 0.};
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(line, Acts::Vector3::UnitY(), 0., 0, zero));
  auto nonfinite   = line;
  nonfinite[0].x() = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(nonfinite, Acts::Vector3::UnitY()));
  CHECK_FALSE(estimateTrackParamsFromSpacePoints(line, nonfinite[0]));
}
