// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Eigen/Cholesky>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include "algorithms/tracking/B0TrackerActsSeeding.h"

namespace {
const auto geometryContext = Acts::GeometryContext::dangerouslyDefaultConstruct();

struct HelixMeasurements {
  std::vector<std::shared_ptr<Acts::PlaneSurface>> surfaces;
  std::vector<eicrecon::b0acts::Measurement> measurements;
  Acts::Vector3 field;
  Acts::Vector3 direction;
  Acts::Vector3 first;
  double qOverP{};
};

HelixMeasurements makeMeasurements(std::size_t count, double charge, bool phiBoundary = false) {
  const double pitch = phiBoundary ? 0. : .13;
  const Acts::RotationMatrix3 rotation =
      phiBoundary
          ? Eigen::AngleAxisd(-.035, Acts::Vector3::UnitY()).toRotationMatrix()
          : Eigen::AngleAxisd(.6, Acts::Vector3(1., 2., 3.).normalized()).toRotationMatrix();
  const double radius = 20. / (charge * .000299792458 * 1.2);
  HelixMeasurements result;
  result.first     = {-150., 30., 6100.};
  result.field     = rotation * Acts::Vector3(0., 1.2 * Acts::UnitConstants::T, 0.);
  result.direction = rotation * Acts::Vector3(0., pitch, 1.).normalized();
  result.qOverP    = charge / (20. * std::sqrt(1. + pitch * pitch));
  for (std::size_t hit = 0; hit < count; ++hit) {
    const double length = 750. * hit / (count - 1);
    const Acts::Vector3 global =
        result.first + rotation * Acts::Vector3(radius * (std::cos(length / radius) - 1.),
                                                pitch * length, radius * std::sin(length / radius));
    const Acts::RotationMatrix3 sensorRotation =
        (Eigen::AngleAxisd(.11 + .02 * hit, Acts::Vector3::UnitX()) *
         Eigen::AngleAxisd(-.18 + .01 * hit, Acts::Vector3::UnitY()) *
         Eigen::AngleAxisd(.27 * hit, Acts::Vector3::UnitZ()))
            .toRotationMatrix();
    const Acts::Vector2 local(1.4 + .2 * hit, -2.1 + .3 * hit);
    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.linear()         = sensorRotation;
    transform.translation()    = global - sensorRotation * Acts::Vector3(local.x(), local.y(), 0.);
    result.surfaces.push_back(Acts::Surface::makeShared<Acts::PlaneSurface>(transform));
    // Unequal sensor errors and nonzero local xy correlation; units are mm^2.
    const double sigma0      = .004 * (1. + .3 * hit);
    const double sigma1      = .009 * (1. + .2 * hit);
    const double correlation = hit % 2 == 0 ? .35 : -.25;
    Acts::SquareMatrix2 covariance;
    covariance << sigma0 * sigma0, correlation * sigma0 * sigma1, correlation * sigma0 * sigma1,
        sigma1 * sigma1;
    result.measurements.push_back({result.surfaces.back().get(), local, covariance});
  }
  return result;
}
} // namespace

TEST_CASE("B0 multipoint estimates bind to tilted sensor measurements", "[tracking][b0acts]") {
  for (const auto count : {3U, 4U}) {
    for (const auto charge : {-1., 1.}) {
      const auto sample   = makeMeasurements(count, charge);
      const auto estimate = eicrecon::b0acts::estimateWithCovariance(sample.measurements,
                                                                     geometryContext, sample.field);
      REQUIRE(estimate.has_value());
      const auto& parameters = estimate->parameters;
      CHECK((parameters.head<2>() - sample.measurements.front().local).norm() < 1e-10);
      CHECK(parameters(4) == Catch::Approx(sample.qOverP).margin(1e-8));
      const Acts::Vector3 direction(std::sin(parameters(3)) * std::cos(parameters(2)),
                                    std::sin(parameters(3)) * std::sin(parameters(2)),
                                    std::cos(parameters(3)));
      CHECK((direction - sample.direction).norm() < 1e-8);
      const auto global =
          sample.surfaces.front()->localToGlobal(geometryContext, parameters.head<2>(), direction);
      CHECK((global - sample.first).norm() < 1e-9);
      const auto local = sample.surfaces.front()->globalToLocal(geometryContext, global, direction);
      REQUIRE(local.ok());
      CHECK((*local - sample.measurements.front().local).norm() < 1e-10);
      CHECK(estimate->covariance.allFinite());
      CHECK(estimate->covariance.isApprox(estimate->covariance.transpose()));
      CHECK(Eigen::LLT<eicrecon::b0acts::Covariance>(estimate->covariance).info() ==
            Eigen::Success);
    }
  }
}

TEST_CASE("B0 seed covariance predicts independent local-coordinate noise", "[tracking][b0acts]") {
  constexpr std::size_t throws = 3000;
  for (const auto count : {3U, 4U}) {
    // The nominal direction lies at the phi branch cut. Wrapping is essential.
    const auto sample     = makeMeasurements(count, 1., true);
    const auto prediction = eicrecon::b0acts::estimateWithCovariance(sample.measurements,
                                                                     geometryContext, sample.field);
    REQUIRE(prediction.has_value());
    Eigen::LLT<eicrecon::b0acts::Covariance> factor(prediction->covariance);
    REQUIRE(factor.info() == Eigen::Success);
    std::vector<Acts::SquareMatrix2> noiseFactors;
    for (const auto& measurement : sample.measurements) {
      noiseFactors.emplace_back(Eigen::LLT<Acts::SquareMatrix2>(measurement.covariance).matrixL());
    }
    std::mt19937_64 random(98131 + count);
    std::normal_distribution<double> normal;
    eicrecon::b0acts::Parameters sum    = eicrecon::b0acts::Parameters::Zero();
    eicrecon::b0acts::Covariance square = eicrecon::b0acts::Covariance::Zero();
    for (std::size_t trial = 0; trial < throws; ++trial) {
      auto noisy = sample.measurements;
      for (std::size_t hit = 0; hit < noisy.size(); ++hit) {
        const Acts::Vector2 draw(normal(random), normal(random));
        noisy[hit].local += noiseFactors[hit] * draw;
      }
      // This validation never calls the covariance/Jacobian implementation on
      // noisy samples: only the independently perturbed parameter estimates.
      const auto observed =
          eicrecon::b0acts::estimateParameters(noisy, geometryContext, sample.field);
      REQUIRE(observed.has_value());
      eicrecon::b0acts::Parameters residual = *observed - prediction->parameters;
      residual(2)                           = std::remainder(residual(2), 2. * std::acos(-1.));
      const eicrecon::b0acts::Parameters whitened = factor.matrixL().solve(residual);
      sum += whitened;
      square += whitened * whitened.transpose();
    }
    const auto mean       = (sum / throws).eval();
    const auto covariance = ((square - throws * mean * mean.transpose()) / (throws - 1.)).eval();
    INFO("Number of hits: " << count);
    INFO("Whitened empirical covariance:\n" << covariance);
    CHECK(mean.cwiseAbs().maxCoeff() < .10);
    // 3,000 independent throws give ~2.6% diagonal statistical uncertainty.
    // Check the entire matrix, so incorrect xy terms and missing correlations
    // cannot pass just because the marginal variances happen to agree.
    CHECK((covariance - eicrecon::b0acts::Covariance::Identity()).cwiseAbs().maxCoeff() < .10);
  }
}

TEST_CASE("B0 multipoint estimates reject malformed geometry and covariance",
          "[tracking][b0acts]") {
  const auto sample   = makeMeasurements(4, 1.);
  auto bad            = sample.measurements;
  bad.front().surface = nullptr;
  CHECK_FALSE(eicrecon::b0acts::estimateParameters(bad, geometryContext, sample.field));
  CHECK_FALSE(eicrecon::b0acts::estimateWithCovariance(bad, geometryContext, sample.field));
  bad                   = sample.measurements;
  bad.front().local.x() = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(eicrecon::b0acts::estimateParameters(bad, geometryContext, sample.field));
  Acts::Transform3 invalidTransform  = Acts::Transform3::Identity();
  invalidTransform.translation().x() = std::numeric_limits<double>::quiet_NaN();
  const auto invalidSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(invalidTransform);
  bad                       = sample.measurements;
  bad.front().surface       = invalidSurface.get();
  CHECK_FALSE(eicrecon::b0acts::estimateParameters(bad, geometryContext, sample.field));
  CHECK_FALSE(eicrecon::b0acts::estimateParameters(sample.measurements, geometryContext,
                                                   Acts::Vector3::Zero()));
  bad = {sample.measurements.front(), sample.measurements.front(), sample.measurements.front()};
  CHECK_FALSE(eicrecon::b0acts::estimateParameters(bad, geometryContext, sample.field));
  bad = sample.measurements;
  bad.resize(2);
  CHECK_FALSE(eicrecon::b0acts::estimateParameters(bad, geometryContext, sample.field));
  for (const int invalid : {0, 1, 2}) {
    bad = sample.measurements;
    if (invalid == 0) {
      bad[1].covariance(0, 0) = 0.;
    }
    if (invalid == 1) {
      bad[1].covariance(0, 1) = 1.;
    }
    if (invalid == 2) {
      bad[1].covariance(0, 0) = -1.;
    }
    CHECK_FALSE(eicrecon::b0acts::estimateWithCovariance(bad, geometryContext, sample.field));
  }
}
