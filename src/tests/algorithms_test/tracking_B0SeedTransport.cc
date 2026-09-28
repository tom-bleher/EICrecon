// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Eigen/Eigenvalues>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <optional>

#include "algorithms/tracking/B0SeedTransport.h"

namespace {
Acts::GeometryContext geometryContext() {
#if Acts_VERSION_MAJOR >= 45
  return Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  return Acts::GeometryContext{};
#endif
}

Acts::BoundTrackParameters makeSeed(double theta, double qOverP) {
  const auto surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{-150., 20., 5800.});
  Acts::BoundVector parameters;
  parameters << 0., 0., 0., theta, qOverP, 20. * Acts::UnitConstants::ns;
  Acts::BoundVector errors;
  errors << 0.1, 0.1 / std::tan(theta), 0.002 / std::sin(theta), 0.002, 0.25 * std::abs(qOverP),
      10. * Acts::UnitConstants::ns;
  const Acts::BoundTrackParameters::CovarianceMatrix covariance = errors.cwiseAbs2().asDiagonal();
  return {surface, parameters, covariance, Acts::ParticleHypothesis::pion()};
}
} // namespace

TEST_CASE("B0 local transport follows a helix for either charge", "[B0SeedTransport]") {
  const auto gctx = geometryContext();
  const Acts::MagneticFieldContext mctx;
  const double by  = Acts::UnitConstants::T;
  const auto field = std::make_shared<Acts::ConstantBField>(Acts::Vector3{0., by, 0.});
  for (double theta : {0.001, 0.020, 0.050}) {
    for (double qOverP : {-1., -0.1, 0.1, 1.}) {
      for (double distance : {5., 10., 20.}) {
        CAPTURE(theta, qOverP, distance);
        const auto start  = makeSeed(theta, qOverP);
        const auto result = eicrecon::transportB0Seed(start, distance, field, gctx, mctx);
        REQUIRE(result);
        const double curvature = qOverP * by;
        const double angle     = theta + curvature * distance;
        const Acts::Vector3 expectedDirection{std::sin(angle), 0., std::cos(angle)};
        const Acts::Vector3 expectedPosition =
            start.position(gctx) + Acts::Vector3{(std::cos(angle) - std::cos(theta)) / curvature,
                                                 0.,
                                                 (std::sin(theta) - std::sin(angle)) / curvature};
        CHECK((result->direction() - expectedDirection).norm() < 1.e-9);
        CHECK((result->position(gctx) - expectedPosition).norm() < 1.e-6);
        const double mass         = start.particleHypothesis().mass();
        const double expectedTime = start.time() - distance * std::hypot(1., mass * qOverP);
        CHECK(result->time() == Catch::Approx(expectedTime).margin(1.e-8));
        CHECK(result->parameters()[Acts::eBoundQOverP] == qOverP);
        CHECK(result->parameters().head<2>().norm() < 1.e-9);
        REQUIRE(result->covariance());
        CHECK(result->covariance()->allFinite());
        const Eigen::SelfAdjointEigenSolver<Acts::BoundTrackParameters::CovarianceMatrix> eig(
            *result->covariance());
        REQUIRE(eig.info() == Eigen::Success);
        CHECK(eig.eigenvalues().minCoeff() > 0.);
      }
    }
  }
}

TEST_CASE("B0 local transport retains the covariance under round trip", "[B0SeedTransport]") {
  const auto gctx = geometryContext();
  const Acts::MagneticFieldContext mctx;
  const auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3{0., Acts::UnitConstants::T, 0.});
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  Propagator propagator{Acts::EigenStepper<>{field}};
  Propagator::Options<> options(gctx, mctx);
  options.pathLimit = 10.;
  for (double theta : {0.020, 0.050}) {
    for (double qOverP : {-1., 1.}) {
      CAPTURE(theta, qOverP);
      const auto start    = makeSeed(theta, qOverP);
      const auto upstream = eicrecon::transportB0Seed(start, 10., field, gctx, mctx);
      REQUIRE(upstream);
      auto result = propagator.propagate(*upstream, options);
      REQUIRE(result.ok());
      REQUIRE(result->endParameters);
      Acts::EigenStepper<> stepper{field};
      auto state = stepper.makeState(options.stepping);
      stepper.initialize(state, *result->endParameters);
      REQUIRE(stepper.prepareCurvilinearState(state));
      const auto rebound = stepper.boundState(state, start.referenceSurface());
      REQUIRE(rebound.ok());
      const auto& end = std::get<0>(*rebound);
      CHECK((end.position(gctx) - start.position(gctx)).norm() < 1.e-5);
      CHECK((end.direction() - start.direction()).norm() < 1.e-8);
      REQUIRE(end.covariance());
      // Compare dimensionless correlations/errors after scaling each coordinate
      // by its original standard deviation (time otherwise dominates the norm).
      const Acts::BoundTrackParameters::CovarianceMatrix scale =
          start.covariance()->diagonal().cwiseSqrt().cwiseInverse().asDiagonal();
      CHECK((scale * (*end.covariance() - *start.covariance()) * scale).norm() < 1.e-4);
    }
  }
}

TEST_CASE("B0 local transport handles zero and invalid distances", "[B0SeedTransport]") {
  const auto gctx = geometryContext();
  const Acts::MagneticFieldContext mctx;
  const auto field     = std::make_shared<Acts::ConstantBField>(Acts::Vector3::Zero());
  const auto start     = makeSeed(0.020, 1.);
  const auto unchanged = eicrecon::transportB0Seed(start, 0., field, gctx, mctx);
  REQUIRE(unchanged);
  CHECK(unchanged->parameters() == start.parameters());
  CHECK(*unchanged->covariance() == *start.covariance());
  const auto straight = eicrecon::transportB0Seed(start, 10., field, gctx, mctx);
  REQUIRE(straight);
  CHECK((straight->position(gctx) - start.position(gctx) + 10. * start.direction()).norm() < 1.e-8);
  CHECK((straight->direction() - start.direction()).norm() < 1.e-12);
  CHECK_FALSE(eicrecon::transportB0Seed(start, -1., field, gctx, mctx));
  CHECK_FALSE(eicrecon::transportB0Seed(start, std::numeric_limits<double>::quiet_NaN(), field,
                                        gctx, mctx));
}

TEST_CASE("B0 local transport rejects invalid seeds at any distance", "[B0SeedTransport]") {
  const auto gctx = geometryContext();
  const Acts::MagneticFieldContext mctx;
  const auto field = std::make_shared<Acts::ConstantBField>(Acts::Vector3::Zero());
  const auto seed  = makeSeed(0.020, 1.);
  auto covariance  = *seed.covariance();
  covariance(Acts::eBoundPhi, Acts::eBoundPhi) = std::numeric_limits<double>::infinity();
  const Acts::BoundTrackParameters invalid{seed.referenceSurface().getSharedPtr(),
                                           seed.parameters(), covariance,
                                           seed.particleHypothesis()};
  const Acts::BoundTrackParameters withoutCovariance{seed.referenceSurface().getSharedPtr(),
                                                     seed.parameters(), std::nullopt,
                                                     seed.particleHypothesis()};
  for (double distance : {0., 10.}) {
    CAPTURE(distance);
    CHECK_FALSE(eicrecon::transportB0Seed(invalid, distance, field, gctx, mctx));
    CHECK_FALSE(eicrecon::transportB0Seed(withoutCovariance, distance, field, gctx, mctx));
  }
}
