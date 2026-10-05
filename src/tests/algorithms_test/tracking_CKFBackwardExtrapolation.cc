// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Propagator/ActorList.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/VoidNavigator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <memory>
#include <vector>

#include "algorithms/tracking/CKFTracking.h"

namespace {
using eicrecon::CKFTracking;
using TestPropagator = Acts::Propagator<Acts::EigenStepper<>, Acts::VoidNavigator>;
using TestOptions    = TestPropagator::template Options<Acts::ActorList<Acts::MaterialInteractor>>;

constexpr double fieldTesla  = 1.7;
constexpr double momentumGeV = 1.0;

struct HelixPoint {
  Acts::Vector3 position;
  Acts::Vector3 direction;
};

/// Exact circular transport in a uniform +y field, recording the trajectory
/// where it crosses the ascending target z positions.
std::vector<HelixPoint> transportHelix(const Acts::Vector3& start, const Acts::Vector3& direction,
                                       double charge, double momentum,
                                       const std::vector<double>& targetZ) {
  std::vector<HelixPoint> points{{start, direction.normalized()}};
  Acts::Vector3 position  = start;
  Acts::Vector3 tangent   = direction.normalized();
  const double stepLength = 1. * Acts::UnitConstants::mm;
  for (const double target : targetZ) {
    while (position.z() < target) {
      // dphi = -q*B*ds/p about +y sends +z movers towards -x for q > 0
      const double angle =
          -0.3 * charge * fieldTesla * (stepLength / Acts::UnitConstants::m) / momentum;
      const double cosAngle = std::cos(angle);
      const double sinAngle = std::sin(angle);
      const double x        = tangent.x() * cosAngle + tangent.z() * sinAngle;
      const double z        = tangent.z() * cosAngle - tangent.x() * sinAngle;
      tangent               = Acts::Vector3{x, tangent.y(), z}.normalized();
      position += tangent * stepLength;
    }
    points.push_back({position, tangent});
  }
  return points;
}

/// Closed-form z of the upstream circle-line crossing for a helix through
/// `point` with `direction`: independent oracle for the stepper arrival.
double upstreamCrossingZ(const Acts::Vector3& point, const Acts::Vector3& direction, double charge,
                         double momentum) {
  const double rho = 1000. * momentum / (0.3 * charge * fieldTesla) * Acts::UnitConstants::mm;
  const Acts::Vector3 unit   = direction.normalized();
  const Acts::Vector3 center = point + rho * Acts::Vector3{-unit.z(), 0., unit.x()};
  const double dz            = std::sqrt(rho * rho - center.x() * center.x());
  return center.z() - dz;
}

struct Fixture {
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  std::shared_ptr<const Acts::ConstantBField> field{std::make_shared<Acts::ConstantBField>(
      Acts::Vector3{0., fieldTesla* Acts::UnitConstants::T, 0.})};
  std::shared_ptr<const Acts::Logger> logger{
      Acts::getDefaultLogger("B0BackwardTest", Acts::Logging::INFO)};
  TestPropagator propagator{Acts::EigenStepper<>(field), Acts::VoidNavigator{}, logger};
  std::shared_ptr<const Acts::PerigeeSurface> perigee{
      Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero())};
  std::shared_ptr<Acts::VectorTrackContainer> tracks{
      std::make_shared<Acts::VectorTrackContainer>()};
  std::shared_ptr<Acts::VectorMultiTrajectory> states{
      std::make_shared<Acts::VectorMultiTrajectory>()};
  ActsExamples::TrackContainer container{tracks, states};

  Fixture() {
#if Acts_VERSION_MAJOR >= 45
    gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
    gctx = Acts::GeometryContext{};
#endif
  }

  TestOptions options() const {
    TestOptions options{gctx, mctx};
    options.maxSteps = 10000;
    return options;
  }

  /// Smoothed three-state soft track across B0-like stations. The first
  /// state tilts slightly outward with the turning point downstream, so the
  /// backward arc is turning-free; the last state tilts inward, so the
  /// straight-line intersection misleads downstream like production.
  ActsExamples::TrackProxy makeTrack(double charge, double momentum) {
    const Acts::Vector3 start{700. * Acts::UnitConstants::mm, 0., 5900. * Acts::UnitConstants::mm};
    const Acts::Vector3 direction{0.05, 0., 1.};
    const std::vector<double> stations{6200. * Acts::UnitConstants::mm,
                                       6500. * Acts::UnitConstants::mm};
    auto points = transportHelix(start, direction, 1., momentum, stations);
    auto track  = container.makeTrack();
    for (auto [position, tangent] : points) {
      if (charge < 0.) {
        position.x() *= -1.;
        tangent.x() *= -1.;
      }
      auto transform          = Acts::Transform3::Identity();
      transform.translation() = position;
      auto state              = track.appendTrackState(Acts::TrackStatePropMask::All);
      state.setReferenceSurface(Acts::Surface::makeShared<Acts::PlaneSurface>(transform));
#if Acts_VERSION_MAJOR >= 45
      state.typeFlags().setIsMeasurement();
#else
      state.typeFlags().set(Acts::TrackStateFlag::MeasurementFlag);
#endif
      Acts::BoundVector smoothed   = Acts::BoundVector::Zero();
      smoothed[Acts::eBoundLoc0]   = 0.;
      smoothed[Acts::eBoundLoc1]   = 0.;
      smoothed[Acts::eBoundPhi]    = std::atan2(tangent.y(), tangent.x());
      smoothed[Acts::eBoundTheta]  = std::acos(tangent.z());
      smoothed[Acts::eBoundQOverP] = charge / (momentum * Acts::UnitConstants::GeV);
      smoothed[Acts::eBoundTime]   = 0.;
      state.smoothed()             = smoothed;
#if Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 1)
      Acts::BoundMatrix covariance = Acts::BoundMatrix::Identity();
#else
      Acts::BoundSquareMatrix covariance = Acts::BoundSquareMatrix::Identity();
#endif
      covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)   = 0.0004;
      covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)   = 0.0004;
      covariance(Acts::eBoundPhi, Acts::eBoundPhi)     = 1.e-8;
      covariance(Acts::eBoundTheta, Acts::eBoundTheta) = 1.e-8;
      covariance(Acts::eBoundQOverP, Acts::eBoundQOverP) =
          0.0025 * smoothed[Acts::eBoundQOverP] * smoothed[Acts::eBoundQOverP];
      covariance(Acts::eBoundTime, Acts::eBoundTime) = 100.;
      state.smoothedCovariance()                     = covariance;
    }
    return track;
  }

  double expectedZ0(double charge, double momentum) const {
    const Acts::Vector3 start{charge * 700. * Acts::UnitConstants::mm, 0.,
                              5900. * Acts::UnitConstants::mm};
    const Acts::Vector3 direction{charge * 0.05, 0., 1.};
    return upstreamCrossingZ(start, direction, charge, momentum);
  }
};

void checkCovariance(const auto& covariance) {
  for (std::size_t i = 0; i < 6; ++i) {
    CHECK(std::isfinite(covariance(i, i)));
    CHECK(covariance(i, i) > 0.);
  }
}
} // namespace

TEST_CASE("B0 backward extrapolation reaches the upstream perigee", "[B0Tracking]") {
  for (const double charge : {-1., 1.}) {
    CAPTURE(charge);
    Fixture fixture;
    auto track  = fixture.makeTrack(charge, momentumGeV);
    auto result = CKFTracking::extrapolateBackwardToReference(
        track, *fixture.perigee, fixture.propagator, fixture.options(), *fixture.logger);
    REQUIRE(result.ok());
    const auto parameters = track.parameters();
    CHECK(std::abs(parameters[Acts::eBoundLoc1] - fixture.expectedZ0(charge, momentumGeV)) <
          10. * Acts::UnitConstants::mm);
    CHECK(parameters[Acts::eBoundQOverP] * charge > 0.);
    CHECK(std::abs(std::abs(parameters[Acts::eBoundQOverP]) * momentumGeV - 1.) < 0.01);
    checkCovariance(track.covariance());
  }
}

TEST_CASE("B0 original strategy misleads on bent tracks", "[B0Tracking]") {
  for (const double charge : {-1., 1.}) {
    CAPTURE(charge);
    Fixture fixture;
    auto track        = fixture.makeTrack(charge, momentumGeV);
    const auto choice = Acts::findTrackStateForExtrapolation(
        fixture.gctx, track, *fixture.perigee, Acts::TrackExtrapolationStrategy::firstOrLast,
        *fixture.logger);
    REQUIRE(choice.ok());
    // The straight-line intersection points downstream: the production failure.
    CHECK(choice->second > 0.);
    const auto old = Acts::extrapolateTrackToReferenceSurface(
        track, *fixture.perigee, fixture.propagator, fixture.options(),
        Acts::TrackExtrapolationStrategy::firstOrLast, *fixture.logger);
    // Either it fails outright or it binds far downstream, never upstream.
    const bool correctUpstream = old.ok() && std::abs(track.parameters()[Acts::eBoundLoc1]) <
                                                 5000. * Acts::UnitConstants::mm;
    CHECK(!correctUpstream);
    // A path limit emulating the reduced tracking world kills it outright.
    auto limited             = fixture.options();
    limited.pathLimit        = 500. * Acts::UnitConstants::mm;
    auto track2              = fixture.makeTrack(charge, momentumGeV);
    const auto limitedResult = Acts::extrapolateTrackToReferenceSurface(
        track2, *fixture.perigee, fixture.propagator, limited,
        Acts::TrackExtrapolationStrategy::firstOrLast, *fixture.logger);
    CHECK(!limitedResult.ok());
  }
}

TEST_CASE("B0 backward extrapolation stays accurate when stiff", "[B0Tracking]") {
  constexpr double stiffGeV = 40.;
  for (const double charge : {-1., 1.}) {
    CAPTURE(charge);
    Fixture fixture;
    // A plane reference: grazing-incidence perigee binding over long paths
    // is an ACTS numerics edge beyond this regression's scope, while the
    // stiff long-path transport accuracy is what matters here.
    auto track                   = fixture.makeTrack(charge, stiffGeV);
    const double zRef            = fixture.expectedZ0(charge, stiffGeV);
    auto planeTransform          = Acts::Transform3::Identity();
    planeTransform.translation() = Acts::Vector3{0., 0., zRef};
    const auto plane             = Acts::Surface::makeShared<Acts::PlaneSurface>(planeTransform);
    auto result = CKFTracking::extrapolateBackwardToReference(track, *plane, fixture.propagator,
                                                              fixture.options(), *fixture.logger);
    REQUIRE(result.ok());
    CHECK(std::abs(track.parameters()[Acts::eBoundLoc0]) < 10. * Acts::UnitConstants::mm);
    CHECK(std::abs(track.parameters()[Acts::eBoundLoc1]) < 10. * Acts::UnitConstants::mm);
  }
}

TEST_CASE("B0 backward extrapolation fails cleanly at the step limit", "[B0Tracking]") {
  Fixture fixture;
  auto track        = fixture.makeTrack(1., momentumGeV);
  auto limited      = fixture.options();
  limited.maxSteps  = 1;
  const auto result = CKFTracking::extrapolateBackwardToReference(
      track, *fixture.perigee, fixture.propagator, limited, *fixture.logger);
  CHECK(!result.ok());
}
