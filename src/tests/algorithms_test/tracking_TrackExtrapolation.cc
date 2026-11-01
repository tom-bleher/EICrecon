// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Direction.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Material/MaterialInteraction.hpp>
#include <Acts/Propagator/PropagatorResult.hpp>
#include <Acts/Utilities/Holders.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/EventData/BoundTrackParameters.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/TrackContainer.hpp>
#include <Acts/EventData/TrackStatePropMask.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Geometry/CylinderVolumeBounds.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MultiRangeBField.hpp>
#include <Acts/Propagator/ActorList.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Propagator/VoidNavigator.hpp>
#include <Acts/Surfaces/CurvilinearSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/RangeXD.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <any>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "algorithms/tracking/TrackExtrapolation.h"

namespace {

using Acts::UnitConstants::GeV;
using Acts::UnitConstants::m;
using Acts::UnitConstants::mm;
using Acts::UnitConstants::T;

using Propagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
using PropagatorOptions =
    Propagator::Options<Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>;
using FreePropagator        = Acts::Propagator<Acts::EigenStepper<>, Acts::VoidNavigator>;
using FreePropagatorOptions = FreePropagator::Options<>;

// Like the B0 tracker: measurements at 5.85 and 6.88 m and, by default, a dipole that reaches past
// the end of the tracking world at z = 7 m. The default dipole starts behind the first measurement,
// so the extrapolation back from it is exact.
constexpr double dipoleBy = 1.3 * T;
constexpr double firstZ   = 5.85 * m;
constexpr double lastZ    = 6.88 * m;

struct Setup {
  Acts::GeometryContext gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  Acts::MagneticFieldContext mctx;
  Propagator propagator;
  PropagatorOptions options;
  FreePropagator free;
  FreePropagatorOptions freeOptions;

  explicit Setup(double worldHalfZ = 7. * m, double dipoleZ = 5.9 * m)
      : propagator(makeStepper(dipoleZ), makeNavigator(worldHalfZ),
                   Acts::getDefaultLogger("Propagator", Acts::Logging::INFO))
      , options(gctx, mctx)
      , free(makeStepper(dipoleZ), Acts::VoidNavigator(),
             Acts::getDefaultLogger("FreePropagator", Acts::Logging::INFO))
      , freeOptions(gctx, mctx) {
    freeOptions.surfaceTolerance = 1e-6 * mm;
    freeOptions.pathLimit        = 20. * m;
  }

  static Acts::EigenStepper<> makeStepper(double dipoleZ) {
    // No field anywhere but in the dipole, the later range
    std::vector<std::pair<Acts::RangeXD<3, double>, Acts::Vector3>> ranges{
        {Acts::RangeXD<3, double>(), Acts::Vector3::Zero()},
        {Acts::RangeXD<3, double>({-1. * m, -1. * m, dipoleZ}, {1. * m, 1. * m, 10. * m}),
         Acts::Vector3{0., dipoleBy, 0.}}};
    return Acts::EigenStepper<>(std::make_shared<Acts::MultiRangeBField>(std::move(ranges)));
  }

  static Acts::Navigator makeNavigator(double worldHalfZ) {
    // The full constructor, which creates the boundary surfaces the navigator ends the world at
    auto world = std::make_shared<Acts::TrackingVolume>(
        Acts::Transform3::Identity(),
        std::make_shared<Acts::CylinderVolumeBounds>(0., 1. * m, worldHalfZ), nullptr, nullptr,
        nullptr, Acts::MutableTrackingVolumeVector{}, "World");
    return Acts::Navigator({.trackingGeometry = std::make_shared<Acts::TrackingGeometry>(world)});
  }

  /// Parameters of a pion at a point on the perigee surface through it
  Acts::BoundTrackParameters start(const Acts::Vector3& origin, double theta, double phi,
                                   double qOverP) const {
    Acts::BoundMatrix cov = Acts::BoundMatrix::Identity() * 1e-6;
    return Acts::BoundTrackParameters::create(
               gctx, Acts::Surface::makeShared<Acts::PerigeeSurface>(origin),
               Acts::Vector4{origin.x(), origin.y(), origin.z(), 0.},
               Acts::Vector3{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                             std::cos(theta)},
               qOverP, cov, Acts::ParticleHypothesis::pion())
        .value();
  }

  /// Smoothed track with measurements on planes at the given z, from exact propagation of the start
  template <typename track_container_t>
  auto makeTrack(track_container_t& tracks, const Acts::BoundTrackParameters& start,
                 const std::vector<double>& zs) const {
    auto track = tracks.makeTrack();
    track.setParticleHypothesis(start.particleHypothesis());
    for (double z : zs) {
      auto plane =
          Acts::CurvilinearSurface(Acts::Vector3{0., 0., z}, Acts::Vector3::UnitZ()).planeSurface();
      auto state = propagator.propagate(start, *plane, options).value().endParameters.value();
      auto ts    = track.appendTrackState(Acts::TrackStatePropMask::Smoothed);
      ts.setReferenceSurface(plane);
      ts.smoothed()           = state.parameters();
      ts.smoothedCovariance() = state.covariance().value();
      ts.typeFlags().setIsMeasurement();
    }
    return track;
  }
};

} // namespace

TEST_CASE("A track bending towards the beam line reaches its origin in two stages",
          "[TrackExtrapolation]") {
  // The dipole starts before the first measurement, where a negative track already points towards
  // the beam line, and the world reaches past where the track crosses it downstream
  Setup setup(12. * m, 5.8 * m);
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{}, Acts::VectorMultiTrajectory{}};
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  auto plane   = Acts::CurvilinearSurface(Acts::Vector3{0., 0., 5.79 * m}, Acts::Vector3::UnitZ())
                     .planeSurface();
  const std::vector<double> zs{firstZ, 6.27 * m, lastZ};

  for (const double qOverP : {-1. / (2. * GeV), 1. / (2. * GeV), -1. / (20. * GeV)}) {
    CAPTURE(qOverP);
    auto origin = setup.start(Acts::Vector3::Zero(), 0.025, std::numbers::pi, qOverP);
    // The ACTS helper expresses a negative track where it crosses the beam line downstream
    if (qOverP < 0. && qOverP < -0.1 / GeV) {
      auto helper = setup.makeTrack(tracks, origin, zs);
      REQUIRE(Acts::extrapolateTrackToReferenceSurface(
                  helper, *perigee, setup.propagator, setup.options,
                  Acts::TrackExtrapolationStrategy::firstOrLast)
                  .ok());
      CHECK(helper.loc1() > lastZ);
    }

    auto track  = setup.makeTrack(tracks, origin, zs);
    auto result = eicrecon::extrapolateTrackInTwoStages(
        track, *plane, *perigee, setup.propagator, setup.options, setup.free, setup.freeOptions);
    REQUIRE(result.ok());
    CHECK(*result == eicrecon::TwoStageExtrapolation::Target);
    CHECK(&track.referenceSurface() == perigee.get());
    // Numerical stepping through the dipole, at the default stepper tolerance. The position along
    // the nearly parallel beam line magnifies an angle error by the distance over tan(theta):
    // 1e-4 rad gives 20 mm here, far below the 150-650 mm to which B0 tracks determine it
    CHECK(track.loc0() == Catch::Approx(0.).margin(1e-3 * mm));
    CHECK(track.loc1() == Catch::Approx(0.).margin(50. * mm));
    CHECK(track.phi() == Catch::Approx(origin.phi()).margin(1e-4));
    CHECK(track.theta() == Catch::Approx(origin.theta()).margin(1e-4));
    CHECK(track.qOverP() == Catch::Approx(origin.qOverP()).epsilon(1e-4));
    CHECK(track.covariance().allFinite());
  }
}

TEST_CASE("A track that cannot reach the perigee is kept on the plane", "[TrackExtrapolation]") {
  Setup setup(12. * m, 5.8 * m);
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{}, Acts::VectorMultiTrajectory{}};
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  auto plane   = Acts::CurvilinearSurface(Acts::Vector3{0., 0., 5.79 * m}, Acts::Vector3::UnitZ())
                     .planeSurface();
  // A decay product from 0.3 m off the beam line, which passes the beam line only 10 m upstream of
  // the origin, beyond the path limit
  auto origin =
      setup.start(Acts::Vector3{-0.3 * m, 0., 5. * m}, 0.02, std::numbers::pi, -1. / (2. * GeV));
  auto track            = setup.makeTrack(tracks, origin, {firstZ, lastZ});
  auto freeOptions      = setup.freeOptions;
  freeOptions.pathLimit = 8. * m;
  auto result = eicrecon::extrapolateTrackInTwoStages(track, *plane, *perigee, setup.propagator,
                                                      setup.options, setup.free, freeOptions);
  REQUIRE(result.ok());
  CHECK(*result == eicrecon::TwoStageExtrapolation::ReferencePlane);
  CHECK(&track.referenceSurface() == plane.get());
  // On the plane the track is where it crossed it, with its momentum there
  const auto atPlane = track.createParametersAtReference();
  CHECK(atPlane.position(setup.gctx).z() == Catch::Approx(5.79 * m).margin(1e-3 * mm));
  CHECK(track.qOverP() == Catch::Approx(origin.qOverP()).epsilon(1e-6));
  CHECK(track.covariance().allFinite());
}

TEST_CASE("The plane must lie upstream of the first measurement", "[TrackExtrapolation]") {
  Setup setup(12. * m, 5.8 * m);
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{}, Acts::VectorMultiTrajectory{}};
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  auto plane   = Acts::CurvilinearSurface(Acts::Vector3{0., 0., 6. * m}, Acts::Vector3::UnitZ())
                     .planeSurface();
  auto origin  = setup.start(Acts::Vector3::Zero(), 0.025, std::numbers::pi, 1. / (2. * GeV));
  auto track   = setup.makeTrack(tracks, origin, {firstZ, lastZ});
  const auto before = track.parameters();
  auto result = eicrecon::extrapolateTrackInTwoStages(track, *plane, *perigee, setup.propagator,
                                                      setup.options, setup.free, setup.freeOptions);
  CHECK_FALSE(result.ok());
  CHECK(track.parameters() == before);
}
