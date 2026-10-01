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
#include <boost/container/detail/std_fwd.hpp>
#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/GenericBoundTrackParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#endif
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

// Like the B0 tracker: measurements at 5.85 and 6.88 m and, by default, a dipole that reaches past
// the end of the tracking world at z = 7 m. The default dipole starts behind the first measurement,
// so the extrapolation back from it is exact.
constexpr double dipoleBy = 1.3 * T;
constexpr double firstZ   = 5.85 * m;
constexpr double lastZ    = 6.88 * m;

struct Setup {
#if Acts_VERSION_MAJOR >= 45
  Acts::GeometryContext gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  Acts::GeometryContext gctx;
#endif
  Acts::MagneticFieldContext mctx;
  Propagator propagator;
  PropagatorOptions options;

  explicit Setup(double worldHalfZ = 7. * m, double dipoleZ = 5.9 * m)
      : propagator(makeStepper(dipoleZ), makeNavigator(worldHalfZ),
                   Acts::getDefaultLogger("Propagator", Acts::Logging::INFO))
      , options(gctx, mctx) {}

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
#if Acts_VERSION_MAJOR >= 45
      ts.typeFlags().setIsMeasurement();
#else
      ts.typeFlags().set(Acts::TrackStateFlag::MeasurementFlag);
#endif
    }
    return track;
  }
};

} // namespace

TEST_CASE("A track bending towards the beam line is extrapolated to its origin",
          "[TrackExtrapolation]") {
  Setup setup;
  // A negative pion bends towards the beam line: its tangent at the last measurement meets
  // the beam line behind the end of the world
  auto origin = setup.start(Acts::Vector3::Zero(), 0.025, std::numbers::pi, -1. / (2. * GeV));
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{}, Acts::VectorMultiTrajectory{}};
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());

  auto helper = setup.makeTrack(tracks, origin, {firstZ, lastZ});
  REQUIRE_FALSE(
      Acts::extrapolateTrackToReferenceSurface(helper, *perigee, setup.propagator, setup.options,
                                               Acts::TrackExtrapolationStrategy::firstOrLast)
          .ok());

  auto track = setup.makeTrack(tracks, origin, {firstZ, lastZ});
  REQUIRE(eicrecon::extrapolateTrackToPerigee(track, *perigee, origin.referenceSurface(),
                                              setup.propagator, setup.options,
                                              Acts::getDummyLogger())
              .ok());
  CHECK(&track.referenceSurface() == perigee.get());
  CHECK(track.loc0() == Catch::Approx(0.).margin(1e-3 * mm));
  CHECK(track.loc1() == Catch::Approx(0.).margin(1e-3 * mm));
  CHECK(track.phi() == Catch::Approx(origin.phi()).margin(1e-6));
  CHECK(track.theta() == Catch::Approx(origin.theta()).margin(1e-6));
  CHECK(track.qOverP() == Catch::Approx(origin.qOverP()).epsilon(1e-6));
  CHECK(track.covariance().allFinite());
}

TEST_CASE("A track bending towards the beam line from its first measurement is extrapolated to its "
          "origin",
          "[TrackExtrapolation]") {
  // The dipole starts before the first measurement, where the track already points towards the beam
  // line, and the world reaches past where the track crosses it downstream
  Setup setup(12. * m, 5.8 * m);
  auto origin = setup.start(Acts::Vector3::Zero(), 0.025, std::numbers::pi, -1. / (2. * GeV));
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{}, Acts::VectorMultiTrajectory{}};
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  const std::vector<double> zs{6.27 * m, 6.6 * m, lastZ};

  // The ACTS helper expresses the track where it crosses the beam line after its measurements
  auto helper = setup.makeTrack(tracks, origin, zs);
  REQUIRE(Acts::extrapolateTrackToReferenceSurface(helper, *perigee, setup.propagator,
                                                   setup.options,
                                                   Acts::TrackExtrapolationStrategy::firstOrLast)
              .ok());
  CHECK(helper.loc1() > lastZ);

  auto track = setup.makeTrack(tracks, origin, zs);
  REQUIRE(eicrecon::extrapolateTrackToPerigee(track, *perigee, origin.referenceSurface(),
                                              setup.propagator, setup.options,
                                              Acts::getDummyLogger())
              .ok());
  CHECK(&track.referenceSurface() == perigee.get());
  // Numerical stepping through the dipole on both sides; the position along the nearly parallel
  // beam line magnifies the transverse precision by 1/tan(theta)
  CHECK(track.loc0() == Catch::Approx(0.).margin(1e-3 * mm));
  CHECK(track.loc1() == Catch::Approx(0.).margin(10. * mm));
  CHECK(track.phi() == Catch::Approx(origin.phi()).margin(1e-4));
  CHECK(track.theta() == Catch::Approx(origin.theta()).margin(1e-4));
  CHECK(track.qOverP() == Catch::Approx(origin.qOverP()).epsilon(1e-4));
}

TEST_CASE("A track the helper extrapolates is left to it", "[TrackExtrapolation]") {
  Setup setup;
  // A positive pion bends away from the beam line
  auto origin = setup.start(Acts::Vector3::Zero(), 0.025, std::numbers::pi, 1. / (2. * GeV));
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{}, Acts::VectorMultiTrajectory{}};
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());

  auto helper = setup.makeTrack(tracks, origin, {firstZ, lastZ});
  REQUIRE(Acts::extrapolateTrackToReferenceSurface(helper, *perigee, setup.propagator,
                                                   setup.options,
                                                   Acts::TrackExtrapolationStrategy::firstOrLast)
              .ok());

  auto track = setup.makeTrack(tracks, origin, {firstZ, lastZ});
  REQUIRE(eicrecon::extrapolateTrackToPerigee(track, *perigee, origin.referenceSurface(),
                                              setup.propagator, setup.options,
                                              Acts::getDummyLogger())
              .ok());
  CHECK(track.parameters() == helper.parameters());
}

TEST_CASE("A displaced track is kept on the fallback surface", "[TrackExtrapolation]") {
  Setup setup;
  // A decay product from 0.3 m off the beam line that passes the beam line behind the world
  auto origin =
      setup.start(Acts::Vector3{-0.3 * m, 0., 5. * m}, 0.02, std::numbers::pi, -1. / (2. * GeV));
  Acts::TrackContainer tracks{Acts::VectorTrackContainer{}, Acts::VectorMultiTrajectory{}};
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());

  auto track = setup.makeTrack(tracks, origin, {firstZ, lastZ});
  REQUIRE_FALSE(
      eicrecon::extrapolateTrackBackward(track, *perigee, setup.propagator, setup.options).ok());
  REQUIRE(eicrecon::extrapolateTrackToPerigee(track, *perigee, origin.referenceSurface(),
                                              setup.propagator, setup.options,
                                              Acts::getDummyLogger())
              .ok());
  CHECK(&track.referenceSurface() == &origin.referenceSurface());
  CHECK(track.loc0() == Catch::Approx(0.).margin(1e-3 * mm));
  CHECK(track.loc1() == Catch::Approx(0.).margin(1e-3 * mm));
  CHECK(track.qOverP() == Catch::Approx(origin.qOverP()).epsilon(1e-6));
}
