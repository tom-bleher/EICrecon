// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldError.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Utilities/Result.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/Vector3f.h>
#include <podio/RelationRange.h>
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "algorithms/tracking/B0TripletSeeding.h"
#include "algorithms/tracking/B0TripletSeedingConfig.h"
#include "algorithms/tracking/B0TripletSeedingHelpers.h"

using eicrecon::B0TripletSeeding;
using eicrecon::B0TripletSeedingConfig;
using namespace eicrecon::b0_triplet_seeding;

namespace {

// Station positions of the B0 tracker on epic main, and a field along +y like that of the B0 magnet
constexpr std::array<double, 4> kStationZ{5900., 6170., 6439., 6709.};
const Acts::Vector3 kField{0., 1.3 * Acts::UnitConstants::T, 0.};

/// Position at z of a track of charge q and momentum p that starts at (0, 0, z0) along +z in kField.
/// The force q v x B points along -x for q > 0, so the track bends towards -x.
Acts::Vector3 helix(double z, double q, double p, double z0 = 5800.) {
  const double radius = p / kField.norm();
  const double dz     = z - z0;
  return {-q * (radius - std::sqrt(radius * radius - dz * dz)), 0., z};
}

Acts::GeometryContext geometryContext() {
#if Acts_VERSION_MAJOR >= 45
  return Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  return Acts::GeometryContext{};
#endif
}

/// A hit in time with a particle at the speed of light from the origin (Acts measures time as a
/// length), with a 30 ps time error
Point point(const Acts::Vector3& position, std::size_t index) {
  return {.position       = position,
          .field          = kField,
          .fieldDirection = kField.normalized(),
          .time           = position.norm(),
          .variance       = 1e-4,
          .index          = index,
          .timeError      = 0.03 * Acts::UnitConstants::ns};
}

void addHit(edm4eic::TrackerHitCollection& hits, const Acts::Vector3& position, float time = 0.F,
            edm4eic::CovDiag3f cov = {1e-4F, 1e-4F, 0.F}) {
  hits.create(hits.size(),
              edm4hep::Vector3f(static_cast<float>(position.x()), static_cast<float>(position.y()),
                                static_cast<float>(position.z())),
              cov, time, 0.F, 1.F, 0.F);
}

/// kField below maxZ, and a failed lookup (as outside a field map) beyond
class TruncatedField final : public Acts::MagneticFieldProvider {
public:
  explicit TruncatedField(double maxZ) : m_maxZ(maxZ) {}
  Cache makeCache(const Acts::MagneticFieldContext& /* mctx */) const override {
    return Cache(std::in_place_type<int>, 0);
  }
  Acts::Result<Acts::Vector3> getField(const Acts::Vector3& position,
                                       Cache& /* cache */) const override {
    if (position.z() > m_maxZ) {
      return Acts::Result<Acts::Vector3>::failure(Acts::MagneticFieldError::OutOfBounds);
    }
    return Acts::Result<Acts::Vector3>::success(kField);
  }

private:
  double m_maxZ;
};

/// The B0 magnet as in DD4hep: a dipole along +y plus a quadrupole gradient, about an axis at
/// x = -145.8 mm, between z = 5800 and 7000 mm (as for epic main and the realistic B0 tracker)
class B0MultipoleField final : public Acts::MagneticFieldProvider {
public:
  Cache makeCache(const Acts::MagneticFieldContext& /* mctx */) const override {
    return Cache(std::in_place_type<int>, 0);
  }
  Acts::Result<Acts::Vector3> getField(const Acts::Vector3& position,
                                       Cache& /* cache */) const override {
    return Acts::Result<Acts::Vector3>::success(at(position));
  }
  static Acts::Vector3 at(const Acts::Vector3& position) {
    constexpr double dipole   = 1.18405 * Acts::UnitConstants::T;
    constexpr double gradient = -8.1224 * Acts::UnitConstants::T / Acts::UnitConstants::m;
    constexpr double axisX    = -145.779;
    if (position.z() < 5800. || position.z() > 7000.) {
      return Acts::Vector3::Zero();
    }
    return {gradient * position.y(), dipole + gradient * (position.x() - axisX), 0.};
  }
};

/// Positions at the given z planes of a track of charge q and momentum p from a start position
/// and direction, integrated through B0MultipoleField with a fourth-order Runge-Kutta-Nystroem
/// stepper (as Acts::EigenStepper). The steps are shortened to end on every plane and on the
/// field edges, so that no step straddles the jump of the field or has to be stepped back: a
/// straight step back keeps the direction of the overshoot, a kink of up to 1e-4 at 2 GeV, and
/// an edge inside a step integrates the field wrongly. Agrees with an adaptive 8th-order
/// integrator to 1.2 um.
std::vector<Acts::Vector3> track(double q, double p, Acts::Vector3 position,
                                 Acts::Vector3 direction, const std::vector<double>& planes) {
  const double qOverP = q / p;
  auto derivative     = [&](const Acts::Vector3& t, const Acts::Vector3& r) -> Acts::Vector3 {
    return qOverP * t.cross(B0MultipoleField::at(r));
  };
  std::vector<double> stops = planes;
  for (const double edge : {5800., 7000.}) {
    if (edge > position.z()) {
      stops.push_back(edge);
    }
  }
  std::ranges::sort(stops);
  std::vector<Acts::Vector3> crossings;
  const double maxStep = 0.5 * Acts::UnitConstants::mm;
  direction.normalize();
  for (const double stop : stops) {
    while (stop - position.z() > 1e-9 * Acts::UnitConstants::mm) {
      const double h         = std::min(maxStep, (stop - position.z()) / direction.z());
      const Acts::Vector3 r0 = position;
      const Acts::Vector3 t0 = direction;
      const Acts::Vector3 k1 = derivative(t0, r0);
      const Acts::Vector3 rm = r0 + 0.5 * h * t0 + h * h / 8. * k1;
      const Acts::Vector3 k2 = derivative(t0 + 0.5 * h * k1, rm);
      const Acts::Vector3 k3 = derivative(t0 + 0.5 * h * k2, rm);
      const Acts::Vector3 k4 = derivative(t0 + h * k3, r0 + h * t0 + 0.5 * h * h * k3);
      position               = r0 + h * t0 + h * h / 6. * (k1 + k2 + k3);
      direction              = (t0 + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4)).normalized();
    }
    if (std::ranges::find(planes, stop) != planes.end()) {
      crossings.push_back(position);
    }
  }
  return crossings;
}

} // namespace

TEST_CASE("B0TripletSeeding groups hits into stations at z gaps", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  addHit(hits, {0., 0., kStationZ[0]});
  addHit(hits, {1., 0., kStationZ[0] + 5.}); // second face of the first station
  addHit(hits, {0., 0., kStationZ[2]});
  addHit(hits, {0., 0., kStationZ[1]});
  addHit(hits, {0., 0., std::numeric_limits<double>::quiet_NaN()});
  addHit(hits, {0., 0., kStationZ[3]}, std::numeric_limits<float>::quiet_NaN());

  auto stations = groupStations(hits, 50. * Acts::UnitConstants::mm);

  // The hits with a NaN position or time are skipped, the rest are sorted by z
  REQUIRE(stations.size() == 3);
  CHECK(stations[0].size() == 2);
  CHECK(stations[1].size() == 1);
  CHECK(stations[2].size() == 1);
  CHECK(stations[1][0].index == 3);
  CHECK(stations[2][0].index == 2);
  CHECK(stations[0][0].variance == Catch::Approx(1e-4));

  // The field is only looked up on the stations that can hold a middle hit
  const Acts::ConstantBField field(kField);
  Acts::MagneticFieldContext mctx;
  auto cache = field.makeCache(mctx);
  CHECK(lookUpField(stations, field, cache) == 0);
  CHECK(stations[0][0].field.isZero());
  CHECK(stations[2][0].field.isZero());
  CHECK(stations[1][0].field.isApprox(kField));
  CHECK(stations[1][0].fieldDirection.isApprox(kField.normalized()));
}

TEST_CASE("B0TripletSeeding takes the hit variance in the sensor plane", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  // The larger in-plane variance is used, and the one along the sensor normal is ignored
  addHit(hits, {0., 0., kStationZ[0]}, 0.F, {1e-4F, 4e-4F, 1.F});
  const auto stations = groupStations(hits, 50. * Acts::UnitConstants::mm);
  REQUIRE(stations.size() == 1);
  CHECK(stations[0][0].variance == Catch::Approx(4e-4));
}

TEST_CASE("B0TripletSeeding counts failed field lookups", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  for (const double z : kStationZ) {
    addHit(hits, {0., 0., z});
  }
  auto stations = groupStations(hits, 50. * Acts::UnitConstants::mm);
  REQUIRE(stations.size() == 4);

  // The field ends between the two middle stations: the hit on the third station has none
  const TruncatedField field(0.5 * (kStationZ[1] + kStationZ[2]));
  Acts::MagneticFieldContext mctx;
  auto cache = field.makeCache(mctx);
  CHECK(lookUpField(stations, field, cache) == 1);
  CHECK(stations[1][0].field.isApprox(kField));
  CHECK(stations[2][0].field.isZero());
  CHECK(stations[2][0].fieldDirection.isZero());
}

TEST_CASE("B0TripletSeeding measures the chord residual along the field", "[B0TripletSeeding]") {
  const auto a = point({0., 0., kStationZ[0]}, 0);
  const auto c = point({0., 0., kStationZ[2]}, 2);
  // An offset of the middle hit across the field (the bending direction) is not a residual
  CHECK(chordResidual(a, point({2., 0., kStationZ[1]}, 1), c) == Catch::Approx(0.).margin(1e-12));
  // An offset along the field is
  CHECK(chordResidual(a, point({0., 2., kStationZ[1]}, 1), c) == Catch::Approx(2.));
}

TEST_CASE("B0TripletSeeding recovers charge and momentum on a helix", "[B0TripletSeeding]") {
  const auto gctx = geometryContext();
  B0TripletSeedingConfig cfg;
  for (const double q : {+1., -1.}) {
    for (const double p : {5., 20., 100.}) {
      const auto a    = point(helix(kStationZ[0], q, p), 0);
      const auto b    = point(helix(kStationZ[1], q, p), 1);
      const auto c    = point(helix(kStationZ[2], q, p), 2);
      const auto seed = makeSeed(cfg, gctx, a, b, c);
      REQUIRE(seed.has_value());
      const double qOverP = seed->parameters[Acts::eBoundQOverP] * Acts::UnitConstants::GeV;
      CHECK(qOverP == Catch::Approx(q / p).epsilon(0.01));
      // The seed sits on a perigee anchorDistance upstream of its first hit
      CHECK(seed->anchor.z() == Catch::Approx(kStationZ[0] - cfg.anchorDistance).margin(0.1));
      CHECK((seed->errors.array() > 0.).all());
    }
  }
}

TEST_CASE("B0TripletSeeding rejects triplets below minMomentum", "[B0TripletSeeding]") {
  const auto gctx = geometryContext();
  B0TripletSeedingConfig cfg;
  cfg.minMomentum = 1. * Acts::UnitConstants::GeV;
  const double p  = 0.5;
  const auto seed =
      makeSeed(cfg, gctx, point(helix(kStationZ[0], 1., p), 0),
               point(helix(kStationZ[1], 1., p), 1), point(helix(kStationZ[2], 1., p), 2));
  CHECK_FALSE(seed.has_value());
}

TEST_CASE("B0TripletSeeding writes seeds in EDM4eic units", "[B0TripletSeeding]") {
  const auto gctx = geometryContext();
  B0TripletSeedingConfig cfg;
  edm4eic::TrackerHitCollection hits;
  std::vector<Point> points;
  for (std::size_t i = 0; i < 3; ++i) {
    const auto position = helix(kStationZ.at(i), 1., 20.);
    addHit(hits, position);
    points.push_back(point(position, i));
  }
  const auto seed = makeSeed(cfg, gctx, points[0], points[1], points[2]);
  REQUIRE(seed.has_value());

  edm4eic::TrackSeedCollection seeds;
  edm4eic::TrackParametersCollection params;
  const Candidate candidate{.residual = 0.25, .points = {&points[0], &points[1], &points[2]}};
  writeSeed(*seed, candidate, hits, seeds, params);

  REQUIRE(seeds.size() == 1);
  REQUIRE(params.size() == 1);
  CHECK(seeds[0].getHits().size() == 3);
  CHECK(seeds[0].getQuality() == Catch::Approx(-0.25));
  CHECK(seeds[0].getPerigee().z == Catch::Approx(kStationZ[0] - 10.).margin(0.1));
  CHECK(params[0].getQOverP() == Catch::Approx(1. / 20.).epsilon(0.01));
  // The time prior of 10 ns is stored as a variance in ns^2
  CHECK(params[0].getCovariance()(5, 5) == Catch::Approx(100.).epsilon(1e-4));

  // Every parameter and its variance is converted from Acts to EDM4eic units
  const auto& par = seed->parameters;
  const auto& err = seed->errors;
  CHECK(params[0].getType() == -1);
  CHECK(params[0].getLoc().a ==
        Catch::Approx(par[Acts::eBoundLoc0] / Acts::UnitConstants::mm * edm4eic::unit::mm)
            .margin(1e-5));
  CHECK(params[0].getLoc().b ==
        Catch::Approx(par[Acts::eBoundLoc1] / Acts::UnitConstants::mm * edm4eic::unit::mm)
            .margin(1e-3));
  CHECK(params[0].getPhi() == Catch::Approx(par[Acts::eBoundPhi]).epsilon(1e-6));
  CHECK(params[0].getTheta() == Catch::Approx(par[Acts::eBoundTheta]).epsilon(1e-6));
  CHECK(params[0].getTime() ==
        Catch::Approx(par[Acts::eBoundTime] / Acts::UnitConstants::ns * edm4eic::unit::ns)
            .epsilon(1e-6));
  const auto cov = params[0].getCovariance();
  CHECK(cov(0, 0) == Catch::Approx(err[Acts::eBoundLoc0] * err[Acts::eBoundLoc0]).epsilon(1e-5));
  CHECK(cov(1, 1) == Catch::Approx(err[Acts::eBoundLoc1] * err[Acts::eBoundLoc1]).epsilon(1e-5));
  CHECK(cov(2, 2) == Catch::Approx(err[Acts::eBoundPhi] * err[Acts::eBoundPhi]).epsilon(1e-5));
  CHECK(cov(3, 3) == Catch::Approx(err[Acts::eBoundTheta] * err[Acts::eBoundTheta]).epsilon(1e-5));
  const double qOverPError = err[Acts::eBoundQOverP] * Acts::UnitConstants::GeV;
  CHECK(cov(4, 4) == Catch::Approx(qOverPError * qOverPError).epsilon(1e-5));
  // The prior is diagonal
  CHECK(cov(0, 1) == 0.F);
  CHECK(cov(4, 5) == 0.F);
  // The seed refers to its parameters and its hits
  CHECK(seeds[0].getParams().getObjectID() == params[0].getObjectID());
  for (std::size_t i = 0; i < 3; ++i) {
    CHECK(seeds[0].getHits()[i].getObjectID() == hits[i].getObjectID());
  }
}

// A valid configuration cannot pass init() here: the test fixture's ActsGeometryProvider has no
// magnetic field, which init() requires
TEST_CASE("B0TripletSeeding rejects an invalid configuration", "[B0TripletSeeding]") {
  B0TripletSeeding algo("B0TripletSeedingTest");
  B0TripletSeedingConfig cfg;
  cfg.maxResidual = -1.;
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("maxResidual"));

  cfg            = B0TripletSeedingConfig{};
  cfg.stationGap = std::numeric_limits<double>::quiet_NaN();
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("stationGap"));

  cfg                = B0TripletSeedingConfig{};
  cfg.anchorDistance = -1.;
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("anchorDistance"));

  cfg          = B0TripletSeedingConfig{};
  cfg.maxSeeds = 0;
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("maxSeeds"));

  cfg           = B0TripletSeedingConfig{};
  cfg.timeError = 0.;
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("timeError"));

  cfg                  = B0TripletSeedingConfig{};
  cfg.timeWindowSigmas = -1.;
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("timeWindowSigmas"));

  cfg                    = B0TripletSeedingConfig{};
  cfg.fieldTermAllowance = -1.;
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("fieldTermAllowance"));

  cfg                = B0TripletSeedingConfig{};
  cfg.timeErrorFloor = 0.;
  algo.applyConfig(cfg);
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("timeErrorFloor"));

  // A valid configuration passes the checks and only fails for the missing field
  algo.applyConfig(B0TripletSeedingConfig{});
  CHECK_THROWS_WITH(algo.init(), Catch::Matchers::ContainsSubstring("magnetic field"));
}

TEST_CASE("B0TripletSeeding keeps triplets within maxResidual", "[B0TripletSeeding]") {
  B0TripletSeedingConfig cfg;
  std::vector<std::vector<Point>> stations;
  for (std::size_t i = 0; i < kStationZ.size(); ++i) {
    stations.push_back({point(helix(kStationZ.at(i), 1., 20.), i)});
  }
  // A track on all four stations gives one candidate per station triple
  CHECK(findCandidates(cfg, stations).size() == 4);

  // A stray hit on the second station, 3 mm along the field off the track. As the first hit of the
  // triplet on the last three stations its residual is 1.5 mm, as a middle hit 3 mm.
  stations[1].push_back(point(helix(kStationZ[1], 1., 20.) + Acts::Vector3{0., 3., 0.}, 4));
  CHECK(findCandidates(cfg, stations).size() == 4);
  cfg.maxResidual = 2. * Acts::UnitConstants::mm;
  CHECK(findCandidates(cfg, stations).size() == 5);
  cfg.maxResidual = 4. * Acts::UnitConstants::mm;
  CHECK(findCandidates(cfg, stations).size() == 7);

  // A hit without a field is never a middle hit
  for (auto& station : stations) {
    for (auto& p : station) {
      p.fieldDirection = Acts::Vector3::Zero();
    }
  }
  CHECK(findCandidates(cfg, stations).empty());
}

TEST_CASE("B0TripletSeeding corrects for the field change along the triplet",
          "[B0TripletSeeding]") {
  const auto gctx = geometryContext();
  B0TripletSeedingConfig cfg;
  const FieldLookup fieldAt = [](const Acts::Vector3& r) { return B0MultipoleField::at(r); };
  // Stations of the realistic B0 tracker. A 2 GeV track 80 mm above the magnet midplane, from
  // upstream of the field along the ion beam (-25 mrad in x), where the quadrupole term is large.
  const std::vector<double> planes{5850., 6272.2, 6601., 6880.};
  for (const double q : {+1., -1.}) {
    const double p  = 2. * Acts::UnitConstants::GeV;
    const auto hits = track(q, p, {-130., 80., 5700.}, {-0.025, 0., 1.}, planes);
    auto pointAt    = [&](std::size_t i) {
      auto point           = eicrecon::b0_triplet_seeding::Point{};
      point.position       = hits.at(i);
      point.field          = B0MultipoleField::at(hits.at(i));
      point.fieldDirection = point.field.normalized();
      point.time           = hits.at(i).norm();
      point.variance       = 1e-4;
      point.index          = i;
      return point;
    };
    // Every station triple: 123, 124, 134, 234
    for (const auto& [i, j, k] :
         std::vector<std::array<std::size_t, 3>>{{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}) {
      CAPTURE(q, i, j, k);
      const auto a        = pointAt(i);
      const auto b        = pointAt(j);
      const auto c        = pointAt(k);
      const auto integral = fieldIntegral(a, b, c, cfg.fieldSamples, fieldAt);
      const auto plain    = makeSeed(cfg, gctx, a, b, c);
      const auto seed     = makeSeed(cfg, gctx, a, b, c, integral.field);
      REQUIRE(plain.has_value());
      REQUIRE(seed.has_value());
      const double qOverP    = q / p;
      const double plainBias = plain->parameters[Acts::eBoundQOverP] / qOverP - 1.;
      const double bias      = seed->parameters[Acts::eBoundQOverP] / qOverP - 1.;
      CAPTURE(plainBias, bias);
      // The effective field removes most of the momentum bias of the field at the middle hit, up
      // to 8 % here. What is left, up to 0.8 % (134, q = +1), comes from sampling the field on the
      // parabola through the hits, which is off the true path by up to a few mm in this gradient
      // (with the field on the true path it is 0.16 %). It is far below the 25 % prior on q/p and
      // below the 1.2-1.9 % spread of the q/p of real triplets from multiple scattering. Where the
      // field at the middle hit is already good (q = -1 bends towards the magnet axis), the bias
      // can be larger than with it, so the improvement is only required where there is one.
      CHECK(std::abs(bias) < 0.01);
      if (std::abs(plainBias) > 0.01) {
        CHECK(std::abs(bias) < 0.2 * std::abs(plainBias));
      }
      // The predicted residual removes the field term from the chord residual, up to 1.8 mm here,
      // to within 0.03 mm. The same parabola leaves up to 0.021 mm (124, q = +1; 0.002 mm with
      // the true path), below the 0.025 mm that the 20 um hit resolution gives the chord residual.
      const double raw = signedChordResidual(a, b, c);
      const double corrected =
          raw - seed->parameters[Acts::eBoundQOverP] * integral.residualPerQOverP;
      CAPTURE(raw, corrected);
      CHECK(std::abs(corrected) < 0.03 * Acts::UnitConstants::mm);
      CHECK(std::abs(corrected) < std::abs(raw));
    }
  }

  // Without the field lookup, or with fieldSamples = 0, the field at the middle hit is used
  std::vector<std::vector<Point>> stations;
  const auto hits = track(1., 2., {-130., 80., 5700.}, {-0.025, 0., 1.}, planes);
  for (std::size_t i = 0; i < hits.size(); ++i) {
    auto point           = Point{};
    point.position       = hits[i];
    point.field          = B0MultipoleField::at(hits[i]);
    point.fieldDirection = point.field.normalized();
    point.time           = hits[i].norm();
    point.variance       = 1e-4;
    point.index          = i;
    stations.push_back({point});
  }
  // With the correction all four triplets pass maxResidual, with residuals below the hit
  // resolution. Without it the triplets that skip a station, 124 and 134 (raw 1.5 and 1.7 mm),
  // fail, and the adjacent ones, 123 and 234, have larger residuals.
  const auto corrected = findCandidates(cfg, stations, fieldAt);
  REQUIRE(corrected.size() == 4);
  auto off         = cfg;
  off.fieldSamples = 0;
  const auto plain = findCandidates(off, stations, fieldAt);
  REQUIRE(plain.size() == 2);
  for (const auto& candidate : plain) {
    CHECK(candidate.field.isApprox(candidate.points[1]->field));
    CHECK(candidate.points[2]->index - candidate.points[0]->index == 2);
  }
  for (const auto& candidate : corrected) {
    CAPTURE(candidate.points[0]->index, candidate.points[1]->index, candidate.points[2]->index);
    CHECK(candidate.residual < 0.03 * Acts::UnitConstants::mm);
    for (const auto& other : plain) {
      if (other.points == candidate.points) {
        CHECK(candidate.residual < other.residual);
      }
    }
  }
}

TEST_CASE("B0TripletSeeding keeps hits within the time window", "[B0TripletSeeding]") {
  B0TripletSeedingConfig cfg;
  const double ns = Acts::UnitConstants::ns;
  // Two hits 1030 mm apart, as the first and the last station of the realistic B0 tracker, with
  // 30 ps time errors. Acts measures time as a length, so a hit in time with light has t = z.
  auto hit = [&](double z, double lag, double error = 0.03) {
    auto p      = point({0., 0., z}, 0);
    p.time      = z + lag * ns;
    p.timeError = error * ns;
    return p;
  };
  const double z0  = 5850.;
  const double z1  = z0 + 1030.;
  const auto first = hit(z0, 0.);
  // A light-speed hit passes, and so does a 1.2 GeV proton: beta = 0.79, a lag of 0.92 ns
  CHECK(inTimeWindow(cfg, first, hit(z1, 0.)));
  CHECK(inTimeWindow(cfg, first, hit(z1, 0.92)));
  // At the default minMomentum of 1 GeV the window ends at the lag of a 1 GeV proton, 1.28 ns,
  // plus 5 sigma of 0.042 ns
  CHECK(inTimeWindow(cfg, first, hit(z1, 1.45)));
  CHECK_FALSE(inTimeWindow(cfg, first, hit(z1, 1.55)));
  // Nothing is faster than light: a hit 0.3 ns early is rejected, unlike one within the errors
  CHECK_FALSE(inTimeWindow(cfg, first, hit(z1, -0.3)));
  CHECK(inTimeWindow(cfg, first, hit(z1, -0.1)));
  // The window scales with the distance: over 422 mm a 1 GeV proton lags only 0.52 ns
  CHECK(inTimeWindow(cfg, first, hit(z0 + 422., 0.5)));
  CHECK_FALSE(inTimeWindow(cfg, first, hit(z0 + 422., 0.9)));
  // A heavier slow particle widens it: a deuteron at 1 GeV per unit charge lags 3.87 ns
  auto deuteron             = cfg;
  deuteron.slowParticleMass = 1.876 * Acts::UnitConstants::GeV;
  CHECK(inTimeWindow(deuteron, first, hit(z1, 3.85)));
  CHECK_FALSE(inTimeWindow(cfg, first, hit(z1, 3.85)));
  // With 8 ns hit time errors, as before the B0 timing was set to 30 ps, the window rejects nothing
  CHECK(inTimeWindow(cfg, hit(z0, 0., 8.), hit(z1, 30., 8.)));
  CHECK(inTimeWindow(cfg, hit(z0, 0., 8.), hit(z1, -30., 8.)));
  // A hit without a valid time error is never rejected
  CHECK(inTimeWindow(cfg, first, hit(z1, 30., std::numeric_limits<double>::quiet_NaN())));
  // and the window can be turned off
  auto off             = cfg;
  off.timeWindowSigmas = 0.;
  CHECK(inTimeWindow(off, first, hit(z1, -30.)));
}

TEST_CASE("B0TripletSeeding rejects triplets out of time", "[B0TripletSeeding]") {
  B0TripletSeedingConfig cfg;
  std::vector<std::vector<Point>> stations;
  for (std::size_t i = 0; i < 3; ++i) {
    stations.push_back({point(helix(kStationZ.at(i), 1., 20.), i)});
  }
  CHECK(findCandidates(cfg, stations).size() == 1);
  // A middle or last hit from another bunch crossing, 10 ns later or earlier, is rejected
  for (const std::size_t shifted : {1, 2}) {
    for (const double offset : {-10., 10.}) {
      auto other = stations;
      other[shifted][0].time += offset * Acts::UnitConstants::ns;
      CHECK(findCandidates(cfg, other).empty());
    }
  }
}

TEST_CASE("B0TripletSeeding ranks candidates and stops at maxSeeds", "[B0TripletSeeding]") {
  std::vector<Point> points;
  for (std::size_t i = 0; i < 6; ++i) {
    points.push_back(point({0., 0., 0.}, i));
  }
  auto candidate = [&](double residual, std::size_t first) {
    return Candidate{.residual = residual,
                     .points   = {&points[first], &points[first + 1], &points[first + 2]}};
  };
  B0TripletSeedingConfig cfg;

  SECTION("in order of residual, ties broken by the hit indices") {
    std::vector<Candidate> candidates{candidate(0.3, 0), candidate(0.1, 2), candidate(0.2, 1),
                                      candidate(0.1, 0)};
    std::vector<double> residuals;
    std::vector<std::size_t> firsts;
    const auto left = selectSeeds(cfg, candidates, [&](const Candidate& c) {
      residuals.push_back(c.residual);
      firsts.push_back(c.points[0]->index);
      return true;
    });
    CHECK(left == 0);
    CHECK(residuals == std::vector<double>{0.1, 0.1, 0.2, 0.3});
    CHECK(firsts == std::vector<std::size_t>{0, 2, 1, 0});
  }

  SECTION("stopping once maxSeeds seeds were made, not counting failed ones") {
    cfg.maxSeeds = 2;
    std::vector<Candidate> candidates{candidate(0.4, 0), candidate(0.1, 0), candidate(0.3, 0),
                                      candidate(0.2, 0)};
    std::vector<double> made;
    const auto left = selectSeeds(cfg, candidates, [&](const Candidate& c) {
      // The best candidate fails to make a seed, as makeSeed can
      if (c.residual < 0.15) {
        return false;
      }
      made.push_back(c.residual);
      return true;
    });
    CHECK(made == std::vector<double>{0.2, 0.3});
    CHECK(left == 1);
  }
}

TEST_CASE("B0TripletSeeding widens the position prior by the anchor sagitta",
          "[B0TripletSeeding]") {
  const auto gctx = geometryContext();
  B0TripletSeedingConfig cfg;
  const double p = 5.;
  // Start the track far upstream of the first station, so that it is curved at the anchor
  const double z0 = 5000.;
  const auto a    = point(helix(kStationZ[0], 1., p, z0), 0);
  const auto b    = point(helix(kStationZ[1], 1., p, z0), 1);
  const auto c    = point(helix(kStationZ[2], 1., p, z0), 2);

  const auto near = makeSeed(cfg, gctx, a, b, c);
  REQUIRE(near.has_value());
  // At the default anchorDistance the sagitta is negligible
  CHECK(near->errors[Acts::eBoundLoc0] == Catch::Approx(cfg.positionError).epsilon(1e-3));

  cfg.anchorDistance = 500. * Acts::UnitConstants::mm;
  const auto far     = makeSeed(cfg, gctx, a, b, c);
  REQUIRE(far.has_value());
  // The anchor misses the track by about the sagitta, and the prior covers it
  const double radius  = p * Acts::UnitConstants::GeV / kField.norm();
  const double sagitta = 0.5 * cfg.anchorDistance * cfg.anchorDistance / radius;
  const double miss    = std::abs(far->anchor.x() - helix(far->anchor.z(), 1., p, z0).x());
  CHECK(miss == Catch::Approx(sagitta).epsilon(0.1));
  CHECK(far->errors[Acts::eBoundLoc0] ==
        Catch::Approx(std::hypot(cfg.positionError, sagitta)).epsilon(0.02));
}
