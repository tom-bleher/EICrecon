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

Point point(const Acts::Vector3& position, std::size_t index) {
  return {.position       = position,
          .field          = kField,
          .fieldDirection = kField.normalized(),
          .time           = 0.,
          .variance       = 1e-4,
          .index          = index};
}

void addHit(edm4eic::TrackerHitCollection& hits, const Acts::Vector3& position, float time = 0.f) {
  hits.create(hits.size(),
              edm4hep::Vector3f(static_cast<float>(position.x()), static_cast<float>(position.y()),
                                static_cast<float>(position.z())),
              edm4eic::CovDiag3f(1e-4f, 1e-4f, 0.f), time, 0.f, 1.f, 0.f);
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
    const auto position = helix(kStationZ[i], 1., 20.);
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
}

// A valid configuration cannot be tested here: the test fixture's ActsGeometryProvider has no
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
}
