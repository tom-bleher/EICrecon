// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Utilities/Result.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/Vector3f.h>
#include <podio/RelationRange.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "algorithms/tracking/B0TripletSeeding.h"
#include "algorithms/tracking/B0TripletSeedingConfig.h"

namespace {

// Hits and field in EDM4eic units (mm, GeV, ns) and Tesla
constexpr double fieldY    = 1.3;
constexpr double hitSigma  = 0.02;
constexpr double stationZ0 = 5900.;
constexpr double stationDz = 250.;

std::shared_ptr<Acts::ConstantBField> field() {
  return std::make_shared<Acts::ConstantBField>(
      Acts::Vector3{0., fieldY * Acts::UnitConstants::T, 0.});
}

/// Curvature in 1/mm of a track across the field; its direction turns by -k per mm
double curvature(double qOverP) { return qOverP * fieldY * Acts::UnitConstants::T; }

/// Hit at a distance dz downstream of the first station for a track starting there
/// with angle theta in the x-z plane, perpendicular to the field
edm4hep::Vector3f helixHitAt(double dz, double qOverP, double theta) {
  const double k     = curvature(qOverP);
  const double angle = std::asin(std::sin(theta) - k * dz);
  return {static_cast<float>((std::cos(angle) - std::cos(theta)) / k), 0.F,
          static_cast<float>(stationZ0 + dz)};
}

/// Hit on the plane of a station
edm4hep::Vector3f helixHit(std::size_t station, double qOverP, double theta) {
  return helixHitAt(station * stationDz, qOverP, theta);
}

/// Exactly collinear hits
edm4hep::Vector3f straightHit(std::size_t station) {
  return {10.F * station, 0.F, static_cast<float>(stationZ0 + station * stationDz)};
}

void addHit(edm4eic::TrackerHitCollection& hits, std::uint64_t cellID,
            const edm4hep::Vector3f& position) {
  const auto variance = static_cast<float>(hitSigma * hitSigma);
  hits.create(cellID, position, edm4eic::CovDiag3f{variance, variance, 0.F}, 20.F, 0.1F, 1.e-4F,
              0.F);
}

/// Dipole field along y with a quadrupole gradient, like the B0 magnet: B = (g y, B0 + g x, 0)
class QuadrupoleField final : public Acts::MagneticFieldProvider {
public:
  static constexpr double dipole   = 1.184 * Acts::UnitConstants::T;
  static constexpr double gradient = -8.12 * Acts::UnitConstants::T / Acts::UnitConstants::m;

  Cache makeCache(const Acts::MagneticFieldContext& /* mctx */) const override {
    // The field needs no cache
    return Cache(std::in_place_type<int>, 0);
  }
  Acts::Result<Acts::Vector3> getField(const Acts::Vector3& position,
                                       Cache& /* cache */) const override {
    return Acts::Result<Acts::Vector3>::success(fieldAt(position));
  }
  static Acts::Vector3 fieldAt(const Acts::Vector3& position) {
    return {gradient * position.y(), dipole + gradient * position.x(), 0.};
  }
};

/// Hits on planes at the given z of a track through the quadrupole field, from Runge-Kutta steps
std::vector<edm4hep::Vector3f> quadrupoleHits(Acts::Vector3 position, Acts::Vector3 direction,
                                              double qOverP, const std::vector<double>& zs) {
  auto derivative = [&](const Acts::Vector3& p, const Acts::Vector3& d) {
    return std::pair{d, Acts::Vector3{qOverP * d.cross(QuadrupoleField::fieldAt(p))}};
  };
  std::vector<edm4hep::Vector3f> hits;
  for (const double z : zs) {
    while (position.z() < z) {
      const double h      = std::min(1., (z - position.z()) / direction.z());
      const auto [p1, d1] = derivative(position, direction);
      const auto [p2, d2] = derivative(position + 0.5 * h * p1, direction + 0.5 * h * d1);
      const auto [p3, d3] = derivative(position + 0.5 * h * p2, direction + 0.5 * h * d2);
      const auto [p4, d4] = derivative(position + h * p3, direction + h * d3);
      position += h / 6. * (p1 + 2. * p2 + 2. * p3 + p4);
      direction = (direction + h / 6. * (d1 + 2. * d2 + 2. * d3 + d4)).normalized();
    }
    hits.emplace_back(static_cast<float>(position.x()), static_cast<float>(position.y()),
                      static_cast<float>(position.z()));
  }
  return hits;
}

struct Result {
  edm4eic::TrackSeedCollection seeds;
  edm4eic::TrackParametersCollection params;
};

Result run(const edm4eic::TrackerHitCollection& hits,
           const eicrecon::B0TripletSeedingConfig& cfg                      = {},
           std::shared_ptr<const Acts::MagneticFieldProvider> fieldProvider = field()) {
  eicrecon::B0TripletSeeding algo("B0TripletSeeding", std::move(fieldProvider));
  algo.applyConfig(cfg);
  algo.init();
  Result result;
  algo.process({&hits}, {&result.seeds, &result.params});
  return result;
}

} // namespace

TEST_CASE("B0 seeds recover the curvature of either charge", "[B0TripletSeeding]") {
  for (const double qOverP : {-0.1, 0.1}) {
    CAPTURE(qOverP);
    edm4eic::TrackerHitCollection hits;
    for (std::size_t station = 0; station < 3; ++station) {
      addHit(hits, station + 1, helixHit(station, qOverP, 0.01));
    }
    const auto result = run(hits);
    REQUIRE(result.seeds.size() == 1);
    REQUIRE(result.params.size() == 1);
    const auto seed = result.seeds[0];
    CHECK(seed.getHits().size() == 3);
    CHECK(seed.getParams().getQOverP() == Catch::Approx(qOverP).epsilon(0.02));
    // The default anchor is 10 mm upstream of the first hit, where the helix
    // direction differs from the one at the first hit by 10 mm * curvature
    CHECK(seed.getPerigee().z == Catch::Approx(stationZ0 - 10.).margin(0.01));
    const double anchorTheta = std::asin(std::sin(0.01) + curvature(qOverP) * 10.);
    CHECK(seed.getParams().getTheta() == Catch::Approx(anchorTheta).margin(1.e-5));
    const auto cov = seed.getParams().getCovariance();
    CHECK(std::ranges::all_of(cov.covariance, [](float v) { return std::isfinite(v); }));
    CHECK(cov(4, 4) > 0.F);
  }
}

TEST_CASE("B0 straight triplet keeps the curvature resolution of its hits", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  for (std::size_t station = 0; station < 3; ++station) {
    addHit(hits, station + 1, straightHit(station));
  }
  const auto result = run(hits);
  REQUIRE(result.params.size() == 1);
  const auto params = result.params[0];
  CHECK(std::abs(params.getQOverP()) < 1.e-6F);
  // Sagitta error of three equal hits over two equal lever arms d
  const double d        = std::hypot(stationDz, 10.);
  const double sagitta  = hitSigma * std::sqrt(1.5);
  const double expected = 2. * sagitta / (d * d) / (fieldY * Acts::UnitConstants::T);
  CHECK(std::sqrt(params.getCovariance()(4, 4)) == Catch::Approx(expected).epsilon(0.01));
}

TEST_CASE("B0 seeds need hits on three stations", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection twoStations;
  addHit(twoStations, 1, helixHit(0, 0.1, 0.01));
  addHit(twoStations, 2, helixHit(1, 0.1, 0.01));
  // Front and back sensors of one station are closer than the station gap
  auto back = helixHit(1, 0.1, 0.01);
  back.z += 20.F;
  addHit(twoStations, 3, back);
  CHECK(run(twoStations).seeds.empty());

  edm4eic::TrackerHitCollection fourStations;
  for (std::size_t station = 0; station < 4; ++station) {
    addHit(fourStations, station + 1, helixHit(station, 0.1, 0.01));
  }
  CHECK(run(fourStations).seeds.size() == 4);
}

TEST_CASE("B0 seeds need the middle hit on the chord along the field", "[B0TripletSeeding]") {
  // Only the offset along the field counts; the offset across it is the bending
  for (const auto& [offset, seeds] : {std::pair{edm4hep::Vector3f{0.F, 0.9F, 0.F}, 1UL},
                                      std::pair{edm4hep::Vector3f{0.F, 1.1F, 0.F}, 0UL},
                                      std::pair{edm4hep::Vector3f{1.1F, 0.F, 0.F}, 1UL}}) {
    CAPTURE(offset.x, offset.y);
    edm4eic::TrackerHitCollection hits;
    addHit(hits, 1, helixHit(0, 0.1, 0.01));
    auto middle = helixHit(1, 0.1, 0.01);
    middle.x += offset.x;
    middle.y += offset.y;
    addHit(hits, 2, middle);
    addHit(hits, 3, helixHit(2, 0.1, 0.01));
    CHECK(run(hits).seeds.size() == seeds);
  }
}

TEST_CASE("B0 seeds need the minimum momentum", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  for (std::size_t station = 0; station < 3; ++station) {
    addHit(hits, station + 1, helixHit(station, 2., 0.01));
  }
  CHECK(run(hits).seeds.empty());
  eicrecon::B0TripletSeedingConfig cfg;
  cfg.minMomentum = 0.4 * Acts::UnitConstants::GeV;
  CHECK(run(hits, cfg).seeds.size() == 1);
}

TEST_CASE("B0 seeds of two tracks keep the smallest residuals", "[B0TripletSeeding]") {
  // A second track 50 mm higher, with its middle hit 0.5 mm off the chord along the field
  edm4eic::TrackerHitCollection hits;
  for (std::size_t station = 0; station < 3; ++station) {
    addHit(hits, station + 1, helixHit(station, 0.1, 0.01));
    auto hit = helixHit(station, -0.1, 0.01);
    hit.y += station == 1 ? 50.5F : 50.F;
    addHit(hits, station + 11, hit);
  }
  CHECK(run(hits).seeds.size() == 2);
  eicrecon::B0TripletSeedingConfig cfg;
  cfg.maxSeeds      = 1;
  const auto result = run(hits, cfg);
  REQUIRE(result.seeds.size() == 1);
  CHECK(result.seeds[0].getQuality() > -0.01F);
}

TEST_CASE("B0 seeds use both sensor planes of a station", "[B0TripletSeeding]") {
  // Front and back sensors 6 mm apart on each of three stations give 2^3 triplets
  edm4eic::TrackerHitCollection hits;
  for (std::size_t station = 0; station < 3; ++station) {
    for (const double plane : {0., 6.}) {
      addHit(hits, hits.size() + 1, helixHitAt(station * stationDz + plane, 0.1, 0.01));
    }
  }
  CHECK(run(hits).seeds.size() == 8);
}

TEST_CASE("B0 seeding ignores nonfinite hits", "[B0TripletSeeding]") {
  const float nan = std::numeric_limits<float>::quiet_NaN();
  edm4eic::TrackerHitCollection hits;
  addHit(hits, 9, {0.F, 0.F, nan});
  for (std::size_t station = 0; station < 3; ++station) {
    addHit(hits, station + 1, helixHit(station, 0.1, 0.01));
  }
  addHit(hits, 10, {nan, 0.F, static_cast<float>(stationZ0 + stationDz)});
  const auto result = run(hits);
  REQUIRE(result.seeds.size() == 1);
  for (const auto& hit : result.seeds[0].getHits()) {
    CHECK(hit.getCellID() <= 3);
  }
}

TEST_CASE("B0 seeding rejects a non-positive configuration", "[B0TripletSeeding]") {
  eicrecon::B0TripletSeeding algo("B0TripletSeeding", field());
  eicrecon::B0TripletSeedingConfig cfg;
  cfg.maxResidual = 0.;
  algo.applyConfig(cfg);
  CHECK_THROWS(algo.init());
}

TEST_CASE("B0 seeds skip a station in a quadrupole field", "[B0TripletSeeding]") {
  // A soft track off the magnet midplane, where the field direction turns along the track, on the
  // stations of the B0 tracker without the third one
  constexpr double qOverP = -1. / (1.5 * Acts::UnitConstants::GeV);
  const auto points = quadrupoleHits({20., 60., 5800.}, Acts::Vector3{0.01, 0.01, 1.}.normalized(),
                                     qOverP, {5850., 6272., 6880.});
  edm4eic::TrackerHitCollection hits;
  for (std::size_t i = 0; i < points.size(); ++i) {
    addHit(hits, i + 1, points[i]);
  }
  const auto result = run(hits, {}, std::make_shared<QuadrupoleField>());
  REQUIRE(result.seeds.size() == 1);
  CHECK(result.params[0].getQOverP() == Catch::Approx(qOverP).epsilon(0.05));
}
