// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/CovDiag3f.h>
#include <edm4hep/Vector3f.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "algorithms/tracking/B0TripletSeeding.h"

namespace {
using eicrecon::B0TripletSeeding;
using eicrecon::B0TripletSeedingConfig;

void addHit(edm4eic::TrackerHitCollection& hits, std::uint64_t id, double x, double y, double z) {
  auto hit = hits.create();
  hit.setCellID(id);
  hit.setPosition({static_cast<float>(x), static_cast<float>(y), static_cast<float>(z)});
  hit.setPositionError({0.0004f, 0.0004f, 0.0004f});
  hit.setTime(8.f);
}

struct Output {
  edm4eic::TrackSeedCollection seeds;
  edm4eic::TrackParametersCollection parameters;
};

void seed(const edm4eic::TrackerHitCollection& hits, const B0TripletSeedingConfig& cfg,
          Output& output, double fieldTesla = 1.) {
#if Acts_VERSION_MAJOR >= 45
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  const Acts::GeometryContext gctx;
#endif
  const Acts::MagneticFieldContext mctx;
  const Acts::ConstantBField field({0., fieldTesla * Acts::UnitConstants::T, 0.});
  B0TripletSeeding::seedHits(cfg, gctx, mctx, field, hits, output.seeds, output.parameters);
}

using Triplet = std::array<std::uint64_t, 3>;
std::vector<Triplet> triplets(const Output& output) {
  std::vector<Triplet> result;
  for (const auto& seed : output.seeds) {
    REQUIRE(seed.getHits().size() == 3);
    result.push_back({seed.getHits()[0].getCellID(), seed.getHits()[1].getCellID(),
                      seed.getHits()[2].getCellID()});
  }
  return result;
}

void checkParameters(const Output& output) {
  REQUIRE(output.seeds.size() == output.parameters.size());
  for (std::size_t i = 0; i < output.parameters.size(); ++i) {
    const auto parameters = output.parameters[i];
    CHECK(std::isfinite(parameters.getPhi()));
    CHECK(std::isfinite(parameters.getTheta()));
    CHECK(std::isfinite(parameters.getQOverP()));
    CHECK(std::isfinite(parameters.getTime()));
    CHECK(parameters.getTheta() > 0.f);
    CHECK(parameters.getTheta() < 0.5f);
    const auto covariance = parameters.getCovariance();
    for (std::size_t j = 0; j < 6; ++j) {
      CHECK(std::isfinite(covariance(j, j)));
      CHECK(covariance(j, j) > 0.f);
    }
    CHECK(output.seeds[i].getParams().getObjectID() == parameters.getObjectID());
  }
}

void addStraightTracks(edm4eic::TrackerHitCollection& hits, bool reverse = false) {
  struct Hit {
    std::uint64_t id;
    double x;
    double z;
  };
  std::vector<Hit> points;
  for (std::size_t station = 0; station < 4; ++station) {
    const double z = 5900. + 300. * station;
    for (std::size_t track = 0; track < 10; ++track) {
      points.push_back({100 * station + track + 1, 0.02 * z + track, z});
    }
  }
  if (reverse) {
    std::ranges::reverse(points);
  }
  for (const auto& point : points) {
    addHit(hits, point.id, point.x, 0., point.z);
  }
}
} // namespace

TEST_CASE("B0 triplet seeding estimates both charge signs in a transverse field",
          "[B0TripletSeeding]") {
  // Exact circular trajectories: q>0 bends towards -x in a +y field. The initial
  // direction is tilted so the perigee phi and longitudinal prior are defined.
  for (const double charge : {-1., 1.}) {
    CAPTURE(charge);
    constexpr double momentum     = 10.;
    const double radius           = charge * 1000. * momentum / 0.299792458;
    constexpr double initialAngle = 0.04;
    edm4eic::TrackerHitCollection hits;
    const std::array<double, 3> zPositions{5850., 6272.2, 6601.};
    for (std::size_t i = 0; i < zPositions.size(); ++i) {
      const double angle =
          std::asin(std::sin(initialAngle) + (zPositions[i] - zPositions[0]) / radius);
      const double x = -200. + radius * (std::cos(angle) - std::cos(initialAngle));
      addHit(hits, i + 1, x, 0., zPositions[i]);
    }
    Output output;
    seed(hits, {}, output);
    REQUIRE(output.seeds.size() == 1);
    checkParameters(output);
    CHECK(output.parameters[0].getQOverP() == Catch::Approx(charge / momentum).epsilon(0.005));
    CHECK(output.seeds[0].getPerigee().z < hits[0].getPosition().z);
  }
}

TEST_CASE("B0 straight triplets retain finite positive curvature uncertainty",
          "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  for (std::size_t station = 0; station < 3; ++station) {
    const double z = 5900. + 300. * station;
    addHit(hits, station + 1, 0.02 * z, 0., z);
  }
  Output output;
  seed(hits, {}, output);
  REQUIRE(output.seeds.size() == 1);
  checkParameters(output);
  CHECK(std::abs(output.parameters[0].getQOverP()) < 1.e-5f);
  CHECK(output.parameters[0].getCovariance()(4, 4) > 0.f);
}

TEST_CASE("B0 triplets need three separated stations and a nonzero field", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  addHit(hits, 1, -200., 0., 5900.);
  addHit(hits, 2, -199., 0., 5905.);
  addHit(hits, 3, -190., 0., 6200.);
  Output output;
  seed(hits, {}, output);
  CHECK(output.seeds.empty());
  CHECK(output.parameters.empty());
  addHit(hits, 4, -175., 0., 6500.);
  seed(hits, {}, output, 0.);
  CHECK(output.seeds.empty());
  CHECK(output.parameters.empty());
}

TEST_CASE("B0 capped seeds are independent of hit insertion order", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  edm4eic::TrackerHitCollection reversed;
  addStraightTracks(hits);
  addStraightTracks(reversed, true);
  B0TripletSeedingConfig cfg;
  cfg.maxSeeds = 7;
  Output first;
  Output second;
  seed(hits, cfg, first, 1.7);
  seed(reversed, cfg, second, 1.7);
  REQUIRE(first.seeds.size() == cfg.maxSeeds);
  REQUIRE(second.seeds.size() == cfg.maxSeeds);
  CHECK(triplets(first) == triplets(second));
  checkParameters(first);
  for (std::size_t i = 0; i < first.parameters.size(); ++i) {
    CHECK(first.parameters[i].getQOverP() == second.parameters[i].getQOverP());
    CHECK(first.parameters[i].getPhi() == second.parameters[i].getPhi());
    CHECK(first.parameters[i].getTheta() == second.parameters[i].getTheta());
  }
  // Repeated processing into populated outputs must respect the event cap.
  seed(hits, cfg, first, 1.7);
  CHECK(first.seeds.size() == cfg.maxSeeds);
  CHECK(first.parameters.size() == cfg.maxSeeds);
}

TEST_CASE("B0 invalid leading candidates do not consume the seed cap", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  addStraightTracks(hits);
  B0TripletSeedingConfig cfg;
  cfg.minMomentum = 1000.;
  cfg.maxSeeds    = 5000;
  Output reference;
  seed(hits, cfg, reference, 1.7);
  REQUIRE(reference.seeds.size() > 7);
  REQUIRE(reference.seeds.size() < 4000);
  cfg.maxSeeds = 7;
  Output capped;
  seed(hits, cfg, capped, 1.7);
  REQUIRE(capped.seeds.size() == cfg.maxSeeds);
  auto expected = triplets(reference);
  expected.resize(cfg.maxSeeds);
  CHECK(triplets(capped) == expected);
  checkParameters(capped);
}

TEST_CASE("B0 seeding rejects invalid configuration and tolerates invalid hits",
          "[B0TripletSeeding]") {
  B0TripletSeeding algorithm("test_b0_invalid_configuration");
  B0TripletSeedingConfig cfg;
  cfg.stationGap = 0.;
  algorithm.applyConfig(cfg);
  CHECK_THROWS(algorithm.init());
  cfg          = {};
  cfg.maxSeeds = 0;
  algorithm.applyConfig(cfg);
  CHECK_THROWS(algorithm.init());
  cfg             = {};
  cfg.maxResidual = std::numeric_limits<double>::quiet_NaN();
  algorithm.applyConfig(cfg);
  CHECK_THROWS(algorithm.init());
  cfg                = {};
  cfg.anchorDistance = 0.;
  algorithm.applyConfig(cfg);
  CHECK_NOTHROW(algorithm.init());

  edm4eic::TrackerHitCollection hits;
  addHit(hits, 1, std::numeric_limits<double>::quiet_NaN(), 0., 5900.);
  addHit(hits, 2, -190., 0., 6200.);
  addHit(hits, 3, -175., 0., 6500.);
  Output output;
  seed(hits, {}, output);
  CHECK(output.seeds.empty());
  CHECK(output.parameters.empty());
}
