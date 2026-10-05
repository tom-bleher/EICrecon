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
#include <random>
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

B0TripletSeeding::Stats seed(const edm4eic::TrackerHitCollection& hits,
                             const B0TripletSeedingConfig& cfg, Output& output,
                             double fieldTesla = 1.) {
#if Acts_VERSION_MAJOR >= 45
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  const Acts::GeometryContext gctx;
#endif
  const Acts::MagneticFieldContext mctx;
  const Acts::ConstantBField field({0., fieldTesla * Acts::UnitConstants::T, 0.});
  return B0TripletSeeding::seedHits(cfg, gctx, mctx, field, hits, output.seeds, output.parameters);
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
  cfg            = {};
  cfg.stationGap = std::numeric_limits<double>::infinity();
  algorithm.applyConfig(cfg);
  CHECK_THROWS(algorithm.init());
  cfg                     = {};
  cfg.qOverPRelativeError = std::numeric_limits<double>::infinity();
  algorithm.applyConfig(cfg);
  CHECK_THROWS(algorithm.init());
  cfg                = {};
  cfg.anchorDistance = 0.;
  algorithm.applyConfig(cfg);
  CHECK_NOTHROW(algorithm.init());

  // Invalid hits are dropped: a maximum-based variance check would keep the
  // negative and NaN covariance components below.
  auto addValid = [](edm4eic::TrackerHitCollection& hits) {
    for (std::size_t station = 0; station < 3; ++station) {
      const double z = 5900. + 300. * station;
      for (std::size_t track = 0; track < 2; ++track) {
        addHit(hits, 100 * station + track + 1, 0.02 * z + track, 0., z);
      }
    }
  };
  edm4eic::TrackerHitCollection valid;
  addValid(valid);
  Output reference;
  seed(valid, {}, reference);
  REQUIRE(!reference.seeds.empty());
  edm4eic::TrackerHitCollection polluted;
  addValid(polluted);
  auto nanPos = polluted.create();
  nanPos.setCellID(101);
  nanPos.setPosition({std::numeric_limits<float>::quiet_NaN(), 0.f, 5900.f});
  nanPos.setPositionError({0.0004f, 0.0004f, 0.0004f});
  nanPos.setTime(8.f);
  auto negVar = polluted.create();
  negVar.setCellID(102);
  negVar.setPosition({124.f, 0.f, 6200.f});
  negVar.setPositionError({0.0004f, -1.f, 0.0004f});
  negVar.setTime(8.f);
  auto nanVar = polluted.create();
  nanVar.setCellID(103);
  nanVar.setPosition({130.f, 0.f, 6500.f});
  nanVar.setPositionError({0.0004f, 0.0004f, std::numeric_limits<float>::quiet_NaN()});
  nanVar.setTime(8.f);
  Output output;
  seed(polluted, {}, output);
  CHECK(triplets(output) == triplets(reference));
  CHECK(output.parameters.size() == reference.parameters.size());
}

TEST_CASE("B0 seeding reports exact diagnostics counters", "[B0TripletSeeding]") {
  edm4eic::TrackerHitCollection hits;
  for (std::size_t station = 0; station < 3; ++station) {
    const double z = 5900. + 300. * station;
    for (std::size_t track = 0; track < 2; ++track) {
      addHit(hits, 100 * station + track + 1, 0.02 * z + track, 0., z);
    }
  }
  Output output;
  const auto stats = seed(hits, {}, output);
  CHECK(stats.stations == 3);
  CHECK(stats.enumerated == 8);
  CHECK(stats.residualPassing == 8);
  CHECK(stats.rankPruned == 0);
  CHECK(stats.estimated == 8);
  CHECK(stats.heapReplacements == 0);
  CHECK(stats.retained + stats.estimationFailures == 8);
  CHECK(stats.retained == output.seeds.size());

  B0TripletSeedingConfig capped;
  capped.maxSeeds = 3;
  Output saturated;
  const auto saturatedStats = seed(hits, capped, saturated);
  CHECK(saturatedStats.retained == 3);
  CHECK(saturatedStats.estimated + saturatedStats.rankPruned == saturatedStats.residualPassing);
  CHECK(saturatedStats.retained ==
        std::min<std::size_t>(3, saturatedStats.estimated - saturatedStats.estimationFailures));
  CHECK(saturatedStats.rankPruned > 0);
}

TEST_CASE("B0 capped seeds match the exhaustive reference under saturation", "[B0TripletSeeding]") {
  struct Hit {
    std::uint64_t id;
    double x;
    double y;
    double z;
  };
  std::vector<Hit> points;
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> near(-0.3, 0.3);
  std::uniform_real_distribution<double> jitter(-0.5, 0.5);
  std::uniform_real_distribution<double> wild(-50., 50.);
  std::uint64_t id = 1;
  for (std::size_t station = 0; station < 4; ++station) {
    const double z = 5900. + 300. * station;
    for (int k = 0; k < 3; ++k) {
      points.push_back({id++, 0.02 * z + jitter(rng), near(rng), z});
    }
    for (int k = 0; k < 3; ++k) {
      points.push_back({id++, wild(rng), wild(rng), z});
    }
  }
  auto build = [](const std::vector<Hit>& order) {
    edm4eic::TrackerHitCollection hits;
    for (const auto& point : order) {
      addHit(hits, point.id, point.x, point.y, point.z);
    }
    return hits;
  };
  B0TripletSeedingConfig full;
  full.maxSeeds = 100000;
  Output reference;
  const auto referenceStats = seed(build(points), full, reference);
  REQUIRE(reference.seeds.size() > 7);
  REQUIRE(reference.seeds.size() < 864);
  CHECK(referenceStats.rankPruned == 0);
  for (std::size_t i = 1; i < reference.seeds.size(); ++i) {
    CHECK(reference.seeds[i - 1].getQuality() >= reference.seeds[i].getQuality());
  }
  B0TripletSeedingConfig capped;
  capped.maxSeeds                            = 7;
  const auto expected                        = triplets(reference);
  std::vector<std::vector<Hit>> permutations = {points, points, points};
  std::shuffle(permutations[1].begin(), permutations[1].end(), rng);
  std::reverse(permutations[2].begin(), permutations[2].end());
  for (const auto& order : permutations) {
    Output output;
    const auto stats = seed(build(order), capped, output);
    CHECK(stats.retained == 7);
    CHECK(stats.estimated + stats.rankPruned == stats.residualPassing);
    CHECK(stats.heapReplacements > 0);
    const auto selected = triplets(output);
    REQUIRE(selected.size() == 7);
    for (std::size_t i = 0; i < 7; ++i) {
      CHECK(selected[i] == expected[i]);
    }
  }
}
