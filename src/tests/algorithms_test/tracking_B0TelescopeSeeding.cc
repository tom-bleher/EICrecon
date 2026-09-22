// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#include "algorithms/tracking/B0TelescopeActs.h"
#include <catch2/catch_test_macros.hpp>

#if EICRECON_HAS_B0_TELESCOPE
#include "algorithms/tracking/B0TelescopeMath.h"
#include "algorithms/tracking/B0TelescopeSeedFinder.h"
#include "algorithms/tracking/B0TrackParameterEstimator.h"
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/SpacePointColumns.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <catch2/catch_approx.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace {
using namespace eicrecon::b0telescope;
using Groups = std::array<std::vector<Point>, 4>;

std::vector<std::array<std::uint32_t, 3>> find(const Groups& groups) {
  SpacePoints points(Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR);
  std::array<std::uint32_t, 5> offsets{};
  for (unsigned station = 0; station < 4; ++station) {
    for (const auto& p : groups[station]) {
      auto sp = points.createSpacePoint();
      sp.xy() = {static_cast<float>(p.x), static_cast<float>(p.y)};
      sp.zr() = {static_cast<float>(p.z), static_cast<float>(std::hypot(p.x, p.y))};
    }
    offsets[station + 1] = points.size();
  }
  Seeds seeds;
  seeds.assignSpacePointContainer(points);
  eicrecon::B0TelescopeSeedingConfig cfg;
  cfg.maxSeedsPerMiddle = 100;
  findTelescopeSeeds(points, offsets, 0., cfg, seeds);
  std::vector<std::array<std::uint32_t, 3>> result;
  for (const auto seed : seeds) {
    const auto ids = seed.spacePointIndices();
    REQUIRE(ids.size() == 3);
    result.push_back({ids[0], ids[1], ids[2]});
  }
  return result;
}
}

TEST_CASE("B0 ACTS telescope groups use three distinct stations", "[tracking][b0telescope]") {
  Groups groups;
  for (unsigned s = 0; s < 4; ++s) {
    groups[s].push_back({-150. + 5. * s, 10., 6100. + 250. * s});
  }
  const auto all = find(groups);
  REQUIRE(all.size() == 4);
  CHECK(all == find(groups)); // No state leaks across events.
  for (unsigned missing = 0; missing < 4; ++missing) {
    auto three = groups;
    three[missing].clear();
    CHECK(find(three).size() == 1);
  }
  groups[2].clear();
  groups[3].clear();
  CHECK(find(groups).empty());
  CHECK(find(Groups{}).empty());
}

TEST_CASE("B0 ACTS telescope does not confuse front/back sensors with stations", "[tracking][b0telescope]") {
  Groups groups;
  for (unsigned s = 0; s < 4; ++s) {
    for (double face : {-3.5, 3.5}) {
      groups[s].push_back({-150. + 5. * s + 0.02 * face, 10., 6100. + 250. * s + face});
    }
  }
  const auto result = find(groups);
  REQUIRE(result.size() == 32);
  for (const auto& ids : result) {
    CHECK(ids[0] / 2 < ids[1] / 2);
    CHECK(ids[1] / 2 < ids[2] / 2);
  }
  groups[2].clear();
  groups[3].clear();
  CHECK(find(groups).empty()); // Four hits in two physical stations are not a seed.
}

TEST_CASE("B0 ACTS telescope ignores radial ordering and has no IP cut", "[tracking][b0telescope]") {
  Groups decreasingRadius;
  for (unsigned s = 0; s < 4; ++s) {
    decreasingRadius[s].push_back({-150. + 5. * s, 10., 6100. + 250. * s});
  }
  CHECK(find(decreasingRadius).size() == 4);
  Groups displaced;
  for (unsigned s = 0; s < 4; ++s) {
    displaced[s].push_back({1000. + 5. * s, 2000., 6100. + 250. * s});
  }
  CHECK(find(displaced).size() == 4);
  for (auto& group : displaced) {
    group.front().z = 14000. - group.front().z;
  }
  CHECK(find(displaced).empty()); // Station labels must agree with forward flight.
}

TEST_CASE("B0 ACTS telescope rejects incompatible cross-track combinations", "[tracking][b0telescope]") {
  Groups groups;
  for (unsigned s = 0; s < 4; ++s) {
    groups[s].push_back({-150., -100., 6100. + 250. * s});
    groups[s].push_back({-150., 100., 6100. + 250. * s});
  }
  CHECK(find(groups).size() == 8);
}

TEST_CASE("B0 estimator transports both charge signs using the actual field direction", "[tracking][b0telescope]") {
  using namespace Acts::UnitLiterals;
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const Acts::MagneticFieldContext mctx;
  const auto field = std::make_shared<Acts::ConstantBField>(Acts::Vector3{0., 1.5_T, 0.});
  for (const double charge : {-1., 1.}) {
    // An exact helix in a uniform transverse field; positions are ordered by flight.
    const double radius = 20._GeV / (charge * 1.5_T);
    const double angle = -0.03;
    std::array<Acts::Vector3, 3> positions;
    for (std::size_t i = 0; i < 3; ++i) {
      const double phase = angle - (6100._mm + 250._mm * i) / radius;
      positions[i] = {radius * (std::cos(phase) - std::cos(angle)), 0.,
                       radius * (std::sin(angle) - std::sin(phase))};
    }
    const auto estimate = estimateTrackParameters(positions, 21._ns, gctx, mctx, field, {});
    REQUIRE(estimate.has_value());
    CHECK(estimate->parameters[Acts::eBoundQOverP] * 20._GeV == Catch::Approx(charge).margin(0.01));
    CHECK(estimate->parameters.allFinite());
    CHECK(estimate->covariance.allFinite());
    CHECK(estimate->covariance.diagonal().minCoeff() > 0.);
  }
}

TEST_CASE("B0 estimator rejects unmeasurable or invalid seeds", "[tracking][b0telescope]") {
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const Acts::MagneticFieldContext mctx;
  const auto zero = std::make_shared<Acts::ConstantBField>(Acts::Vector3::Zero());
  std::array<Acts::Vector3, 3> positions{
      Acts::Vector3{-150., 0., 6100.}, Acts::Vector3{-155., 0., 6350.}, Acts::Vector3{-160., 0., 6600.}};
  CHECK_FALSE(estimateTrackParameters(positions, 0., gctx, mctx, zero, {}));
  CHECK_FALSE(estimateTrackParameters(positions, 0., gctx, mctx, nullptr, {}));
  positions[0].x() = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(estimateTrackParameters(positions, 0., gctx, mctx, zero, {}));
}
#else
TEST_CASE("B0 telescope seeding requires the modular ACTS API", "[tracking][b0telescope]") {
  SUCCEED("Unavailable on ACTS < 45.3; requesting telescope mode fails explicitly at startup");
}
#endif
