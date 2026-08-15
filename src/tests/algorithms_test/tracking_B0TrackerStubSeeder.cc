// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>

#include "algorithms/tracking/B0TrackerStubSeeder.h"

using Catch::Approx;
using eicrecon::b0stub::chargeFromCurvature;
using eicrecon::b0stub::fitStub;
using eicrecon::b0stub::momentumFromCurvature;
using eicrecon::b0stub::perigeeFromRay;
using eicrecon::b0stub::Point3;
using eicrecon::b0stub::StubFit;

namespace {

/// Ion-frame z of the four official B0 disks [mm].
constexpr std::array<double, 4> kStationZ{5902.0, 6172.0, 6442.0, 6712.0};

/// Sample a parabola/line trajectory at the B0 stations.
std::vector<Point3> sampleTrack(double c0, double c1, double c2, double b0, double b1) {
  std::vector<Point3> pts;
  for (double z : kStationZ) {
    pts.push_back({c0 + c1 * z + c2 * z * z, b0 + b1 * z, z});
  }
  return pts;
}

} // namespace

TEST_CASE("B0 stub fit recovers the generating trajectory", "[B0TrackerStubSeeder]") {
  const double c0 = 12.5;
  const double c1 = -3.0e-3;
  const double c2 = -4.3878e-6;
  const double b0 = -1.5;
  const double b1 = 2.0e-3;

  const StubFit fit = fitStub(sampleTrack(c0, c1, c2, b0, b1));

  REQUIRE(fit.valid);
  // The fit is over-determined (4 points, 3 parameters) but exact for noiseless
  // input, so the residuals must vanish and the coefficients must come back.
  CHECK(fit.rmsX == Approx(0.0).margin(1e-6));
  CHECK(fit.rmsY == Approx(0.0).margin(1e-9));
  CHECK(fit.c2 == Approx(c2).epsilon(1e-6));
  CHECK(fit.b1 == Approx(b1).epsilon(1e-9));
  // Evaluate rather than compare c0/c1 directly: they are the values of a
  // polynomial extrapolated back to z = 0, six metres outside the fit range.
  for (double z : kStationZ) {
    CHECK(fit.x(z) == Approx(c0 + c1 * z + c2 * z * z).margin(1e-6));
    CHECK(fit.y(z) == Approx(b0 + b1 * z).margin(1e-9));
  }
}

TEST_CASE("B0 stub fit needs at least three points", "[B0TrackerStubSeeder]") {
  std::vector<Point3> pts{{0.0, 0.0, 5902.0}, {1.0, 0.0, 6172.0}};
  CHECK_FALSE(fitStub(pts).valid);
}

TEST_CASE("B0 stub fit rejects degenerate z", "[B0TrackerStubSeeder]") {
  std::vector<Point3> pts{{0.0, 0.0, 6000.0}, {1.0, 0.0, 6000.0}, {2.0, 0.0, 6000.0}};
  CHECK_FALSE(fitStub(pts).valid);
}

TEST_CASE("B0 sagitta returns the generating momentum and charge", "[B0TrackerStubSeeder]") {
  const double By = 1.184; // B0pf nominal [T]

  // d(tx)/dz = -q * 2.998e-4 * By / p, so c2 = 0.5 * that.
  auto curvature = [&](double p, int q) { return -0.5 * q * 2.998e-4 * By / p; };

  SECTION("41 GeV proton") {
    const double p    = 41.0;
    const double c2   = curvature(p, +1);
    const StubFit fit = fitStub(sampleTrack(12.5, -3.0e-3, c2, 0.0, 0.0));
    REQUIRE(fit.valid);
    CHECK(momentumFromCurvature(fit.c2, By) == Approx(p).epsilon(1e-4));
    CHECK(chargeFromCurvature(fit.c2, By) == +1);
  }

  SECTION("5 GeV pi-") {
    const double p    = 5.0;
    const double c2   = curvature(p, -1);
    const StubFit fit = fitStub(sampleTrack(12.5, -3.0e-3, c2, 0.0, 0.0));
    REQUIRE(fit.valid);
    CHECK(momentumFromCurvature(fit.c2, By) == Approx(p).epsilon(1e-4));
    CHECK(chargeFromCurvature(fit.c2, By) == -1);
  }

  SECTION("charge inference follows the field sign") {
    const double c2 = curvature(41.0, +1);
    CHECK(chargeFromCurvature(c2, By) == +1);
    CHECK(chargeFromCurvature(c2, -By) == -1);
  }

  SECTION("a straight track has no measurable momentum") {
    CHECK(momentumFromCurvature(0.0, By) == 0.0);
  }
}

TEST_CASE("B0 perigee parameters of a ray through the origin", "[B0TrackerStubSeeder]") {
  // A ray leaving the origin at 25 mrad in the x-z plane: it passes through the
  // perigee line, so both local coordinates must vanish.
  const double slope = 0.025;
  const Point3 dir{slope, 0.0, 1.0};
  const auto par = perigeeFromRay({0.0, 0.0, 0.0}, dir, {0.0, 0.0, 0.0});

  CHECK(par.loc0 == Approx(0.0).margin(1e-9));
  CHECK(par.loc1 == Approx(0.0).margin(1e-9));
  CHECK(par.phi == Approx(0.0).margin(1e-9));
  CHECK(par.theta == Approx(std::atan(slope)).epsilon(1e-9));
}

TEST_CASE("B0 perigee parameters of a displaced ray", "[B0TrackerStubSeeder]") {
  // Parallel to z, offset by +3 mm in y: the impact parameter is 3 mm and the
  // track is at phi = pi/2 relative to the point of closest approach.
  const auto par = perigeeFromRay({0.0, 3.0, 1000.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
  CHECK(std::abs(par.loc0) == Approx(3.0).epsilon(1e-9));
  CHECK(par.theta == Approx(0.0).margin(1e-9));
}

TEST_CASE("B0 perigee follows the surface it is expressed on", "[B0TrackerStubSeeder]") {
  // The same ray, expressed on a perigee surface at the field entrance rather
  // than at the origin: loc1 is measured from the surface centre.
  const Point3 dir{0.025, 0.0, 1.0};
  const auto atOrigin   = perigeeFromRay({0.0, 0.0, 0.0}, dir, {0.0, 0.0, 0.0});
  const auto atEntrance = perigeeFromRay({0.0, 0.0, 0.0}, dir, {0.0, 0.0, 5800.0});

  CHECK(atOrigin.theta == Approx(atEntrance.theta).epsilon(1e-12));
  CHECK(atOrigin.phi == Approx(atEntrance.phi).epsilon(1e-12));
  CHECK(atOrigin.loc1 == Approx(0.0).margin(1e-9));
  CHECK(atEntrance.loc1 == Approx(-5800.0).epsilon(1e-9));
}
