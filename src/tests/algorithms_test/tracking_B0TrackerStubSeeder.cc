// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

#include "algorithms/tracking/B0TrackerStubSeeder.h"

using Catch::Approx;
using eicrecon::b0stub::assignStation;
using eicrecon::b0stub::clusterStations;
using eicrecon::b0stub::covarianceModelAdditions;
using eicrecon::b0stub::endpointCompatibility;
using eicrecon::b0stub::FieldIntegralFit;
using eicrecon::b0stub::FieldSample;
using eicrecon::b0stub::fitFieldIntegral;
using eicrecon::b0stub::fitStub;
using eicrecon::b0stub::groupHitsByStation;
using eicrecon::b0stub::integrateFieldSamples;
using eicrecon::b0stub::perigeeFromRay;
using eicrecon::b0stub::Point3;
using eicrecon::b0stub::seedCovarianceFromFit;
using eicrecon::b0stub::StubFit;

namespace {

/// Ion-frame z of the four official B0 disks [mm].
constexpr std::array<double, 4> kOfficialStationZ{5902.0, 6172.0, 6442.0, 6712.0};
/// Realistic-geometry station offsets 0 / 269.93 / 469.54 / 809.93 mm, anchored
/// at the official first-disk z so the lever arm stays 810 mm.
constexpr std::array<double, 4> kRealisticStationZ{5902.0, 6171.93, 6371.54, 6711.93};

constexpr std::array<std::array<double, 4>, 2> kGeometries{kOfficialStationZ, kRealisticStationZ};

/// Sample a parabola/line trajectory at the given stations.
std::vector<Point3> sampleTrack(double c0, double c1, double c2, double b0, double b1,
                                const std::array<double, 4>& stations,
                                double varianceX = 0.0, double varianceY = 0.0) {
  std::vector<Point3> pts;
  for (double z : stations) {
    pts.push_back({c0 + c1 * z + c2 * z * z, b0 + b1 * z, z, varianceX, varianceY});
  }
  return pts;
}

std::vector<Point3> sampleTrack(double c0, double c1, double c2, double b0, double b1,
                                double varianceX = 0.0, double varianceY = 0.0) {
  return sampleTrack(c0, c1, c2, b0, b1, kOfficialStationZ, varianceX, varianceY);
}

} // namespace

TEST_CASE("B0 stub fit recovers the generating trajectory", "[B0TrackerStubSeeder]") {
  const double c0 = 12.5;
  const double c1 = -3.0e-3;
  const double c2 = -4.3878e-6;
  const double b0 = -1.5;
  const double b1 = 2.0e-3;

  for (const auto& stations : kGeometries) {
    const StubFit fit = fitStub(sampleTrack(c0, c1, c2, b0, b1, stations));

    REQUIRE(fit.valid);
    // The fit is over-determined (4 points, 3 parameters) but exact for noiseless
    // input, so the residuals must vanish and the coefficients must come back.
    CHECK(fit.rmsX == Approx(0.0).margin(1e-6));
    CHECK(fit.rmsY == Approx(0.0).margin(1e-9));
    CHECK(fit.c2 == Approx(c2).epsilon(1e-6));
    CHECK(fit.b1 == Approx(b1).epsilon(1e-9));
    // Evaluate rather than compare c0/c1 directly: they are the values of a
    // polynomial extrapolated back to z = 0, six metres outside the fit range.
    for (double z : stations) {
      CHECK(fit.x(z) == Approx(c0 + c1 * z + c2 * z * z).margin(1e-6));
      CHECK(fit.y(z) == Approx(b0 + b1 * z).margin(1e-9));
    }
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

TEST_CASE("B0 stub fit rejects a rank-deficient quadratic", "[B0TrackerStubSeeder]") {
  std::vector<Point3> pts{
      {0.0, 0.0, 5902.0}, {1.0, 0.0, 5902.0}, {2.0, 0.0, 6172.0}, {3.0, 0.0, 6172.0}};
  CHECK_FALSE(fitStub(pts).valid);
}

TEST_CASE("B0 weighted fit derives coefficient covariance", "[B0TrackerStubSeeder]") {
  constexpr double varianceX = 4.0e-4;
  constexpr double varianceY = 9.0e-4;
  const StubFit fit =
      fitStub(sampleTrack(12.5, -3.0e-3, -4.0e-6, -1.5, 2.0e-3, varianceX, varianceY));
  REQUIRE(fit.valid);
  REQUIRE(fit.covarianceValid);
  for (int i = 0; i < 3; ++i) {
    CHECK(fit.covarianceX[3 * i + i] > 0.0);
  }
  for (int i = 0; i < 2; ++i) {
    CHECK(fit.covarianceY[2 * i + i] > 0.0);
  }

  const StubFit scaled =
      fitStub(sampleTrack(12.5, -3.0e-3, -4.0e-6, -1.5, 2.0e-3, 4.0 * varianceX, 4.0 * varianceY));
  REQUIRE(scaled.covarianceValid);
  CHECK(scaled.covarianceX[8] == Approx(4.0 * fit.covarianceX[8]).epsilon(1e-10));
  CHECK(scaled.covarianceY[3] == Approx(4.0 * fit.covarianceY[3]).epsilon(1e-10));
}

TEST_CASE("B0 three-station fit retains measurement covariance", "[B0TrackerStubSeeder]") {
  auto points = sampleTrack(12.5, -3.0e-3, -4.0e-6, -1.5, 2.0e-3, 4.0e-4, 4.0e-4);
  points.pop_back();
  const StubFit fit = fitStub(points);
  REQUIRE(fit.valid);
  REQUIRE(fit.covarianceValid);
  CHECK(fit.covarianceX[8] > 0.0);
  CHECK(fit.rmsX == Approx(0.0).margin(1e-6));
}

TEST_CASE("B0 field quadrature is exact for a linear field", "[B0TrackerStubSeeder]") {
  constexpr double z0       = 5800.0;
  constexpr double field0   = 1.1;
  constexpr double gradient = 2.0e-4;
  std::vector<FieldSample> samples;
  for (double z : {z0, 5902.0, 6172.0, 6712.0}) {
    samples.push_back({z, field0 + gradient * (z - z0)});
  }
  const auto moments = integrateFieldSamples(samples);
  REQUIRE(moments.size() == samples.size());
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const double dz = samples[i].z - z0;
    CHECK(moments[i].first == Approx(field0 * dz + 0.5 * gradient * dz * dz).epsilon(1e-12));
    CHECK(moments[i].second ==
          Approx(0.5 * field0 * dz * dz + gradient * dz * dz * dz / 6.0).epsilon(1e-12));
  }
}

TEST_CASE("B0 field-integral fit recovers signed q over p", "[B0TrackerStubSeeder]") {
  constexpr double zEntrance  = 5800.0;
  constexpr double fieldY     = 1.184;
  constexpr double qOverP     = 1.0 / 41.0;
  constexpr double xEntrance  = -2.5;
  constexpr double txEntrance = -1.0e-3;
  constexpr double kBend      = 2.998e-4;
  for (const auto& stations : kGeometries) {
    std::vector<Point3> points;
    std::vector<double> integrals;
    for (double z : stations) {
      const double dz       = z - zEntrance;
      const double integral = 0.5 * fieldY * dz * dz;
      points.push_back({xEntrance + txEntrance * dz - kBend * qOverP * integral, -1.5 + 2.0e-3 * z,
                        z, 4.0e-4, 9.0e-4});
      integrals.push_back(integral);
    }
    const FieldIntegralFit fit = fitFieldIntegral(points, integrals, zEntrance);
    REQUIRE(fit.valid);
    REQUIRE(fit.covarianceValid);
    CHECK(fit.xReference == Approx(xEntrance).margin(1e-10));
    CHECK(fit.txReference == Approx(txEntrance).margin(1e-12));
    CHECK(fit.qOverP == Approx(qOverP).epsilon(1e-10));
    CHECK(fit.rmsX == Approx(0.0).margin(1e-10));

    auto scaledPoints = points;
    for (auto& point : scaledPoints) {
      point.varianceX *= 4.0;
    }
    const auto scaled = fitFieldIntegral(scaledPoints, integrals, zEntrance);
    REQUIRE(scaled.covarianceValid);
    CHECK(scaled.covariance[8] == Approx(4.0 * fit.covariance[8]).epsilon(1e-10));
  }
}

TEST_CASE("B0 endpoint pruning keeps prompt pairs and rejects cross-pairs",
          "[B0TrackerStubSeeder]") {
  const Point3 promptFirst{0.0, 12.0, 5902.0};
  const Point3 promptLast{0.0, 12.0 * 6712.0 / 5902.0, 6712.0};
  const auto prompt = endpointCompatibility(promptFirst, promptLast, 0.10, 5.0, true);
  REQUIRE(prompt.valid);
  CHECK(prompt.beamResidual == Approx(0.0).margin(1e-12));

  // This pair has a modest slope but extrapolates far from the beamline.
  const auto cross =
      endpointCompatibility({0.0, 100.0, 5902.0}, {0.0, 105.0, 6712.0}, 0.10, 5.0, true);
  CHECK_FALSE(cross.valid);
  const auto unconstrained =
      endpointCompatibility({0.0, 100.0, 5902.0}, {0.0, 105.0, 6712.0}, 0.10, 5.0, false);
  CHECK(unconstrained.valid);

  CHECK_FALSE(
      endpointCompatibility({0.0, 0.0, 5902.0}, {0.0, 100.0, 6000.0}, 0.10, 5.0, true).valid);
  CHECK_FALSE(endpointCompatibility(promptLast, promptFirst, 0.10, 5.0, true).valid);
}

TEST_CASE("B0 calibrated covariance additions include material and field terms",
          "[B0TrackerStubSeeder]") {
  constexpr double qOverP = 0.04;
  const auto additions =
      covarianceModelAdditions(qOverP, 1.0e-5, 0.15, 2.0e-9, 1.0e-3, 4.0e-8, 0.02, 0.01);
  CHECK(additions[0] == Approx(1.0e-5 + std::pow(0.15 * qOverP, 2)).epsilon(1e-12));
  CHECK(additions[1] == Approx(2.0e-9 + std::pow(1.0e-3 * qOverP, 2)).epsilon(1e-12));
  CHECK(additions[2] ==
        Approx(4.0e-8 + (0.02 * 0.02 + 0.01 * 0.01) * qOverP * qOverP).epsilon(1e-12));

  const auto opposite =
      covarianceModelAdditions(-qOverP, 1.0e-5, 0.15, 2.0e-9, 1.0e-3, 4.0e-8, 0.02, 0.01);
  CHECK(opposite == additions);
}

TEST_CASE("B0 field-fit covariance propagates to correlated seed parameters",
          "[B0TrackerStubSeeder]") {
  constexpr double zEntrance = 5800.0;
  constexpr double fieldY    = 1.184;
  constexpr double qOverP    = 1.0 / 41.0;
  constexpr double kBend     = 2.998e-4;
  std::vector<Point3> points;
  std::vector<double> integrals;
  for (double z : kOfficialStationZ) {
    const double dz       = z - zEntrance;
    const double integral = 0.5 * fieldY * dz * dz;
    points.push_back(
        {-2.5 - 1.0e-3 * dz - kBend * qOverP * integral, -1.5 + 2.0e-3 * z, z, 4.0e-4, 4.0e-4});
    integrals.push_back(integral);
  }
  const StubFit nonBendFit       = fitStub(points);
  const FieldIntegralFit bendFit = fitFieldIntegral(points, integrals, zEntrance);
  REQUIRE(nonBendFit.covarianceValid);
  REQUIRE(bendFit.covarianceValid);
  const auto flat = seedCovarianceFromFit(bendFit, nonBendFit, -0.025, true);

  Eigen::Matrix<double, 5, 5> covariance;
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      covariance(row, col) = flat[5 * row + col];
      CHECK(std::isfinite(flat[5 * row + col]));
      CHECK(flat[5 * row + col] == Approx(flat[5 * col + row]).margin(1e-15));
    }
  }
  CHECK(covariance(2, 2) > 0.0);
  CHECK(covariance(3, 3) > 0.0);
  CHECK(covariance(4, 4) > 0.0);
  CHECK(std::abs(covariance(3, 4)) > 0.0);
  CHECK(covariance(4, 4) == Approx(bendFit.covariance[8]).epsilon(1e-8));
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 5, 5>> solver(covariance);
  REQUIRE(solver.info() == Eigen::Success);
  CHECK(solver.eigenvalues().minCoeff() >= Approx(0.0).margin(1e-12));
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

TEST_CASE("B0 station clustering keeps official and realistic disks split",
          "[B0TrackerStubSeeder]") {
  constexpr double gap = 50.0;
  for (const auto& stations : kGeometries) {
    const auto clustered = clusterStations({stations.begin(), stations.end()}, gap);
    REQUIRE(clustered.size() == 4);
    for (std::size_t i = 0; i < stations.size(); ++i) {
      CHECK(clustered[i].zMean == Approx(stations[i]).margin(1e-9));
      CHECK(assignStation(stations[i], clustered, gap) == static_cast<int>(i));
    }
  }
  CHECK(kRealisticStationZ[1] - kRealisticStationZ[0] == Approx(269.93).margin(1e-9));
  CHECK(kRealisticStationZ[2] - kRealisticStationZ[1] == Approx(199.61).margin(1e-9));
  CHECK(kRealisticStationZ[3] - kRealisticStationZ[2] == Approx(340.39).margin(1e-9));
}

TEST_CASE("B0 station clustering merges realistic front/back faces", "[B0TrackerStubSeeder]") {
  constexpr double gap = 50.0;
  std::vector<double> surfaceZ;
  for (double z : kRealisticStationZ) {
    surfaceZ.push_back(z - 3.5);
    surfaceZ.push_back(z + 3.5);
  }
  const auto clustered = clusterStations(surfaceZ, gap);
  REQUIRE(clustered.size() == 4);
  for (std::size_t i = 0; i < kRealisticStationZ.size(); ++i) {
    CHECK(clustered[i].zMean == Approx(kRealisticStationZ[i]).margin(1e-9));
    CHECK(clustered[i].zMax - clustered[i].zMin == Approx(7.0).margin(1e-9));
  }
}

TEST_CASE("B0 cached stations drop mid-gap hits instead of inventing a station",
          "[B0TrackerStubSeeder]") {
  constexpr double gap     = 50.0;
  const auto stations      = clusterStations({kRealisticStationZ.begin(), kRealisticStationZ.end()}, gap);
  const double midGap      = 0.5 * (kRealisticStationZ[0] + kRealisticStationZ[1]);
  std::vector<double> hitZ = {kRealisticStationZ[0], midGap, kRealisticStationZ[1],
                              kRealisticStationZ[2], kRealisticStationZ[3]};
  const auto cached        = groupHitsByStation(hitZ, stations, gap);
  REQUIRE(cached.size() == 4);
  CHECK(cached.at(0).size() == 1);
  CHECK(cached.at(1).size() == 1);
  CHECK(assignStation(midGap, stations, gap) == -1);

  const auto eventContent = groupHitsByStation(hitZ, {}, gap);
  CHECK(eventContent.size() == 5);
}
