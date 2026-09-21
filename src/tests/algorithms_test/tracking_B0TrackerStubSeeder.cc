// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "algorithms/tracking/B0TrackerStubSeeder.h"

using Catch::Approx;
using eicrecon::b0stub::assignStation;
using eicrecon::b0stub::beamSpotCovarianceAdditions;
using eicrecon::b0stub::BendFit;
using eicrecon::b0stub::bendFitFromParabola;
using eicrecon::b0stub::clusterStations;
using eicrecon::b0stub::endpointCompatibility;
using eicrecon::b0stub::fitStub;
using eicrecon::b0stub::fitWithNonBendCorrection;
using eicrecon::b0stub::groupHitsByStation;
using eicrecon::b0stub::perigeeFromRay;
using eicrecon::b0stub::Point3;
using eicrecon::b0stub::removeNonBendCurvature;
using eicrecon::b0stub::scatteringCovarianceAdditions;
using eicrecon::b0stub::seedCovarianceFromFit;
using eicrecon::b0stub::seedCovarianceWithFallbacks;
using eicrecon::b0stub::seedParametersFromFit;
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
                                const std::array<double, 4>& stations, double varianceX = 0.0,
                                double varianceY = 0.0) {
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
  std::vector<Point3> pts{{.x = 0.0, .y = 0.0, .z = 5902.0}, {.x = 1.0, .y = 0.0, .z = 6172.0}};
  CHECK_FALSE(fitStub(pts).valid);
}

TEST_CASE("B0 stub fit rejects degenerate z", "[B0TrackerStubSeeder]") {
  std::vector<Point3> pts{{.x = 0.0, .y = 0.0, .z = 6000.0},
                          {.x = 1.0, .y = 0.0, .z = 6000.0},
                          {.x = 2.0, .y = 0.0, .z = 6000.0}};
  CHECK_FALSE(fitStub(pts).valid);
}

TEST_CASE("B0 stub fit rejects a rank-deficient quadratic", "[B0TrackerStubSeeder]") {
  std::vector<Point3> pts{{.x = 0.0, .y = 0.0, .z = 5902.0},
                          {.x = 1.0, .y = 0.0, .z = 5902.0},
                          {.x = 2.0, .y = 0.0, .z = 6172.0},
                          {.x = 3.0, .y = 0.0, .z = 6172.0}};
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
    CHECK(fit.covarianceX.at(3 * i + i) > 0.0);
  }
  for (int i = 0; i < 2; ++i) {
    CHECK(fit.covarianceY.at(2 * i + i) > 0.0);
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

TEST_CASE("B0 bend fit recovers signed q over p from the parabola", "[B0TrackerStubSeeder]") {
  constexpr double zEntrance  = 5800.0;
  constexpr double fieldY     = 1.184;
  constexpr double xEntrance  = -2.5;
  constexpr double txEntrance = -1.0e-3;
  constexpr double kBend      = 2.998e-4;
  for (const double qOverP : {1.0 / 41.0, -1.0 / 8.0}) {
    for (const auto& stations : kGeometries) {
      std::vector<Point3> points;
      for (double z : stations) {
        const double dz = z - zEntrance;
        points.push_back({xEntrance + txEntrance * dz - 0.5 * kBend * qOverP * fieldY * dz * dz,
                          -1.5 + 2.0e-3 * z, z, 4.0e-4, 9.0e-4});
      }
      const StubFit parabola = fitStub(points);
      REQUIRE(parabola.valid);
      const BendFit fit = bendFitFromParabola(parabola, fieldY, zEntrance);
      REQUIRE(fit.valid);
      REQUIRE(fit.covarianceValid);
      CHECK(fit.xReference == Approx(xEntrance).margin(1e-8));
      CHECK(fit.txReference == Approx(txEntrance).margin(1e-11));
      CHECK(fit.qOverP == Approx(qOverP).epsilon(1e-8));
      CHECK(fit.rmsX == Approx(0.0).margin(1e-8));

      auto scaledPoints = points;
      for (auto& point : scaledPoints) {
        point.varianceX *= 4.0;
      }
      const auto scaled = bendFitFromParabola(fitStub(scaledPoints), fieldY, zEntrance);
      REQUIRE(scaled.covarianceValid);
      CHECK(scaled.covariance[8] == Approx(4.0 * fit.covariance[8]).epsilon(1e-10));
      // Flipping the field flips the inferred charge.
      CHECK(bendFitFromParabola(parabola, -fieldY, zEntrance).qOverP ==
            Approx(-qOverP).epsilon(1e-8));
    }
  }
  CHECK_FALSE(
      bendFitFromParabola(fitStub(sampleTrack(0.0, 0.0, 0.0, 0.0, 0.0)), 0.0, zEntrance).valid);
}

TEST_CASE("B0 bend fit steps back through a lower-field gap exactly", "[B0TrackerStubSeeder]") {
  constexpr double zEntrance  = 5800.0;
  constexpr double fieldY     = 1.30;
  constexpr double fieldYGap  = 1.19;
  constexpr double qOverP     = 1.0 / 20.0;
  constexpr double xEntrance  = -2.5;
  constexpr double txEntrance = -1.0e-3;
  constexpr double kBend      = 2.998e-4;
  const double zFirst         = kOfficialStationZ.front();
  const double gap            = zFirst - zEntrance;
  // Piecewise-uniform truth: fieldYGap up to the first station, fieldY beyond.
  const double tx1 = txEntrance - kBend * qOverP * fieldYGap * gap;
  const double x1  = xEntrance + txEntrance * gap - 0.5 * kBend * qOverP * fieldYGap * gap * gap;
  std::vector<Point3> points;
  for (double z : kOfficialStationZ) {
    const double dz = z - zFirst;
    points.push_back(
        {x1 + tx1 * dz - 0.5 * kBend * qOverP * fieldY * dz * dz, 0.0, z, 4.0e-4, 4.0e-4});
  }
  const StubFit parabola = fitStub(points);
  REQUIRE(parabola.valid);
  const BendFit stepped = bendFitFromParabola(parabola, fieldY, zEntrance, zFirst, fieldYGap);
  REQUIRE(stepped.valid);
  REQUIRE(stepped.covarianceValid);
  CHECK(stepped.qOverP == Approx(qOverP).epsilon(1e-8));
  CHECK(stepped.xReference == Approx(xEntrance).margin(1e-8));
  CHECK(stepped.txReference == Approx(txEntrance).margin(1e-11));
  // Ignoring the gap field biases the entrance state.
  const BendFit uniform = bendFitFromParabola(parabola, fieldY, zEntrance);
  CHECK(std::abs(uniform.txReference - txEntrance) > 1e-7);
  // With an equal gap field the two forms agree, covariance included.
  const BendFit same = bendFitFromParabola(parabola, fieldY, zEntrance, zFirst, fieldY);
  CHECK(same.xReference == Approx(uniform.xReference).margin(1e-9));
  CHECK(same.txReference == Approx(uniform.txReference).margin(1e-12));
  for (int i = 0; i < 9; ++i) {
    CHECK(same.covariance.at(i) == Approx(uniform.covariance.at(i)).epsilon(1e-9).margin(1e-30));
  }
}

TEST_CASE("B0 endpoint pruning keeps prompt pairs and rejects cross-pairs",
          "[B0TrackerStubSeeder]") {
  const Point3 promptFirst{.x = 0.0, .y = 12.0, .z = 5902.0};
  const Point3 promptLast{.x = 0.0, .y = 12.0 * 6712.0 / 5902.0, .z = 6712.0};
  const auto prompt = endpointCompatibility(promptFirst, promptLast, 0.10, 5.0, true);
  REQUIRE(prompt.valid);
  CHECK(prompt.beamResidual == Approx(0.0).margin(1e-12));

  // This pair has a modest slope but extrapolates far from the beamline.
  const auto cross = endpointCompatibility({.x = 0.0, .y = 100.0, .z = 5902.0},
                                           {.x = 0.0, .y = 105.0, .z = 6712.0}, 0.10, 5.0, true);
  CHECK_FALSE(cross.valid);
  const auto unconstrained = endpointCompatibility(
      {.x = 0.0, .y = 100.0, .z = 5902.0}, {.x = 0.0, .y = 105.0, .z = 6712.0}, 0.10, 5.0, false);
  CHECK(unconstrained.valid);

  CHECK_FALSE(
      endpointCompatibility({0.0, 0.0, 5902.0}, {0.0, 100.0, 6000.0}, 0.10, 5.0, true).valid);
  // The reversed order is the point: a pair with the last hit first is not prompt.
  // NOLINTNEXTLINE(readability-suspicious-call-argument)
  CHECK_FALSE(endpointCompatibility(promptLast, promptFirst, 0.10, 5.0, true).valid);
}

TEST_CASE("B0 scattering covariance additions scale with q/p and project onto phi",
          "[B0TrackerStubSeeder]") {
  constexpr double qOverP = 0.04;
  constexpr double theta  = 0.03;
  const auto additions    = scatteringCovarianceAdditions(qOverP, theta, 1.0e-3, 0.02);
  const double angle2     = std::pow(1.0e-3 * qOverP, 2);
  CHECK(additions[1] == Approx(angle2).epsilon(1e-12));
  CHECK(additions[0] == Approx(angle2 / std::pow(std::sin(theta), 2)).epsilon(1e-12));
  CHECK(additions[2] == Approx(std::pow(0.02 * qOverP, 2)).epsilon(1e-12));
  CHECK(scatteringCovarianceAdditions(-qOverP, theta, 1.0e-3, 0.02) == additions);
}

namespace {

/// Analytic B0pf-like combined-function field in the ion frame: hard-edged in
/// z, By = B0 + G x, Bx = G y (the DD4hep MultipoleMagnet). Units T, mm.
struct BoxQuadrupoleField {
  double zEntrance;
  double zExit;
  double dipole;
  double gradient;
  std::array<double, 3> operator()(double x, double y, double z) const {
    if (z < zEntrance || z > zExit) {
      return {0.0, 0.0, 0.0};
    }
    return {gradient * y, dipole + gradient * x, 0.0};
  }
};

/// Lorentz-force RK4 propagation of a charged particle from the origin, in the
/// ion frame, returning the positions at the requested stations.
std::vector<Point3> propagateTruth(const BoxQuadrupoleField& field, double qOverP, double tx0,
                                   double ty0, const std::array<double, 4>& stations,
                                   double variance) {
  constexpr double kBend = 2.998e-4;
  // State: position (0-2) and unit direction (3-5) in the ion frame.
  using State = Eigen::Matrix<double, 6, 1>;
  State state;
  state << 0.0, 0.0, 0.0, tx0, ty0, 1.0;
  state.tail<3>().normalize();
  const auto derivative = [&](const State& s) {
    const auto b = field(s(0), s(1), s(2));
    // d(dir)/ds = kappa (q/p) dir x B
    State d;
    d << s(3), s(4), s(5), kBend * qOverP * (s(4) * b[2] - s(5) * b[1]),
        kBend * qOverP * (s(5) * b[0] - s(3) * b[2]), kBend * qOverP * (s(3) * b[1] - s(4) * b[0]);
    return d;
  };
  std::vector<Point3> out;
  const double h   = 1.0;
  std::size_t next = 0;
  while (next < stations.size()) {
    const State k1       = derivative(state);
    const State k2       = derivative(state + 0.5 * h * k1);
    const State k3       = derivative(state + 0.5 * h * k2);
    const State k4       = derivative(state + h * k3);
    const State previous = state;
    state += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    while (next < stations.size() && state(2) >= stations.at(next)) {
      const double f = (stations.at(next) - previous(2)) / (state(2) - previous(2));
      out.push_back({previous(0) + f * (state(0) - previous(0)),
                     previous(1) + f * (state(1) - previous(1)), stations.at(next), variance,
                     variance});
      ++next;
    }
  }
  return out;
}

} // namespace

TEST_CASE("B0 seed recovers charge, momentum and direction of an RK4 truth track",
          "[B0TrackerStubSeeder]") {
  constexpr double zEntrance     = 5800.0;
  constexpr double crossingAngle = -0.025;
  const BoxQuadrupoleField field{
      .zEntrance = zEntrance, .zExit = zEntrance + 1200.0, .dipole = 1.184, .gradient = -8.12e-3};
  constexpr double kBend = 2.998e-4;

  for (const double charge : {1.0, -1.0}) {
    for (const double momentum : {41.0, 15.0}) {
      const double qOverP = charge / momentum;
      const double tx0    = -2.0e-3;
      const double ty0    = 1.2e-2; // ~75 mm at B0: the Bx = G y term is well visible
      const auto points   = propagateTruth(field, qOverP, tx0, ty0, kOfficialStationZ, 4.0e-4);
      REQUIRE(points.size() == 4);

      // Sample the field once on the fitted trajectory at mid z, like the seeder.
      const StubFit parabola = fitStub(points);
      REQUIRE(parabola.valid);
      const double zMid  = 0.5 * (points.front().z + points.back().z);
      const auto bMid    = field(parabola.x(zMid), parabola.y(zMid), zMid);
      const auto bendFit = bendFitFromParabola(parabola, bMid[1], zEntrance);
      REQUIRE(bendFit.valid);
      // Signed q/p from the actual Lorentz-force trajectory: a wrong sign
      // convention or a wrong kappa would fail here. The tolerance is the
      // uniform-field approximation of the quadrupole-shaped field.
      CHECK(bendFit.qOverP == Approx(qOverP).epsilon(3e-3));
      // Parabola residual of the ~5 % field variation across the stations:
      // a few microns against the 20 um hit resolution.
      CHECK(bendFit.rmsX < 1.0e-2);

      // The quadrupole bends the non-bend plane: a straight line no longer
      // fits, until the kappa (q/p) Bx term is removed.
      const StubFit rawLine = fitStub(points);
      const StubFit corrected =
          fitStub(removeNonBendCurvature(points, bMid[0], bendFit.qOverP, zEntrance));
      REQUIRE(corrected.valid);
      CHECK(rawLine.rmsY > 0.05);
      // Residual from Bx varying along the track; well below the 0.5 mm cut.
      CHECK(corrected.rmsY < 2.0e-2);

      // Seed direction on the origin perigee equals the true initial direction
      // in the lab frame; the prompt track has zero impact parameters.
      const auto seed        = seedParametersFromFit(bendFit, corrected, crossingAngle, true);
      const double ca        = std::cos(crossingAngle);
      const double sa        = std::sin(crossingAngle);
      const double dx        = tx0 * ca + sa;
      const double dz        = -tx0 * sa + ca;
      const double truePhi   = std::atan2(ty0, dx);
      const double trueTheta = std::acos(dz / std::sqrt(dx * dx + ty0 * ty0 + dz * dz));
      CHECK(seed[0] == Approx(0.0).margin(1e-9));
      CHECK(seed[1] == Approx(0.0).margin(1e-9));
      CHECK(seed[2] == Approx(truePhi).margin(2e-3));
      CHECK(seed[3] == Approx(trueTheta).margin(5e-6));
      CHECK(seed[4] == Approx(qOverP).epsilon(3e-3));
      (void)kBend;
    }
  }
}

TEST_CASE("B0 bend-fit covariance propagates to correlated seed parameters",
          "[B0TrackerStubSeeder]") {
  constexpr double zEntrance = 5800.0;
  constexpr double fieldY    = 1.184;
  constexpr double qOverP    = 1.0 / 41.0;
  constexpr double kBend     = 2.998e-4;
  std::vector<Point3> points;
  for (double z : kOfficialStationZ) {
    const double dz = z - zEntrance;
    points.push_back({-2.5 - 1.0e-3 * dz - 0.5 * kBend * qOverP * fieldY * dz * dz,
                      -1.5 + 2.0e-3 * z, z, 4.0e-4, 4.0e-4});
  }
  const StubFit nonBendFit = fitStub(points);
  const BendFit bendFit    = bendFitFromParabola(nonBendFit, fieldY, zEntrance);
  REQUIRE(nonBendFit.covarianceValid);
  REQUIRE(bendFit.covarianceValid);
  const auto flat = seedCovarianceFromFit(bendFit, nonBendFit, -0.025, true);

  Eigen::Matrix<double, 5, 5> covariance;
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      covariance(row, col) = flat.at(5 * row + col);
      CHECK(std::isfinite(flat.at(5 * row + col)));
      CHECK(flat.at(5 * row + col) == Approx(flat.at(5 * col + row)).margin(1e-15));
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

TEST_CASE("B0 corrected fit covariance matches original hit perturbations",
          "[B0TrackerStubSeeder]") {
  constexpr double zEntrance = 5800.0;
  constexpr double fieldY    = 1.184;
  constexpr double kBend     = 2.998e-4;
  using Matrix5              = Eigen::Matrix<double, 5, 5>;
  using Vector5              = Eigen::Matrix<double, 5, 1>;
  for (const int nHits : {3, 4}) {
    for (const double fieldX : {-0.6, 0.0, 0.6}) {
      for (const double qOverP : {-1.0 / 41.0, 1.0 / 41.0}) {
        for (const bool constrained : {false, true}) {
          CAPTURE(nHits, fieldX, qOverP, constrained);
          std::vector<Point3> points;
          for (int i = 0; i < nHits; ++i) {
            const double z  = kOfficialStationZ.at(i);
            const double dz = z - zEntrance;
            points.push_back(
                {-2.5 - 2.5 / zEntrance * dz - 0.5 * kBend * qOverP * fieldY * dz * dz,
                 75.0 + 75.0 / zEntrance * dz + 0.5 * kBend * qOverP * fieldX * dz * dz, z,
                 4.0e-4 * (i + 1), 4.0e-4 / (i + 1)});
          }
          const auto bend      = bendFitFromParabola(fitStub(points), fieldY, zEntrance);
          const auto corrected = fitWithNonBendCorrection(points, fieldX, bend.qOverP, zEntrance);
          REQUIRE(bend.covarianceValid);
          REQUIRE(corrected.covarianceValid);
          const auto flat = seedCovarianceFromFit(bend, corrected, -0.025, constrained);
          const Matrix5 reported =
              Eigen::Map<const Eigen::Matrix<double, 5, 5, Eigen::RowMajor>>(flat.data());
          // Differentiate the complete estimator from original, independent
          // sensor coordinates: each x perturbation must also redo the y correction.
          const auto state = [&](const std::vector<Point3>& hits) {
            const auto b          = bendFitFromParabola(fitStub(hits), fieldY, zEntrance);
            const auto y          = fitWithNonBendCorrection(hits, fieldX, b.qOverP, zEntrance);
            const auto parameters = seedParametersFromFit(b, y, -0.025, constrained);
            return Vector5(Eigen::Map<const Vector5>(parameters.data()));
          };
          Matrix5 expected = Matrix5::Zero();
          for (int i = 0; i < nHits; ++i) {
            for (const bool xCoordinate : {false, true}) {
              auto plus             = points;
              auto minus            = points;
              constexpr double step = 1.0e-4; // mm
              (xCoordinate ? plus.at(i).x : plus.at(i).y) += step;
              (xCoordinate ? minus.at(i).x : minus.at(i).y) -= step;
              Vector5 difference       = state(plus) - state(minus);
              difference(2)            = std::remainder(difference(2), 2.0 * std::acos(-1.0));
              const Vector5 derivative = difference / (2.0 * step);
              const double variance = xCoordinate ? points.at(i).varianceX : points.at(i).varianceY;
              expected += variance * derivative * derivative.transpose();
            }
          }
          for (int row = 0; row < 5; ++row) {
            for (int col = 0; col < 5; ++col) {
              const double scale = std::sqrt(expected(row, row) * expected(col, col));
              CHECK(reported(row, col) ==
                    Approx(expected(row, col)).margin(2.0e-5 * scale + 1.0e-15));
            }
          }
        }
      }
    }
  }
}

TEST_CASE("B0 beam spot preserves missing measurement covariance fallbacks",
          "[B0TrackerStubSeeder]") {
  auto points              = sampleTrack(12.5, -3.0e-3, -4.0e-6, -1.5, 2.0e-3, 4.0e-4, 4.0e-4);
  bool haveVariances       = true;
  bool constrainToBeamline = true;
  double beamSize          = 1.0;
  SECTION("valid measurement covariance is unchanged") {}
  SECTION("missing hit variances retain angular fallbacks and beam correlations") {
    haveVariances = false;
  }
  SECTION("a fixed vertex retains the location fallbacks") { beamSize = 0.0; }
  SECTION("unconstrained seeds with missing variances retain all fallbacks") {
    haveVariances       = false;
    constrainToBeamline = false;
  }
  if (!haveVariances) {
    for (auto& point : points) {
      point.varianceX = 0.0;
      point.varianceY = 0.0;
    }
  }
  const auto fit  = fitStub(points);
  const auto bend = bendFitFromParabola(fit, 1.184, 5800.0);
  REQUIRE(fit.valid);
  REQUIRE(bend.valid);
  REQUIRE(fit.covarianceValid == haveVariances);
  REQUIRE(bend.covarianceValid == haveVariances);
  const auto measurement = seedCovarianceFromFit(bend, fit, -0.025, constrainToBeamline);
  const auto beam = beamSpotCovarianceAdditions(bend, fit, -0.025, constrainToBeamline,
                                                beamSize * 0.17, beamSize * 0.02, beamSize * 37.0);
  const std::array<double, 5> fallback{4.0, 1600.0, 0.01, 1.0e-5, 2.0e-4};
  const auto combined = seedCovarianceWithFallbacks(measurement, beam, fallback);
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      const auto index = 5 * row + col;
      double expected  = measurement.at(index) + beam.at(index);
      if (row == col) {
        if (row >= 2 && !haveVariances) {
          expected += fallback.at(row);
        } else if (expected == 0.0) {
          expected = fallback.at(row);
        }
      }
      CHECK(combined.at(index) == Approx(expected).margin(1e-15));
    }
  }
}

TEST_CASE("B0 perigee parameters of a ray through the origin", "[B0TrackerStubSeeder]") {
  // A ray leaving the origin at 25 mrad in the x-z plane: it passes through the
  // perigee line, so both local coordinates must vanish.
  const double slope = 0.025;
  const Point3 dir{.x = slope, .y = 0.0, .z = 1.0};
  const auto par =
      perigeeFromRay({.x = 0.0, .y = 0.0, .z = 0.0}, dir, {.x = 0.0, .y = 0.0, .z = 0.0});

  CHECK(par.loc0 == Approx(0.0).margin(1e-9));
  CHECK(par.loc1 == Approx(0.0).margin(1e-9));
  CHECK(par.phi == Approx(0.0).margin(1e-9));
  CHECK(par.theta == Approx(std::atan(slope)).epsilon(1e-9));
}

TEST_CASE("B0 perigee parameters of a displaced ray", "[B0TrackerStubSeeder]") {
  // Parallel to z, offset by +3 mm in y: the impact parameter is 3 mm and the
  // track is at phi = pi/2 relative to the point of closest approach.
  const auto par = perigeeFromRay({.x = 0.0, .y = 3.0, .z = 1000.0}, {.x = 0.0, .y = 0.0, .z = 1.0},
                                  {.x = 0.0, .y = 0.0, .z = 0.0});
  CHECK(std::abs(par.loc0) == Approx(3.0).epsilon(1e-9));
  CHECK(par.theta == Approx(0.0).margin(1e-9));
}

TEST_CASE("B0 perigee follows the surface it is expressed on", "[B0TrackerStubSeeder]") {
  // The same ray, expressed on a perigee surface at the field entrance rather
  // than at the origin: loc1 is measured from the surface centre.
  const Point3 dir{.x = 0.025, .y = 0.0, .z = 1.0};
  const auto atOrigin =
      perigeeFromRay({.x = 0.0, .y = 0.0, .z = 0.0}, dir, {.x = 0.0, .y = 0.0, .z = 0.0});
  const auto atEntrance =
      perigeeFromRay({.x = 0.0, .y = 0.0, .z = 0.0}, dir, {.x = 0.0, .y = 0.0, .z = 5800.0});

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
      CHECK(clustered[i].zMean == Approx(stations.at(i)).margin(1e-9));
      CHECK(assignStation(stations.at(i), clustered, gap) == static_cast<int>(i));
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
    CHECK(clustered[i].zMean == Approx(kRealisticStationZ.at(i)).margin(1e-9));
    CHECK(clustered[i].zMax - clustered[i].zMin == Approx(7.0).margin(1e-9));
  }
}

TEST_CASE("B0 cached stations drop mid-gap hits instead of inventing a station",
          "[B0TrackerStubSeeder]") {
  constexpr double gap = 50.0;
  const auto stations =
      clusterStations({kRealisticStationZ.begin(), kRealisticStationZ.end()}, gap);
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

TEST_CASE("B0 charge hypotheses cover the unresolved-curvature limit", "[B0TrackerStubSeeder]") {
  using eicrecon::b0stub::chargeHypotheses;
  const double minSignificance = 3.0;

  SECTION("a significant negative curvature resolves to one negative hypothesis") {
    const auto charges = chargeHypotheses(-0.1, 1.0e-4, minSignificance, -1, false, 0);
    REQUIRE(charges == std::vector<int>{-1});
  }
  SECTION("a significant positive curvature resolves to one positive hypothesis") {
    const auto charges = chargeHypotheses(0.1, 1.0e-4, minSignificance, +1, false, 0);
    REQUIRE(charges == std::vector<int>{1});
  }
  SECTION("a low-significance curvature emits both hypotheses") {
    const auto charges = chargeHypotheses(0.001, 1.0e-4, minSignificance, +1, false, 0);
    REQUIRE(charges == std::vector<int>{-1, 1});
  }
  SECTION("exactly zero curvature emits both hypotheses") {
    // Regression: the old `qOverP != 0` guard let this fall through to a
    // ternary on zero, silently emitting only the positive hypothesis at the
    // point of maximal charge ambiguity.
    const auto charges = chargeHypotheses(0.0, 1.0e-4, minSignificance, +1, false, 0);
    REQUIRE(charges == std::vector<int>{-1, 1});
  }
  SECTION("a zero variance is unresolved, not infinitely significant") {
    const auto charges = chargeHypotheses(0.1, 0.0, minSignificance, +1, false, 0);
    REQUIRE(charges == std::vector<int>{-1, 1});
  }
  SECTION("a negative variance is unresolved") {
    const auto charges = chargeHypotheses(0.1, -1.0, minSignificance, +1, false, 0);
    REQUIRE(charges == std::vector<int>{-1, 1});
  }
  SECTION("a non-finite curvature or variance is unresolved") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(chargeHypotheses(nan, 1.0e-4, minSignificance, +1, false, 0) ==
            std::vector<int>{-1, 1});
    REQUIRE(chargeHypotheses(0.1, nan, minSignificance, +1, false, 0) == std::vector<int>{-1, 1});
  }
  SECTION("an explicitly configured charge overrides the curvature") {
    REQUIRE(chargeHypotheses(0.0, 1.0e-4, minSignificance, +1, false, -1) == std::vector<int>{-1});
    REQUIRE(chargeHypotheses(-0.1, 1.0e-4, minSignificance, -1, false, +1) == std::vector<int>{1});
  }
  SECTION("testBothCharges wins over everything") {
    REQUIRE(chargeHypotheses(0.1, 1.0e-9, minSignificance, +1, true, -1) ==
            std::vector<int>{-1, 1});
  }
}

TEST_CASE("B0 beam-spot covariance is dominated by the longitudinal term",
          "[B0TrackerStubSeeder]") {
  constexpr double kBend = 2.998e-4;
  const double zEntrance = 5800.0;
  const double fieldY    = 1.184;
  const double qOverP    = 1.0 / 41.0;
  std::vector<Point3> points;
  for (double z : kOfficialStationZ) {
    const double dz = z - zEntrance;
    points.push_back({.x         = -2.5 - 1.0e-3 * dz - 0.5 * kBend * qOverP * fieldY * dz * dz,
                      .y         = -1.5 + 2.0e-3 * z,
                      .z         = z,
                      .varianceX = 4.0e-4,
                      .varianceY = 4.0e-4});
  }
  const StubFit nonBendFit = fitStub(points);
  const BendFit bendFit    = bendFitFromParabola(nonBendFit, fieldY, zEntrance);
  REQUIRE(bendFit.covarianceValid);
  const double theta = seedParametersFromFit(bendFit, nonBendFit, -0.025, true)[3];

  const auto zero = beamSpotCovarianceAdditions(bendFit, nonBendFit, -0.025, true, 0.0, 0.0, 0.0);
  CHECK(std::ranges::all_of(zero, [](double v) { return v == 0.0; }));

  // The unconstrained seed measures the direction; the vertex never enters.
  const auto unconstrained =
      beamSpotCovarianceAdditions(bendFit, nonBendFit, -0.025, false, 0.17, 0.02, 37.0);
  CHECK(std::ranges::all_of(unconstrained, [](double v) { return v == 0.0; }));

  // A purely longitudinal spread: theta scales with the lever arm, loc1 is the
  // displacement itself, and the two are perfectly correlated because one
  // variate drives both.
  const double sigmaZ = 37.0;
  const auto longitudinal =
      beamSpotCovarianceAdditions(bendFit, nonBendFit, -0.025, true, 0.0, 0.0, sigmaZ);
  const double sigmaTheta = std::sqrt(longitudinal.at(5 * 3 + 3));
  const double sigmaLoc1  = std::sqrt(longitudinal.at(5 * 1 + 1));
  CHECK(sigmaTheta == Approx(theta * sigmaZ / zEntrance).epsilon(0.05));
  CHECK(sigmaLoc1 == Approx(sigmaZ).epsilon(0.05));
  CHECK(std::abs(longitudinal.at(5 * 1 + 3)) == Approx(sigmaLoc1 * sigmaTheta).epsilon(1e-6));
  CHECK(longitudinal.at(5 * 0 + 0) < 1.0e-12);
  CHECK(longitudinal.at(5 * 2 + 2) < 1.0e-12);

  // A transverse spread moves the impact parameter and the azimuth instead.
  const auto transverse =
      beamSpotCovarianceAdditions(bendFit, nonBendFit, -0.025, true, 0.17, 0.02, 0.0);
  CHECK(std::sqrt(transverse.at(0)) > 0.5 * 0.02);
  CHECK(std::sqrt(transverse.at(0)) < 2.0 * 0.17);
  CHECK(transverse.at(5 * 2 + 2) > 0.0);
  CHECK(transverse.at(5 * 3 + 3) < longitudinal.at(5 * 3 + 3));

  // Quadratic in the sizes, symmetric, positive semi-definite, and q/p free.
  const auto doubled =
      beamSpotCovarianceAdditions(bendFit, nonBendFit, -0.025, true, 0.34, 0.04, 74.0);
  Eigen::Matrix<double, 5, 5> matrix;
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      matrix(row, col) = doubled.at(5 * row + col);
      CHECK(doubled.at(5 * row + col) == Approx(doubled.at(5 * col + row)).margin(1e-15));
      CHECK(doubled.at(5 * row + col) ==
            Approx(4.0 * (longitudinal.at(5 * row + col) + transverse.at(5 * row + col)))
                .epsilon(0.02)
                .margin(1e-15));
    }
    CHECK(doubled.at(5 * row + 4) == 0.0);
    CHECK(doubled.at(5 * 4 + row) == 0.0);
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 5, 5>> solver(matrix);
  REQUIRE(solver.info() == Eigen::Success);
  CHECK(solver.eigenvalues().minCoeff() >= Approx(0.0).margin(1e-12));
}
