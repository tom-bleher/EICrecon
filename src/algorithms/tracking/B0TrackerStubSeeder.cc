// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
//
// Stub seeder for the B0 tracker.
//
// Why this exists: the B0 sits inside the B0pf dipole (B along y). The
// orthogonal seeder's solenoid helix model returns invalid parameters there
// (charge sign is a coin flip, |p| spans 30 MeV - 77 GeV on 41 GeV protons),
// so the CKF can never attach hits. Here we instead exploit the spectrometer
// layout directly:
//
//   1. rotate hits into the ion frame (crossing angle about y), where the
//      dipole bend lies purely in the x-z plane;
//   2. fit y(z) with a line and x(z) with a parabola; with the field (Bx, By)
//      sampled on the fitted path at mid-track, the parabola curvature is the
//      signed q/p and the line, after removing the quadrupole bending
//      kappa (q/p) Bx, is the field-free upstream trajectory;
//   3. step the state back from the first station to the magnet face with a
//      second field sample in that gap, then express the result on the origin
//      perigee surface, which is what CKFTracking hard-codes as seed reference.
//
// The B0pf is modelled in DD4hep as a hard-edged, analytic dipole plus
// quadrupole; the field varies by ~15 % along a track through the four
// stations, yet the mid-track sample recovers q/p to ~0.1 % against RK4
// truth, far below the 2-4 % hit-resolution term, and on a proton-gun sample
// the two-sample model reproduces the seed resolution of a full sampled field
// integral (q/p 3.0 %, direction ~18 urad, pulls ~1).
//
// The seed is only a starting point: the CKF refits with the full field and
// material, so the stub q/p must not be read as a final measurement. Because
// the seed is expressed at the origin it is a prompt-track solution; see
// constrainToBeamline in the config for why displaced decays need more. The
// IP-to-B0pf path is assumed field-free, which holds for the ip6_extended
// configurations; with a solenoid present the constrained direction picks up
// a vertical bias of roughly 0.3 * integral(B_z) * sin(alpha) / p, and init()
// warns about it.

#include "B0TrackerStubSeeder.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DD4hep/VolumeManager.h>
#include <DDRec/CellIDPositionConverter.h>
#include <Eigen/Dense>
#include <algorithms/geo.h>
#include <edm4eic/Cov6f.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "algorithms/interfaces/ActsSvc.h"

namespace eicrecon {

namespace {
  /// Row-major views of the flat covariance arrays in StubFit and BendFit.
  using RowMajor2 = Eigen::Matrix<double, 2, 2, Eigen::RowMajor>;
  using RowMajor3 = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
  using RowMajor5 = Eigen::Matrix<double, 5, 5, Eigen::RowMajor>;
} // namespace

namespace b0stub {

  /// Bend constant: x'' [1/mm] = kBend * (q/p) [1/GeV] * B [T]
  constexpr double kBend = 2.998e-4;

  StubFit fitStub(const std::vector<Point3>& pts) {
    StubFit fit;
    const std::size_t n = pts.size();
    if (n < 3) {
      return fit;
    }

    // Centre and scale z before fitting. B0 measurements sit near z = 6 m;
    // fitting powers through z^4 in normal equations is needlessly ill
    // conditioned. QR on u = (z - zRef)/zScale gives the same polynomial with
    // stable coefficients.
    double zRef = 0.0;
    for (const auto& p : pts) {
      zRef += p.z;
    }
    zRef /= static_cast<double>(n);
    double zScale = 0.0;
    for (const auto& p : pts) {
      zScale = std::max(zScale, std::abs(p.z - zRef));
    }
    if (!(zScale > 0.0)) {
      return fit;
    }
    fit.zRef   = zRef;
    fit.zScale = zScale;

    Eigen::MatrixXd designX(n, 3);
    Eigen::MatrixXd designY(n, 2);
    Eigen::VectorXd valuesX(n);
    Eigen::VectorXd valuesY(n);
    bool haveVariances = true;
    for (std::size_t row = 0; row < n; ++row) {
      const double u  = (pts[row].z - zRef) / zScale;
      designX(row, 0) = 1.0;
      designX(row, 1) = u;
      designX(row, 2) = u * u;
      valuesX(row)    = pts[row].x;
      designY(row, 0) = 1.0;
      designY(row, 1) = u;
      valuesY(row)    = pts[row].y;
      haveVariances &= std::isfinite(pts[row].varianceX) && pts[row].varianceX > 0.0 &&
                       std::isfinite(pts[row].varianceY) && pts[row].varianceY > 0.0;
    }

    Eigen::MatrixXd weightedX       = designX;
    Eigen::MatrixXd weightedY       = designY;
    Eigen::VectorXd weightedValuesX = valuesX;
    Eigen::VectorXd weightedValuesY = valuesY;
    if (haveVariances) {
      for (std::size_t row = 0; row < n; ++row) {
        const double wx = 1.0 / std::sqrt(pts[row].varianceX);
        const double wy = 1.0 / std::sqrt(pts[row].varianceY);
        weightedX.row(row) *= wx;
        weightedY.row(row) *= wy;
        weightedValuesX(row) *= wx;
        weightedValuesY(row) *= wy;
      }
    }

    const auto decompositionX = weightedX.colPivHouseholderQr();
    const auto decompositionY = weightedY.colPivHouseholderQr();
    if (decompositionX.rank() < 3 || decompositionY.rank() < 2) {
      return fit;
    }
    const Eigen::Vector3d ax                       = decompositionX.solve(weightedValuesX);
    const Eigen::Vector2d ay                       = decompositionY.solve(weightedValuesY);
    Eigen::Map<Eigen::Vector3d>(fit.basisX.data()) = ax;
    Eigen::Map<Eigen::Vector2d>(fit.basisY.data()) = ay;

    // Convert back from the shifted/scaled variable to plain powers of z.
    fit.c2 = ax(2) / (zScale * zScale);
    fit.c1 = ax(1) / zScale - 2.0 * zRef * fit.c2;
    fit.c0 = ax(0) - ax(1) * zRef / zScale + ax(2) * zRef * zRef / (zScale * zScale);
    fit.b1 = ay(1) / zScale;
    fit.b0 = ay(0) - ay(1) * zRef / zScale;

    if (haveVariances) {
      const Eigen::Matrix3d normalX = weightedX.transpose() * weightedX;
      const Eigen::Matrix2d normalY = weightedY.transpose() * weightedY;
      const Eigen::Matrix3d covAx   = normalX.ldlt().solve(Eigen::Matrix3d::Identity());
      const Eigen::Matrix2d covAy   = normalY.ldlt().solve(Eigen::Matrix2d::Identity());

      if (covAx.allFinite() && covAy.allFinite() && (covAx.diagonal().array() > 0.0).all() &&
          (covAy.diagonal().array() > 0.0).all()) {
        Eigen::Map<RowMajor3>(fit.covarianceX.data()) = covAx;
        Eigen::Map<RowMajor2>(fit.covarianceY.data()) = covAy;
        fit.covarianceValid                           = true;
      }
    }

    double sumRx2 = 0.0;
    double sumRy2 = 0.0;
    for (const auto& p : pts) {
      const double rx = p.x - fit.x(p.z);
      const double ry = p.y - fit.y(p.z);
      sumRx2 += rx * rx;
      sumRy2 += ry * ry;
    }
    fit.rmsX = std::sqrt(sumRx2 / static_cast<double>(n));
    fit.rmsY = std::sqrt(sumRy2 / static_cast<double>(n));

    fit.valid = std::isfinite(fit.c0) && std::isfinite(fit.c1) && std::isfinite(fit.c2) &&
                std::isfinite(fit.b0) && std::isfinite(fit.b1) && std::isfinite(fit.rmsX) &&
                std::isfinite(fit.rmsY);
    return fit;
  }

  BendFit bendFitFromParabola(const StubFit& fit, double fieldY, double zReference, double zFirst,
                              double fieldYGap) {
    BendFit bend;
    bend.zReference = zReference;
    if (!fit.valid || !std::isfinite(fieldY) || fieldY == 0.0 || !std::isfinite(fieldYGap) ||
        !std::isfinite(zReference) || !std::isfinite(zFirst) || !(fit.zScale > 0.0)) {
      return bend;
    }
    // x'' = -kappa (q/p) B_y with x = c0 + c1 z + c2 z^2, so q/p = -2 c2 / (kappa B_y).
    const double curvatureToQOverP = -2.0 / (kBend * fieldY);
    const double s                 = fit.zScale;
    const double u1                = (zFirst - fit.zRef) / s;
    const double gap               = zFirst - zReference;
    // State at the first hit from the parabola, then a uniform-field step back
    // through the gap: tx_ref = tx1 + kappa (q/p) B_gap L,
    // x_ref = x1 - tx1 L - kappa (q/p) B_gap L^2 / 2.
    const double kGap = kBend * fieldYGap * gap;
    bend.qOverP       = curvatureToQOverP * fit.c2;
    const double x1   = fit.x(zFirst);
    const double tx1  = fit.tx(zFirst);
    bend.txReference  = tx1 + kGap * bend.qOverP;
    bend.xReference   = x1 - tx1 * gap - 0.5 * kGap * gap * bend.qOverP;
    bend.rmsX         = fit.rmsX;

    if (fit.covarianceValid) {
      // Linear map from the scaled basis coefficients (a0, a1, a2) with
      // u = (z - zRef)/zScale to (x_ref, tx_ref, q/p).
      const double dQdA2 = curvatureToQOverP / (s * s);
      Eigen::Matrix3d jacobian;
      jacobian << 1.0, u1 - gap / s, u1 * u1 - 2.0 * u1 * gap / s - 0.5 * kGap * gap * dQdA2, //
          0.0, 1.0 / s, 2.0 * u1 / s + kGap * dQdA2,                                          //
          0.0, 0.0, dQdA2;
      const Eigen::Matrix3d covAx      = Eigen::Map<const RowMajor3>(fit.covarianceX.data());
      const Eigen::Matrix3d covariance = jacobian * covAx * jacobian.transpose();
      if (covariance.allFinite() && (covariance.diagonal().array() > 0.0).all()) {
        Eigen::Map<RowMajor3>(bend.covariance.data()) = covariance;
        bend.covarianceValid                          = true;
      }
    }

    bend.valid = std::isfinite(bend.xReference) && std::isfinite(bend.txReference) &&
                 std::isfinite(bend.qOverP) && std::isfinite(bend.rmsX);
    return bend;
  }

  BendFit bendFitFromParabola(const StubFit& fit, double fieldY, double zReference) {
    return bendFitFromParabola(fit, fieldY, zReference, zReference, fieldY);
  }

  double nonBendCurvatureBound(double fieldX, double fieldY, double fieldZ, double maxAbsSlope,
                               double pMin) {
    if (!(pMin > 0.0) || !std::isfinite(pMin) || !std::isfinite(fieldX) || !std::isfinite(fieldY) ||
        !std::isfinite(fieldZ) || !(maxAbsSlope > 0.0) || !std::isfinite(maxAbsSlope)) {
      return std::numeric_limits<double>::infinity();
    }
    // y'' = kappa (q/p) sqrt(1+tx^2+ty^2)
    //       * [(1+ty^2) Bx - tx ty By - tx Bz].
    const double slope2 = maxAbsSlope * maxAbsSlope;
    return kBend / pMin * std::sqrt(1.0 + 2.0 * slope2) *
           ((1.0 + slope2) * std::abs(fieldX) + slope2 * std::abs(fieldY) +
            maxAbsSlope * std::abs(fieldZ));
  }

  EndpointCompatibility endpointCompatibility(const Point3& first, const Point3& last,
                                              double maxAbsSlope, double maxBeamResidual,
                                              bool constrainToBeamline, double maxAbsCurvatureY,
                                              double zFieldEntrance) {
    EndpointCompatibility result;
    const double dz = last.z - first.z;
    if (!(dz > 0.0) || !std::isfinite(first.y) || !std::isfinite(first.z) ||
        !std::isfinite(last.y) || !std::isfinite(last.z) || !(maxAbsSlope > 0.0) ||
        !(maxAbsCurvatureY >= 0.0) || !std::isfinite(zFieldEntrance) || !(first.z > 0.0)) {
      return result;
    }
    const double zEntrance   = std::clamp(zFieldEntrance, 0.0, first.z);
    const double slopeWindow = maxAbsSlope + maxAbsCurvatureY * (last.z - zEntrance);
    result.slopeY = (last.y - first.y) / dz;
    if (!std::isfinite(result.slopeY) || std::abs(result.slopeY) > slopeWindow) {
      return result;
    }

    double beamWindow = maxBeamResidual;
    if (constrainToBeamline && maxBeamResidual > 0.0) {
      // Integrate a bounded acceleration from the field entrance. The
      // resulting endpoint residual is bounded by
      // A/2 * [(zLast-zEntrance)^2 - zLast/zFirst*(zFirst-zEntrance)^2].
      beamWindow += 0.5 * maxAbsCurvatureY * dz * (last.z - zEntrance * zEntrance / first.z);
      result.beamResidual = last.y - first.y * last.z / first.z;
      if (!std::isfinite(result.beamResidual) || std::abs(result.beamResidual) > beamWindow) {
        return result;
      }
    }
    const double slopeScore = std::abs(result.slopeY) / slopeWindow;
    const double beamScore  = constrainToBeamline && maxBeamResidual > 0.0
                                  ? std::abs(result.beamResidual) / beamWindow
                                  : 0.0;
    result.score            = beamScore + slopeScore;
    result.valid            = true;
    return result;
  }

  std::vector<Point3> removeNonBendCurvature(const std::vector<Point3>& pts, double fieldX,
                                             double qOverP, double zReference) {
    std::vector<Point3> corrected = pts;
    if (!std::isfinite(fieldX) || !std::isfinite(qOverP) || !std::isfinite(zReference)) {
      return corrected;
    }
    // y'' = +kappa (q/p) Bx, the counterpart of x'' = -kappa (q/p) By.
    for (auto& point : corrected) {
      const double dz = point.z - zReference;
      point.y -= 0.5 * kBend * qOverP * fieldX * dz * dz;
    }
    return corrected;
  }

  StubFit fitWithNonBendCorrection(const std::vector<Point3>& pts, double fieldX, double qOverP,
                                   double zReference) {
    auto fit = fitStub(removeNonBendCurvature(pts, fieldX, qOverP, zReference));
    if (!fit.covarianceValid || !std::isfinite(fieldX) || !std::isfinite(qOverP) ||
        !std::isfinite(zReference)) {
      return fit;
    }
    // The weighted line estimator is (D^T W D)^-1 D^T W y. Apply the same
    // estimator to dy_i/d(q/p); the inverse normal matrix is covarianceY.
    Eigen::Vector2d response = Eigen::Vector2d::Zero();
    for (const auto& point : pts) {
      const double dz         = point.z - zReference;
      const double derivative = -0.5 * kBend * fieldX * dz * dz;
      const double u          = (point.z - fit.zRef) / fit.zScale;
      response += Eigen::Vector2d(1.0, u) * (derivative / point.varianceY);
    }
    Eigen::Map<Eigen::Vector2d>(fit.derivativeYQOverP.data()) =
        Eigen::Map<const RowMajor2>(fit.covarianceY.data()) * response;
    return fit;
  }

  std::array<double, 3> scatteringCovarianceAdditions(double qOverP, double theta,
                                                      double scatteringScale,
                                                      double qOverPRelativeUncertainty) {
    const double angle   = scatteringScale * std::abs(qOverP);
    const double sinTh   = std::max(std::abs(std::sin(theta)), 1.0e-6);
    const double relQopS = qOverPRelativeUncertainty * qOverP;
    return {angle * angle / (sinTh * sinTh), angle * angle, relQopS * relQopS};
  }

  PerigeeParams perigeeFromRay(const Point3& ref, const Point3& dir, const Point3& perigee) {
    PerigeeParams out;

    const double dnorm = std::sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
    if (!(dnorm > 0.0)) {
      return out;
    }
    const double dx = dir.x / dnorm;
    const double dy = dir.y / dnorm;
    const double dz = dir.z / dnorm;

    // Point of closest approach to the perigee line (parallel to z through
    // `perigee`), measured from `ref`.
    const double rx    = ref.x - perigee.x;
    const double ry    = ref.y - perigee.y;
    const double denom = dx * dx + dy * dy;
    const double t     = denom > 0.0 ? -(rx * dx + ry * dy) / denom : 0.0;

    const double xPca = rx + t * dx;
    const double yPca = ry + t * dy;
    const double zPca = ref.z + t * dz - perigee.z;

    out.phi   = std::atan2(dy, dx);
    out.theta = std::acos(std::clamp(dz, -1.0, 1.0));
    // ACTS perigee local frame: loc0 = signed transverse impact, loc1 = z
    out.loc0 = -xPca * std::sin(out.phi) + yPca * std::cos(out.phi);
    out.loc1 = zPca;
    return out;
  }

  std::vector<StationInterval> clusterStations(const std::vector<double>& zValues, double gap) {
    std::vector<double> ordered;
    ordered.reserve(zValues.size());
    for (const double z : zValues) {
      if (std::isfinite(z)) {
        ordered.push_back(z);
      }
    }
    std::vector<StationInterval> stations;
    if (ordered.empty() || !(gap > 0.0)) {
      return stations;
    }
    std::ranges::sort(ordered);

    StationInterval current{
        .zMin = ordered.front(), .zMax = ordered.front(), .zMean = ordered.front()};
    double sum    = ordered.front();
    std::size_t n = 1;
    for (std::size_t i = 1; i < ordered.size(); ++i) {
      if (ordered[i] - current.zMax > gap) {
        current.zMean = sum / static_cast<double>(n);
        stations.push_back(current);
        current = {.zMin = ordered[i], .zMax = ordered[i], .zMean = ordered[i]};
        sum     = ordered[i];
        n       = 1;
      } else {
        current.zMax = ordered[i];
        sum += ordered[i];
        ++n;
      }
    }
    current.zMean = sum / static_cast<double>(n);
    stations.push_back(current);
    return stations;
  }

  int assignStation(double z, const std::vector<StationInterval>& stations, double gap) {
    if (!std::isfinite(z) || !(gap > 0.0)) {
      return -1;
    }
    int best        = -1;
    double bestDist = 0.0;
    for (std::size_t i = 0; i < stations.size(); ++i) {
      if (z < stations[i].zMin - gap || z > stations[i].zMax + gap) {
        continue;
      }
      const double dist = std::abs(z - stations[i].zMean);
      if (best < 0 || dist < bestDist) {
        best     = static_cast<int>(i);
        bestDist = dist;
      }
    }
    return best;
  }

  std::map<unsigned int, std::vector<std::size_t>>
  groupHitsByStation(const std::vector<double>& hitZ, const std::vector<StationInterval>& stations,
                     double gap) {
    std::map<unsigned int, std::vector<std::size_t>> byStation;
    if (hitZ.empty()) {
      return byStation;
    }
    if (!stations.empty()) {
      for (std::size_t i = 0; i < hitZ.size(); ++i) {
        const int station = assignStation(hitZ[i], stations, gap);
        if (station >= 0) {
          byStation[static_cast<unsigned int>(station)].push_back(i);
        }
      }
      return byStation;
    }

    std::vector<std::size_t> order(hitZ.size());
    std::iota(order.begin(), order.end(), 0);
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return hitZ[a] < hitZ[b]; });
    unsigned int station = 0;
    double lastZ         = hitZ[order.front()];
    for (const std::size_t idx : order) {
      if (hitZ[idx] - lastZ > gap) {
        ++station;
      }
      lastZ = hitZ[idx];
      byStation[station].push_back(idx);
    }
    return byStation;
  }

} // namespace b0stub

namespace {

  struct StubCandidate {
    std::vector<std::size_t> hitIndices;
    b0stub::StubFit fit;
    b0stub::BendFit bendFit;
    double txFirst{0.0};
    double qOverP{0.0};
    double qOverPVariance{0.0};
    int inferredCharge{0};
  };

  using SeedVector = Eigen::Matrix<double, 5, 1>;
  using SeedMatrix = Eigen::Matrix<double, 5, 5>;

  struct EntranceFit {
    double x{};
    double tx{};
    double qOverP{};
    double y{};
    double ty{};
  };

  EntranceFit entranceFit(const b0stub::BendFit& bendFit, const b0stub::StubFit& nonBendFit) {
    const double u = (bendFit.zReference - nonBendFit.zRef) / nonBendFit.zScale;
    return {.x      = bendFit.xReference,
            .tx     = bendFit.txReference,
            .qOverP = bendFit.qOverP,
            .y      = nonBendFit.basisY[0] + nonBendFit.basisY[1] * u,
            .ty     = nonBendFit.basisY[1] / nonBendFit.zScale};
  }

  /// Beamline-constrained seed state for a track that starts at `vertex`
  /// (lab frame, mm) and passes through the field-entrance point. A zero
  /// vertex is the shipped constrained mode; the vertex is varied only to
  /// derive the beam-spot covariance.
  SeedVector seedStateFromVertex(const EntranceFit& entrance, double zFieldEntrance,
                                 double crossingAngle, const b0stub::Point3& vertex) {
    const double ca = std::cos(crossingAngle);
    const double sa = std::sin(crossingAngle);

    // Field-entrance point in the lab frame: the ion-frame point
    // (entrance.x, entrance.y, zFieldEntrance) rotated by the crossing angle.
    const b0stub::Point3 point{.x = entrance.x * ca + zFieldEntrance * sa,
                               .y = entrance.y,
                               .z = -entrance.x * sa + zFieldEntrance * ca};
    const b0stub::Point3 dir{
        .x = point.x - vertex.x, .y = point.y - vertex.y, .z = point.z - vertex.z};
    const auto perigee = b0stub::perigeeFromRay(vertex, dir, {.x = 0.0, .y = 0.0, .z = 0.0});

    SeedVector state;
    state << perigee.loc0, perigee.loc1, perigee.phi, perigee.theta, entrance.qOverP;
    return state;
  }

  SeedVector seedStateFromEntrance(const EntranceFit& entrance, double zFieldEntrance,
                                   double crossingAngle, bool constrainToBeamline) {
    if (constrainToBeamline) {
      return seedStateFromVertex(entrance, zFieldEntrance, crossingAngle,
                                 {.x = 0.0, .y = 0.0, .z = 0.0});
    }
    const double ca = std::cos(crossingAngle);
    const double sa = std::sin(crossingAngle);

    const double x0 = entrance.x - entrance.tx * zFieldEntrance;
    const double y0 = entrance.y - entrance.ty * zFieldEntrance;

    const b0stub::Point3 dir{
        .x = entrance.tx * ca + sa, .y = entrance.ty, .z = -entrance.tx * sa + ca};
    const b0stub::Point3 ref{.x = x0 * ca, .y = y0, .z = -x0 * sa};
    const auto perigee = b0stub::perigeeFromRay(ref, dir, {.x = 0.0, .y = 0.0, .z = 0.0});

    SeedVector state;
    state << perigee.loc0, perigee.loc1, perigee.phi, perigee.theta, entrance.qOverP;
    return state;
  }

  SeedVector seedStateFromFit(const b0stub::BendFit& bendFit, const b0stub::StubFit& nonBendFit,
                              double crossingAngle, bool constrainToBeamline) {
    return seedStateFromEntrance(entranceFit(bendFit, nonBendFit), bendFit.zReference,
                                 crossingAngle, constrainToBeamline);
  }

  SeedMatrix propagateFitCovariance(const b0stub::BendFit& bendFit,
                                    const b0stub::StubFit& nonBendFit, double crossingAngle,
                                    bool constrainToBeamline) {
    SeedMatrix result = SeedMatrix::Zero();
    if (!bendFit.covarianceValid || !nonBendFit.covarianceValid) {
      return result;
    }

    SeedMatrix fitCovariance          = SeedMatrix::Zero();
    fitCovariance.block<3, 3>(0, 0)   = Eigen::Map<const RowMajor3>(bendFit.covariance.data());
    const Eigen::Matrix2d yCovariance = Eigen::Map<const RowMajor2>(nonBendFit.covarianceY.data());
    const double uEntrance            = (bendFit.zReference - nonBendFit.zRef) / nonBendFit.zScale;
    Eigen::Matrix2d yTransform;
    yTransform << 1.0, uEntrance, 0.0, 1.0 / nonBendFit.zScale;
    fitCovariance.block<2, 2>(3, 3) = yTransform * yCovariance * yTransform.transpose();

    // Corrected y coordinates share the momentum inferred from x. Retain
    // that dependence, including the induced x/tx-to-y/ty cross terms.
    SeedMatrix correction = SeedMatrix::Identity();
    correction.block<2, 1>(3, 2) =
        yTransform * Eigen::Map<const Eigen::Vector2d>(nonBendFit.derivativeYQOverP.data());
    fitCovariance = (correction * fitCovariance * correction.transpose()).eval();

    const EntranceFit nominal = entranceFit(bendFit, nonBendFit);
    const std::array<double, 5> values{nominal.x, nominal.tx, nominal.qOverP, nominal.y,
                                       nominal.ty};
    SeedMatrix jacobian = SeedMatrix::Zero();
    for (int column = 0; column < 5; ++column) {
      const double sigma = std::sqrt(std::max(0.0, fitCovariance(column, column)));
      const double step = std::max({1.0e-10, sigma * 1.0e-3, std::abs(values.at(column)) * 1.0e-7});
      auto plus         = nominal;
      auto minus        = nominal;
      double* plusValue = nullptr;
      double* minusValue = nullptr;
      switch (column) {
      case 0:
        plusValue  = &plus.x;
        minusValue = &minus.x;
        break;
      case 1:
        plusValue  = &plus.tx;
        minusValue = &minus.tx;
        break;
      case 2:
        plusValue  = &plus.qOverP;
        minusValue = &minus.qOverP;
        break;
      case 3:
        plusValue  = &plus.y;
        minusValue = &minus.y;
        break;
      default:
        plusValue  = &plus.ty;
        minusValue = &minus.ty;
        break;
      }
      *plusValue += step;
      *minusValue -= step;
      const SeedVector upper =
          seedStateFromEntrance(plus, bendFit.zReference, crossingAngle, constrainToBeamline);
      const SeedVector lower =
          seedStateFromEntrance(minus, bendFit.zReference, crossingAngle, constrainToBeamline);
      SeedVector difference = upper - lower;
      difference(2)         = std::remainder(upper(2) - lower(2), 2.0 * std::acos(-1.0));
      jacobian.col(column)  = difference / (2.0 * step);
    }

    result = jacobian * fitCovariance * jacobian.transpose();
    if (!result.allFinite()) {
      return SeedMatrix::Zero();
    }
    return 0.5 * (result + result.transpose());
  }

} // namespace

namespace b0stub {

  std::array<double, 5> seedParametersFromFit(const BendFit& bendFit, const StubFit& nonBendFit,
                                              double crossingAngle, bool constrainToBeamline) {
    const SeedVector state =
        seedStateFromFit(bendFit, nonBendFit, crossingAngle, constrainToBeamline);
    return {state(0), state(1), state(2), state(3), state(4)};
  }

  std::array<double, 25> seedCovarianceFromFit(const BendFit& bendFit, const StubFit& nonBendFit,
                                               double crossingAngle, bool constrainToBeamline) {
    const SeedMatrix covariance =
        propagateFitCovariance(bendFit, nonBendFit, crossingAngle, constrainToBeamline);
    std::array<double, 25> result{};
    Eigen::Map<RowMajor5>(result.data()) = covariance;
    return result;
  }

  std::array<double, 25>
  seedCovarianceWithFallbacks(const std::array<double, 25>& fitCovariance,
                              const std::array<double, 25>& beamSpotCovariance,
                              const std::array<double, 5>& fallbackVariances) {
    SeedMatrix covariance    = Eigen::Map<const RowMajor5>(fitCovariance.data());
    const auto applyFallback = [&](int i) {
      if (!std::isfinite(covariance(i, i)) || covariance(i, i) <= 0.0) {
        covariance.row(i).setZero();
        covariance.col(i).setZero();
        covariance(i, i) = fallbackVariances.at(i);
      }
    };
    // Vertex uncertainty cannot replace missing measurement uncertainty on
    // the direction or momentum. Preserve it in addition to these fallbacks.
    for (int i = 2; i < 5; ++i) {
      applyFallback(i);
    }
    covariance += Eigen::Map<const RowMajor5>(beamSpotCovariance.data());
    for (int i = 0; i < 5; ++i) {
      applyFallback(i);
    }
    std::array<double, 25> result{};
    Eigen::Map<RowMajor5>(result.data()) = covariance;
    return result;
  }

  std::array<double, 25> beamSpotCovarianceAdditions(const BendFit& bendFit,
                                                     const StubFit& nonBendFit,
                                                     double crossingAngle, bool constrainToBeamline,
                                                     double sigmaX, double sigmaY, double sigmaZ) {
    std::array<double, 25> result{};
    if (!constrainToBeamline) {
      return result;
    }
    const EntranceFit entrance = entranceFit(bendFit, nonBendFit);
    const std::array<double, 3> sigmas{sigmaX, sigmaY, sigmaZ};
    SeedMatrix covariance = SeedMatrix::Zero();
    for (std::size_t axis = 0; axis < sigmas.size(); ++axis) {
      const double sigma = sigmas.at(axis);
      if (!(sigma > 0.0) || !std::isfinite(sigma)) {
        continue;
      }
      Point3 plus{};
      Point3 minus{};
      // One displacement of +-1 sigma along each lab axis in turn. The map is
      // linear to well within the beam-spot size, so the half-difference is
      // the column of the Jacobian already scaled by sigma.
      std::array<double*, 3> plusAxes{&plus.x, &plus.y, &plus.z};
      std::array<double*, 3> minusAxes{&minus.x, &minus.y, &minus.z};
      *plusAxes.at(axis)  = sigma;
      *minusAxes.at(axis) = -sigma;
      const SeedVector upper =
          seedStateFromVertex(entrance, bendFit.zReference, crossingAngle, plus);
      const SeedVector lower =
          seedStateFromVertex(entrance, bendFit.zReference, crossingAngle, minus);
      SeedVector shift = 0.5 * (upper - lower);
      shift(2)         = 0.5 * std::remainder(upper(2) - lower(2), 2.0 * std::acos(-1.0));
      shift(4)         = 0.0; // q/p comes from the curvature, not from the vertex
      covariance += shift * shift.transpose();
    }
    if (!covariance.allFinite()) {
      return {};
    }
    Eigen::Map<RowMajor5>(result.data()) = covariance;
    return result;
  }

  std::vector<int> chargeHypotheses(double qOverP, double qOverPVariance,
                                    double minCurvatureSignificance, int inferredCharge,
                                    bool testBothCharges, int configuredCharge) {
    if (testBothCharges) {
      return {-1, 1};
    }
    if (configuredCharge == -1 || configuredCharge == 1) {
      return {configuredCharge};
    }
    // A charge sign is only resolved by a curvature that is significantly
    // non-zero. Exactly zero curvature, and a variance that is non-positive or
    // non-finite, are both maximally ambiguous.
    const bool measurable     = std::isfinite(qOverP) && qOverP != 0.0 &&
                                std::isfinite(qOverPVariance) && qOverPVariance > 0.0;
    const double significance = measurable ? std::abs(qOverP) / std::sqrt(qOverPVariance) : 0.0;
    if (significance < minCurvatureSignificance) {
      return {-1, 1};
    }
    return {inferredCharge};
  }

} // namespace b0stub

void B0TrackerStubSeeder::init() {
  m_acts_context = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (m_acts_context == nullptr) {
    throw std::runtime_error("B0TrackerStubSeeder: ACTS geometry service is unavailable");
  }
  const auto* detector = m_acts_context->dd4hepDetector();
  if (detector == nullptr) {
    throw std::runtime_error("B0TrackerStubSeeder: DD4hep detector is unavailable");
  }
  // A geometry without the far-forward definitions has no B0 to seed. Disable
  // the seeder and emit empty collections rather than failing every event.
  m_enabled = true;
  try {
    m_crossing_angle = detector->constant<double>("CrossingAngle") / dd4hep::rad;
  } catch (const std::exception& error) {
    warning("B0TrackerStubSeeder: DD4hep constant CrossingAngle is unavailable ({}); "
            "the seeder is disabled and will produce no seeds",
            error.what());
    m_enabled = false;
    return;
  }
  if (!std::isfinite(m_crossing_angle)) {
    throw std::runtime_error("B0TrackerStubSeeder: DD4hep CrossingAngle is not finite");
  }

  const double ca = std::cos(m_crossing_angle);
  const double sa = std::sin(m_crossing_angle);
  if (m_cfg.zFieldEntrance > 0.0) {
    m_z_field_entrance = m_cfg.zFieldEntrance;
  } else {
    try {
      const double zCenter = detector->constant<double>("B0PF_CenterPosition") / dd4hep::mm;
      const double length  = detector->constant<double>("B0PF_Length") / dd4hep::mm;
      double xLab          = 0.0;
      try {
        xLab = detector->constant<double>("B0PF_XPosition") / dd4hep::mm;
      } catch (const std::exception&) {
        xLab = 0.0;
      }
      const double zLab  = zCenter - 0.5 * length;
      m_z_face_lab       = zLab;
      m_z_field_entrance = xLab * sa + zLab * ca;
    } catch (const std::exception& error) {
      warning("B0TrackerStubSeeder: cannot derive the B0pf entrance from the geometry "
              "({}); set zFieldEntrance explicitly or provide B0PF_CenterPosition and "
              "B0PF_Length. The seeder is disabled and will produce no seeds",
              error.what());
      m_enabled = false;
      return;
    }
  }
  if (!std::isfinite(m_z_field_entrance)) {
    throw std::runtime_error("B0TrackerStubSeeder: B0pf field entrance z is not finite");
  }

  // The seed is back-extrapolated along a straight line from the B0pf
  // entrance to the origin. Warn if the geometry puts a field on that path.
  {
    const auto fieldProvider = m_acts_context->getFieldProvider();
    auto fieldCache  = fieldProvider->makeCache(m_acts_context->getActsMagneticFieldContext());
    const auto field = fieldProvider->getField(Acts::Vector3::Zero(), fieldCache);
    if (field.ok() && field.value().allFinite()) {
      const double magnitude = field.value().norm() / Acts::UnitConstants::T;
      if (magnitude > 0.01) {
        warning("B0TrackerStubSeeder: |B| = {:.2f} T at the origin. The seeder assumes a "
                "field-free path from the IP to the B0pf, so a solenoid biases the "
                "beamline-constrained seed direction by ~0.3*int(B_z)*sin(alpha)/p; "
                "seed pulls will be wrong for low-momentum tracks.",
                magnitude);
      }
    }
  }

  m_stations.clear();
  m_volume_to_station.clear();
  std::vector<double> surfaceZ;
  std::vector<std::uint64_t> surfaceVolumes;
  auto volman = detector->volumeManager();
  for (const auto& [volId, surface] : m_acts_context->surfaceMap()) {
    if (surface == nullptr) {
      continue;
    }
    std::string path;
    try {
      const auto* volumeContext = volman.lookupContext(volId);
      if (volumeContext == nullptr) {
        continue;
      }
      path = volumeContext->element.path();
    } catch (const std::exception&) {
      continue;
    }
    if (path.find("B0Tracker") == std::string::npos ||
        path.find("Companion") != std::string::npos) {
      continue;
    }
    const auto center = surface->center(m_acts_context->getActsGeometryContext());
    const double zIon =
        center.x() / Acts::UnitConstants::mm * sa + center.z() / Acts::UnitConstants::mm * ca;
    if (!std::isfinite(zIon)) {
      continue;
    }
    surfaceZ.push_back(zIon);
    surfaceVolumes.push_back(volId);
  }
  m_stations = b0stub::clusterStations(surfaceZ, m_cfg.stationZGap);
  if (m_stations.empty()) {
    warning("B0TrackerStubSeeder: no B0 ACTS surfaces found; station grouping "
            "will fall back to per-event hit clustering");
  } else {
    for (std::size_t i = 0; i < surfaceZ.size(); ++i) {
      const int station = b0stub::assignStation(surfaceZ[i], m_stations, m_cfg.stationZGap);
      if (station >= 0) {
        m_volume_to_station[surfaceVolumes[i]] = static_cast<unsigned int>(station);
      }
    }
    debug("B0TrackerStubSeeder: cached {} B0 stations from {} ACTS surfaces; "
          "field entrance z = {:.3f} mm",
          m_stations.size(), surfaceZ.size(), m_z_field_entrance);
  }
}

void B0TrackerStubSeeder::process(const Input& input, const Output& output) const {
  const auto [hits]                 = input;
  auto [seeds, track_params_output] = output;

  if (!m_enabled || hits->empty()) {
    return;
  }

  const double ca = std::cos(m_crossing_angle);
  const double sa = std::sin(m_crossing_angle);

  // ------------------------------------------------------------------
  // Rotate hit positions into the ion frame and group them by station.
  // ------------------------------------------------------------------
  const auto* converter   = m_volume_to_station.empty()
                                ? nullptr
                                : algorithms::GeoSvc::instance().cellIDPositionConverter();
  const auto lookupVolume = [&](const edm4eic::TrackerHit& hit) -> std::uint64_t {
    if (converter == nullptr) {
      return 0;
    }
    try {
      const auto* volumeContext = converter->findContext(hit.getCellID());
      return volumeContext == nullptr ? 0 : volumeContext->identifier;
    } catch (const std::exception&) {
      return 0;
    }
  };
  std::vector<b0stub::Point3> ion;
  ion.reserve(hits->size());
  std::vector<edm4eic::TrackerHit> ionHits;
  ionHits.reserve(hits->size());
  std::vector<std::uint64_t> ionVolumes;
  ionVolumes.reserve(hits->size());

  for (const auto& hit : *hits) {
    const auto& p          = hit.getPosition();
    const auto& covariance = hit.getPositionError();
    const auto volumeId    = lookupVolume(hit);
    // The hit position covariance is diagonal in the sensor's local frame.
    // Every B0 sensor uses a square pitch -- 70 um in the realistic geometry,
    // and the same in x and z upstream -- so the two in-plane variances are
    // equal and carry over to the ion-frame bend and non-bend axes whatever
    // the module's azimuth. A rectangular pitch would have to be rotated.
    const b0stub::Point3 q{.x         = p.x * ca - p.z * sa,
                           .y         = p.y,
                           .z         = p.x * sa + p.z * ca,
                           .varianceX = covariance.xx,
                           .varianceY = covariance.yy};
    if (!(std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z))) {
      debug("Skipping B0 hit with non-finite position");
      continue;
    }
    ion.push_back(q);
    ionHits.push_back(hit);
    ionVolumes.push_back(volumeId);
  }
  if (ion.empty()) {
    return;
  }

  // Prefer the stations cached from B0 ACTS surfaces at init. Official B0
  // uses layer 1-4 for four disks; the realistic geometry uses layer 1-8
  // (front/back per disk), so cellID layer/2 is not portable. The cached
  // z-gap clustering is. Hits that do not match a surface station (or, if
  // no surfaces were found, the per-event fallback) are dropped.
  std::map<unsigned int, std::vector<std::size_t>> byStation;
  if (!m_stations.empty()) {
    for (std::size_t i = 0; i < ion.size(); ++i) {
      int station = -1;
      if (const auto found = m_volume_to_station.find(ionVolumes[i]);
          found != m_volume_to_station.end()) {
        station = static_cast<int>(found->second);
      }
      if (station < 0) {
        station = b0stub::assignStation(ion[i].z, m_stations, m_cfg.stationZGap);
      }
      if (station >= 0) {
        byStation[static_cast<unsigned int>(station)].push_back(i);
      }
    }
  } else {
    std::vector<double> hitZ(ion.size());
    for (std::size_t i = 0; i < ion.size(); ++i) {
      hitZ[i] = ion[i].z;
    }
    byStation = b0stub::groupHitsByStation(hitZ, {}, m_cfg.stationZGap);
  }
  for (auto& [stationIndex, indices] : byStation) {
    std::ranges::sort(indices, [&](std::size_t a, std::size_t b) {
      return std::tie(ion[a].x, ion[a].y, ion[a].z) < std::tie(ion[b].x, ion[b].y, ion[b].z);
    });
  }

  // The auxiliary and field-integral fits need >= 3 points, so 3 stations is the hard floor
  // regardless of configuration (a 2-station mode would need a line fit with
  // a momentum prior instead).
  const unsigned int minStations = std::max(m_cfg.minStations, 3U);
  if (byStation.size() < minStations) {
    return;
  }

  // Sample the same ACTS/DD4hep field provider used by CKFTracking. The B0
  // field has a dipole component and a gradient, so an analytic nominal B is
  // not portable across beam-energy configurations or hit positions.
  const auto fieldProvider = m_acts_context->getFieldProvider();
  auto fieldCache = fieldProvider->makeCache(m_acts_context->getActsMagneticFieldContext());

  // The quadrupole bends y before a candidate's q/p is known. Use the
  // minimum allowed momentum to make conservative preliminary roads, rather
  // than applying the prompt straight-line cut to uncorrected hits. Endpoint
  // field samples estimate the field bound for the smooth B0 multipole; they
  // are not a guaranteed bound for an arbitrary field map. The fitted,
  // signed-curvature correction below still has to pass the original cuts.
  std::vector<double> curvatureBounds(ion.size(), std::numeric_limits<double>::infinity());
  if (m_cfg.pMin > 0.0) {
    for (std::size_t i = 0; i < ion.size(); ++i) {
      const auto& point = ion[i];
      const auto field  = fieldProvider->getField(
          {point.x * ca + point.z * sa, point.y, -point.x * sa + point.z * ca}, fieldCache);
      if (field.ok()) {
        const Acts::Vector3 lab = field.value() / Acts::UnitConstants::T;
        curvatureBounds[i] = b0stub::nonBendCurvatureBound(lab.x() * ca - lab.z() * sa, lab.y(),
                                                           lab.x() * sa + lab.z() * ca,
                                                           m_cfg.maxAbsTransverseSlope, m_cfg.pMin);
      }
    }
  }

  // ------------------------------------------------------------------
  // Enumerate candidates: one hit per station, capped combinatorics.
  // In addition to the full station set, try leave-one-station-out
  // subsets (>= minStations): a proton that scattered hard at station k
  // has a clean subset excluding the post-kink stations, whose fit is
  // unbiased, while the all-station fit averages over the kink.
  // ------------------------------------------------------------------
  std::vector<std::vector<std::size_t>> allStationHits;
  allStationHits.reserve(byStation.size());
  for (const auto& [station, idxs] : byStation) {
    allStationHits.push_back(idxs);
  }

  std::vector<std::vector<std::vector<std::size_t>>> subsets;
  subsets.push_back(allStationHits);
  if (allStationHits.size() > minStations) {
    for (std::size_t drop = 0; drop < allStationHits.size(); ++drop) {
      std::vector<std::vector<std::size_t>> sub;
      for (std::size_t s = 0; s < allStationHits.size(); ++s) {
        if (s != drop) {
          sub.push_back(allStationHits[s]);
        }
      }
      subsets.push_back(sub);
    }
  }

  std::vector<StubCandidate> candidates;
  // Split the combination budget across subsets so a busy event cannot
  // exhaust it inside the full-station set and starve the leave-one-out
  // subsets that exist to recover kinked tracks.
  const unsigned int perSubsetBudget =
      std::max(1U, m_cfg.maxCombinations / static_cast<unsigned int>(subsets.size()));

  auto makeCompatible = [&](StubCandidate& cand) -> bool {
    // The y residual is only tested after the quadrupole bending has been
    // removed below.
    if (!cand.fit.valid || cand.fit.rmsX > m_cfg.maxXResidual) {
      return false;
    }

    const auto zMinMax = std::ranges::minmax_element(
        cand.hitIndices, [&](std::size_t a, std::size_t b) { return ion[a].z < ion[b].z; });
    const double zFirst = ion[*zMinMax.min].z;
    const double zLast  = ion[*zMinMax.max].z;
    // Ion-frame z at which this track crosses the magnet entrance face. The
    // face is a plane of constant lab z, so its ion z shifts by x_lab * sin(a)
    // from track to track; ignoring that mis-integrates the hard field edge
    // and biases the seed direction by tens of microradians.
    double zEntrance = m_z_field_entrance;
    if (m_z_face_lab > 0.0) {
      for (int pass = 0; pass < 2; ++pass) {
        const double xLab = cand.fit.x(zEntrance) * ca + zEntrance * sa;
        zEntrance         = xLab * sa + m_z_face_lab * ca;
      }
    }
    if (!std::isfinite(zEntrance) || !(zEntrance < zFirst)) {
      return false;
    }

    std::vector<b0stub::Point3> points;
    points.reserve(cand.hitIndices.size());
    for (const std::size_t idx : cand.hitIndices) {
      points.push_back(ion[idx]);
    }

    // Sample the same ACTS/DD4hep field provider used by CKFTracking once, on
    // the fitted trajectory at the candidate's mid z. The B0pf is a hard-edged
    // analytic dipole+quadrupole whose field is close to linear along a track,
    // so the mid-point sample is the mean field to second order: it recovers
    // q/p to ~0.1 % against RK4 truth, well below the hit-resolution term.
    const double zMid = 0.5 * (zFirst + zLast);
    const double xMid = cand.fit.x(zMid);
    const double yMid = cand.fit.y(zMid);
    const Acts::Vector3 global{xMid * ca + zMid * sa, yMid, -xMid * sa + zMid * ca};
    const auto field = fieldProvider->getField(global, fieldCache);
    if (!field.ok() || !field.value().allFinite()) {
      return false;
    }
    const Acts::Vector3 lab = field.value() / Acts::UnitConstants::T;
    // Rotate the field into the ion frame, like the positions.
    const double fieldX = lab.x() * ca - lab.z() * sa;
    const double fieldY = lab.y();
    if (std::abs(fieldY) < m_cfg.minAbsFieldY) {
      return false;
    }
    // The magnet is not rotated with the beam, so the dipole component along
    // the ion axis rises through the magnet: the gap between the entrance face
    // and the first station sits ~8 % below the mid-track field. One sample in
    // the middle of the gap keeps the back-extrapolated direction unbiased.
    const double zGap = 0.5 * (zEntrance + zFirst);
    const double xGap = cand.fit.x(zGap);
    const double yGap = cand.fit.y(zGap);
    const auto gapField =
        fieldProvider->getField({xGap * ca + zGap * sa, yGap, -xGap * sa + zGap * ca}, fieldCache);
    if (!gapField.ok() || !gapField.value().allFinite()) {
      return false;
    }
    const double fieldYGap = gapField.value().y() / Acts::UnitConstants::T;

    cand.bendFit = b0stub::bendFitFromParabola(cand.fit, fieldY, zEntrance, zFirst, fieldYGap);
    if (!cand.bendFit.valid) {
      return false;
    }

    // Remove the quadrupole bending from the non-bend coordinates and refit
    // the line, so that it describes the field-free upstream trajectory.
    cand.fit = b0stub::fitWithNonBendCorrection(points, fieldX, cand.bendFit.qOverP, zEntrance);
    if (!cand.fit.valid || cand.fit.rmsY > m_cfg.maxYResidual ||
        std::abs(cand.fit.b1) > m_cfg.maxAbsTransverseSlope) {
      return false;
    }
    // Now q/p is measured and the y fit describes the upstream ray. Enforce
    // the configured beamline tolerance without the preliminary envelope.
    if (!b0stub::endpointCompatibility(
             {.y = cand.fit.y(zFirst), .z = zFirst}, {.y = cand.fit.y(zLast), .z = zLast},
             m_cfg.maxAbsTransverseSlope, m_cfg.maxYBeamlineResidual, m_cfg.constrainToBeamline)
             .valid) {
      return false;
    }

    cand.txFirst = cand.fit.tx(zFirst);
    if (std::abs(cand.txFirst) > m_cfg.maxAbsTransverseSlope) {
      return false;
    }

    // Exactly zero curvature would be emitted as an infinite-momentum seed
    // under both charge hypotheses; there is nothing to seed with.
    cand.qOverP = cand.bendFit.qOverP;
    if (!std::isfinite(cand.qOverP) || cand.qOverP == 0.0 ||
        (m_cfg.pMin > 0.0 && std::abs(cand.qOverP) > 1.0 / m_cfg.pMin)) {
      return false;
    }
    cand.inferredCharge = cand.qOverP < 0.0 ? -1 : 1;
    if (cand.bendFit.covarianceValid) {
      cand.qOverPVariance = cand.bendFit.covariance[8];
    } else {
      cand.qOverPVariance = m_cfg.qOverPVariance;
    }
    if (!std::isfinite(cand.qOverPVariance) || cand.qOverPVariance <= 0.0) {
      cand.qOverPVariance = m_cfg.qOverPVariance;
    }
    return true;
  };

  unsigned int truncated        = 0;
  unsigned int endpointPairs    = 0;
  unsigned int endpointRejected = 0;
  for (const auto& stationHits : subsets) {
    // Anchor on the outermost two stations and keep intermediate hits inside
    // a chord road enlarged by the allowed non-bend curvature.
    // This removes most cross-track combinations before any fit, so the
    // budget below binds only in genuinely dense events rather than
    // truncating the enumeration in lexicographic order.
    const std::size_t nStations = stationHits.size();
    const auto& firstStation    = stationHits.front();
    const auto& lastStation     = stationHits.back();

    struct EndpointPair {
      std::size_t first{};
      std::size_t last{};
      double maxAbsCurvatureY{};
      b0stub::EndpointCompatibility compatibility;
    };
    std::vector<EndpointPair> compatibleEndpoints;
    compatibleEndpoints.reserve(firstStation.size() * lastStation.size());
    for (const std::size_t hFirst : firstStation) {
      for (const std::size_t hLast : lastStation) {
        ++endpointPairs;
        const double curvature = std::max(curvatureBounds[hFirst], curvatureBounds[hLast]);
        double zEntrance       = m_z_field_entrance;
        if (m_z_face_lab > 0.0) {
          // Earliest crossing of the lab-z entrance plane allowed by the
          // first hit and |tx| <= maxAbsTransverseSlope. This avoids shrinking
          // the curvature envelope because of the tilted ion coordinates.
          const double tiltSlope = std::abs(sa) * m_cfg.maxAbsTransverseSlope;
          zEntrance =
              (m_z_face_lab + sa * ion[hFirst].x - tiltSlope * ion[hFirst].z) / (ca - tiltSlope);
        }
        const auto compatibility = b0stub::endpointCompatibility(
            ion[hFirst], ion[hLast], m_cfg.maxAbsTransverseSlope, m_cfg.maxYBeamlineResidual,
            m_cfg.constrainToBeamline, curvature, zEntrance);
        if (!compatibility.valid) {
          ++endpointRejected;
          continue;
        }
        compatibleEndpoints.push_back({hFirst, hLast, curvature, compatibility});
      }
    }
    std::ranges::sort(compatibleEndpoints, [&](const EndpointPair& a, const EndpointPair& b) {
      return std::tie(a.compatibility.score, ion[a.first].x, ion[a.first].y, ion[a.last].x,
                      ion[a.last].y) < std::tie(b.compatibility.score, ion[b.first].x,
                                                ion[b.first].y, ion[b.last].x, ion[b.last].y);
    });

    unsigned int tried = 0;
    bool budgetHit     = false;
    for (const auto& endpoints : compatibleEndpoints) {
      if (budgetHit) {
        break;
      }
      const std::size_t hFirst = endpoints.first;
      const std::size_t hLast  = endpoints.last;
      const double zA          = ion[hFirst].z;
      const double slopeY      = endpoints.compatibility.slopeY;

      // Candidate hits per intermediate station, inside the y road.
      std::vector<std::vector<std::size_t>> middle;
      middle.reserve(nStations - 2);
      bool empty = false;
      for (std::size_t s = 1; s + 1 < nStations; ++s) {
        std::vector<std::pair<double, std::size_t>> ranked;
        for (std::size_t idx : stationHits[s]) {
          const double yPred    = ion[hFirst].y + slopeY * (ion[idx].z - zA);
          const double residual = std::abs(ion[idx].y - yPred);
          // A curve with |y''| <= A deviates from its endpoint chord by at
          // most A/2 * (z-zFirst)*(zLast-z) between the endpoints.
          const double allowance =
              0.5 * endpoints.maxAbsCurvatureY * (ion[idx].z - zA) * (ion[hLast].z - ion[idx].z);
          if (residual <= m_cfg.yRoadWidth + allowance) {
            ranked.emplace_back(residual, idx);
          }
        }
        std::ranges::sort(ranked, [&](const auto& a, const auto& b) {
          return std::tie(a.first, ion[a.second].x, ion[a.second].y, ion[a.second].z) <
                 std::tie(b.first, ion[b.second].x, ion[b.second].y, ion[b.second].z);
        });
        std::vector<std::size_t> keep;
        keep.reserve(ranked.size());
        for (const auto& [residual, idx] : ranked) {
          keep.push_back(idx);
        }
        if (keep.empty()) {
          empty = true;
          break;
        }
        middle.push_back(std::move(keep));
      }
      if (empty) {
        continue;
      }

      // Enumerate the surviving middle combinations.
      std::vector<std::size_t> cursor(middle.size(), 0);
      while (true) {
        if (tried >= perSubsetBudget) {
          budgetHit = true;
          ++truncated;
          break;
        }
        std::vector<std::size_t> idxs;
        idxs.reserve(nStations);
        idxs.push_back(hFirst);
        for (std::size_t s = 0; s < middle.size(); ++s) {
          idxs.push_back(middle[s][cursor[s]]);
        }
        idxs.push_back(hLast);

        std::vector<b0stub::Point3> pts;
        pts.reserve(idxs.size());
        for (std::size_t idx : idxs) {
          pts.push_back(ion[idx]);
        }
        StubCandidate cand;
        cand.hitIndices = std::move(idxs);
        cand.fit        = b0stub::fitStub(pts);
        if (makeCompatible(cand)) {
          candidates.push_back(std::move(cand));
        }
        ++tried;

        std::size_t s = 0;
        for (; s < middle.size(); ++s) {
          if (++cursor[s] < middle[s].size()) {
            break;
          }
          cursor[s] = 0;
        }
        if (s == middle.size()) {
          break; // enumerated this pair
        }
      }
    }
  }
  if (truncated > 0) {
    debug("B0 stub seeding rejected {}/{} endpoint pairs and hit the combination budget in {} "
          "of {} station subsets; some combinations were not tried",
          endpointRejected, endpointPairs, truncated, subsets.size());
  }

  if (candidates.empty()) {
    return;
  }

  // ------------------------------------------------------------------
  // Rank by residual, deduplicate by shared hits, emit up to maxSeeds.
  // ------------------------------------------------------------------
  // More stations first (a 3-point auxiliary fit interpolates exactly, so raw
  // residuals would always favour subsets over genuine full-station
  // candidates), then residual within equal station count.
  std::ranges::sort(candidates, [](const StubCandidate& a, const StubCandidate& b) {
    if (a.hitIndices.size() != b.hitIndices.size()) {
      return a.hitIndices.size() > b.hitIndices.size();
    }
    return a.fit.rmsX * a.fit.rmsX + a.fit.rmsY * a.fit.rmsY <
           b.fit.rmsX * b.fit.rmsX + b.fit.rmsY * b.fit.rmsY;
  });

  // Two hits are "the same" for sharing purposes if they are one hit, or if
  // they sit at the same station within sharedHitDistance of each other -
  // the front/back sensor pair of one disk seeing the same track.
  const double dupDist2 = static_cast<double>(m_cfg.sharedHitDistance) * m_cfg.sharedHitDistance;
  const auto sameHit    = [&](std::size_t i, std::size_t j) {
    if (i == j) {
      return true;
    }
    if (dupDist2 <= 0.0 || std::abs(ion[i].z - ion[j].z) > m_cfg.stationZGap) {
      return false;
    }
    const double dx = ion[i].x - ion[j].x;
    const double dy = ion[i].y - ion[j].y;
    return dx * dx + dy * dy <= dupDist2;
  };

  std::vector<const StubCandidate*> accepted;
  for (const auto& cand : candidates) {
    if (accepted.size() >= m_cfg.maxSeeds) {
      break;
    }
    bool overlaps = false;
    for (const auto* acc : accepted) {
      unsigned int shared = 0;
      for (std::size_t i : cand.hitIndices) {
        for (std::size_t j : acc->hitIndices) {
          if (sameHit(i, j)) {
            ++shared;
            break;
          }
        }
      }
      if (shared > m_cfg.maxSharedHits) {
        overlaps = true;
        break;
      }
    }
    if (!overlaps) {
      accepted.push_back(&cand);
    }
  }

  // ------------------------------------------------------------------
  // Convert accepted candidates to origin-perigee seed parameters.
  // ------------------------------------------------------------------
  unsigned int emitted = 0;
  for (const auto* cand : accepted) {
    if (emitted >= m_cfg.maxSeeds) {
      break;
    }
    const SeedVector state =
        seedStateFromFit(cand->bendFit, cand->fit, m_crossing_angle, m_cfg.constrainToBeamline);
    const SeedMatrix fitCovariance = propagateFitCovariance(
        cand->bendFit, cand->fit, m_crossing_angle, m_cfg.constrainToBeamline);

    const std::vector<int> charges =
        b0stub::chargeHypotheses(cand->qOverP, cand->qOverPVariance, m_cfg.minCurvatureSignificance,
                                 cand->inferredCharge, m_cfg.testBothCharges, m_cfg.charge);

    for (const int charge : charges) {
      if (emitted >= m_cfg.maxSeeds) {
        break;
      }
      const double emittedQOverP = charge * std::abs(cand->qOverP);
      const double qOverPSign    = cand->qOverP == 0.0 ? 1.0 : emittedQOverP / cand->qOverP;
      SeedMatrix covariance      = fitCovariance;
      covariance.row(4) *= qOverPSign;
      covariance.col(4) *= qOverPSign;

      // The constrained seed assumes the origin; the interaction-vertex
      // spread is an uncertainty on the parameters it reports.
      const auto beamSpot = b0stub::beamSpotCovarianceAdditions(
          cand->bendFit, cand->fit, m_crossing_angle, m_cfg.constrainToBeamline,
          m_cfg.beamSpotSizeX, m_cfg.beamSpotSizeY, m_cfg.beamSpotSizeZ);
      const std::array<double, 5> fallbackVariances{m_cfg.locaVariance, m_cfg.locbVariance,
                                                    m_cfg.phiVariance, m_cfg.thetaVariance,
                                                    m_cfg.qOverPVariance};
      std::array<double, 25> signedFitCovariance{};
      Eigen::Map<RowMajor5>(signedFitCovariance.data()) = covariance;
      const auto combined =
          b0stub::seedCovarianceWithFallbacks(signedFitCovariance, beamSpot, fallbackVariances);
      covariance            = Eigen::Map<const RowMajor5>(combined.data());
      const auto scattering = b0stub::scatteringCovarianceAdditions(
          emittedQOverP, state(3), m_cfg.scatteringScale, m_cfg.qOverPRelativeUncertainty);
      covariance(2, 2) += scattering[0];
      covariance(3, 3) += scattering[1];
      covariance(4, 4) += scattering[2];

      auto trackparam = track_params_output->create();
      trackparam.setType(-1); // seed
      trackparam.setLoc({static_cast<float>(state(0)), static_cast<float>(state(1))});
      trackparam.setPhi(static_cast<float>(state(2)));
      trackparam.setTheta(static_cast<float>(state(3)));
      trackparam.setQOverP(static_cast<float>(emittedQOverP));
      trackparam.setTime(0.0F); // vertex time at the origin perigee

      edm4eic::Cov6f cov;
      for (int row = 0; row < 5; ++row) {
        for (int col = row; col < 5; ++col) {
          cov(row, col) = static_cast<float>(covariance(row, col));
        }
      }
      cov(5, 5) = m_cfg.timeVariance;
      trackparam.setCovariance(cov);

      auto seed = seeds->create();
      seed.setPerigee({0.F, 0.F, 0.F});
      // Larger is better for ACTS seed quality; use the negative mean residual.
      seed.setQuality(
          static_cast<float>(-(cand->fit.rmsX * cand->fit.rmsX + cand->fit.rmsY * cand->fit.rmsY)));
      seed.setParams(trackparam);
      for (std::size_t idx : cand->hitIndices) {
        seed.addToHits(ionHits[idx]);
      }
      ++emitted;

      trace("B0 stub seed: q={} q/p={:.5f} 1/GeV theta={:.4f} phi={:.3f} "
            "loc=({:.2f},{:.2f}) rms=({:.3f},{:.3f})",
            charge, emittedQOverP, state(3), state(2), state(0), state(1), cand->fit.rmsX,
            cand->fit.rmsY);
    }
  }
}

} // namespace eicrecon
