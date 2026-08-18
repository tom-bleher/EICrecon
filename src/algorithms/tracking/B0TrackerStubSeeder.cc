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
//   2. fit y(z) with a line and use a bend-plane parabola only to initialize
//      sampling of the full field along each candidate;
//   3. fit x(z) in the sampled field-integral basis for signed q/p and local
//      direction, then back-extrapolate and express the result on the origin perigee
//      surface, which is what CKFTracking hard-codes as seed reference.
//
// The seed is only a starting point: the CKF refits with the full field and
// material, so the stub q/p must not be read as a final measurement. Because
// the seed is expressed at the origin it is a prompt-track solution; see
// constrainToBeamline in the config for why displaced decays need more.

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
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "algorithms/interfaces/ActsSvc.h"

namespace eicrecon {

namespace b0stub {

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
    const Eigen::Vector3d ax = decompositionX.solve(weightedValuesX);
    const Eigen::Vector2d ay = decompositionY.solve(weightedValuesY);
    for (int i = 0; i < 3; ++i) {
      fit.basisX[i] = ax(i);
    }
    for (int i = 0; i < 2; ++i) {
      fit.basisY[i] = ay(i);
    }

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
        for (int row = 0; row < 3; ++row) {
          for (int col = 0; col < 3; ++col) {
            fit.covarianceX[3 * row + col] = covAx(row, col);
          }
        }
        for (int row = 0; row < 2; ++row) {
          for (int col = 0; col < 2; ++col) {
            fit.covarianceY[2 * row + col] = covAy(row, col);
          }
        }
        fit.covarianceValid = true;
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

  std::vector<FieldIntegral> integrateFieldSamples(const std::vector<FieldSample>& samples) {
    std::vector<FieldIntegral> result(samples.size());
    if (samples.empty()) {
      return result;
    }
    if (!std::isfinite(samples.front().z) || !std::isfinite(samples.front().fieldY)) {
      return {};
    }
    for (std::size_t i = 1; i < samples.size(); ++i) {
      const double h      = samples[i].z - samples[i - 1].z;
      const double fieldA = samples[i - 1].fieldY;
      const double fieldB = samples[i].fieldY;
      if (!(h > 0.0) || !std::isfinite(fieldA) || !std::isfinite(fieldB)) {
        return {};
      }
      result[i].second =
          result[i - 1].second + h * result[i - 1].first + h * h * (fieldA / 3.0 + fieldB / 6.0);
      result[i].first = result[i - 1].first + 0.5 * h * (fieldA + fieldB);
    }
    return result;
  }

  FieldIntegralFit fitFieldIntegral(const std::vector<Point3>& pts,
                                    const std::vector<double>& secondIntegrals, double zReference) {
    FieldIntegralFit fit;
    fit.zReference      = zReference;
    const std::size_t n = pts.size();
    if (n < 3 || secondIntegrals.size() != n || !std::isfinite(zReference)) {
      return fit;
    }

    double zScale          = 0.0;
    double bendScale       = 0.0;
    constexpr double kBend = 2.998e-4;
    for (std::size_t i = 0; i < n; ++i) {
      zScale    = std::max(zScale, std::abs(pts[i].z - zReference));
      bendScale = std::max(bendScale, std::abs(kBend * secondIntegrals[i]));
    }
    if (!(zScale > 0.0) || !(bendScale > 0.0)) {
      return fit;
    }

    Eigen::MatrixXd design(n, 3);
    Eigen::VectorXd values(n);
    bool haveVariances = true;
    for (std::size_t row = 0; row < n; ++row) {
      if (!(std::isfinite(pts[row].x) && std::isfinite(pts[row].z) &&
            std::isfinite(secondIntegrals[row]))) {
        return fit;
      }
      design(row, 0) = 1.0;
      design(row, 1) = (pts[row].z - zReference) / zScale;
      design(row, 2) = -kBend * secondIntegrals[row] / bendScale;
      values(row)    = pts[row].x;
      haveVariances &= std::isfinite(pts[row].varianceX) && pts[row].varianceX > 0.0;
    }

    Eigen::MatrixXd weightedDesign = design;
    Eigen::VectorXd weightedValues = values;
    if (haveVariances) {
      for (std::size_t row = 0; row < n; ++row) {
        const double weight = 1.0 / std::sqrt(pts[row].varianceX);
        weightedDesign.row(row) *= weight;
        weightedValues(row) *= weight;
      }
    }

    const auto decomposition = weightedDesign.colPivHouseholderQr();
    if (decomposition.rank() < 3) {
      return fit;
    }
    const Eigen::Vector3d beta = decomposition.solve(weightedValues);
    fit.xReference             = beta(0);
    fit.txReference            = beta(1) / zScale;
    fit.qOverP                 = beta(2) / bendScale;

    double sumResidual2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double predicted = fit.xReference + fit.txReference * (pts[i].z - zReference) -
                               kBend * fit.qOverP * secondIntegrals[i];
      const double residual = pts[i].x - predicted;
      sumResidual2 += residual * residual;
    }
    fit.rmsX = std::sqrt(sumResidual2 / static_cast<double>(n));

    if (haveVariances) {
      const Eigen::Matrix3d normal           = weightedDesign.transpose() * weightedDesign;
      const Eigen::Matrix3d scaledCovariance = normal.ldlt().solve(Eigen::Matrix3d::Identity());
      Eigen::Matrix3d transform              = Eigen::Matrix3d::Identity();
      transform(1, 1)                        = 1.0 / zScale;
      transform(2, 2)                        = 1.0 / bendScale;
      const Eigen::Matrix3d covariance       = transform * scaledCovariance * transform.transpose();
      if (covariance.allFinite() && (covariance.diagonal().array() > 0.0).all()) {
        for (int row = 0; row < 3; ++row) {
          for (int col = 0; col < 3; ++col) {
            fit.covariance[3 * row + col] = covariance(row, col);
          }
        }
        fit.covarianceValid = true;
      }
    }

    fit.valid = std::isfinite(fit.xReference) && std::isfinite(fit.txReference) &&
                std::isfinite(fit.qOverP) && std::isfinite(fit.rmsX);
    return fit;
  }

  EndpointCompatibility endpointCompatibility(const Point3& first, const Point3& last,
                                              double maxAbsSlope, double maxBeamResidual,
                                              bool constrainToBeamline) {
    EndpointCompatibility result;
    const double dz = last.z - first.z;
    if (!(dz > 0.0) || !std::isfinite(first.y) || !std::isfinite(first.z) ||
        !std::isfinite(last.y) || !std::isfinite(last.z) || !(maxAbsSlope > 0.0)) {
      return result;
    }
    result.slopeY = (last.y - first.y) / dz;
    if (!std::isfinite(result.slopeY) || std::abs(result.slopeY) > maxAbsSlope) {
      return result;
    }

    if (constrainToBeamline && maxBeamResidual > 0.0) {
      if (first.z == 0.0) {
        return result;
      }
      result.beamResidual = last.y - first.y * last.z / first.z;
      if (!std::isfinite(result.beamResidual) || std::abs(result.beamResidual) > maxBeamResidual) {
        return result;
      }
    }
    const double slopeScore = std::abs(result.slopeY) / maxAbsSlope;
    const double beamScore  = constrainToBeamline && maxBeamResidual > 0.0
                                  ? std::abs(result.beamResidual) / maxBeamResidual
                                  : 0.0;
    result.score            = beamScore + slopeScore;
    result.valid            = true;
    return result;
  }

  std::array<double, 3> covarianceModelAdditions(double qOverP, double phiModelVariance,
                                                 double phiQOverPScale, double thetaModelVariance,
                                                 double thetaQOverPScale,
                                                 double qOverPModelVariance,
                                                 double qOverPRelativeUncertainty,
                                                 double fieldRelativeUncertainty) {
    const double qOverP2 = qOverP * qOverP;
    return {
        std::max(phiModelVariance, 0.0) + phiQOverPScale * phiQOverPScale * qOverP2,
        std::max(thetaModelVariance, 0.0) + thetaQOverPScale * thetaQOverPScale * qOverP2,
        std::max(qOverPModelVariance, 0.0) +
            (qOverPRelativeUncertainty * qOverPRelativeUncertainty +
             fieldRelativeUncertainty * fieldRelativeUncertainty) *
                qOverP2,
    };
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
    std::sort(ordered.begin(), ordered.end());

    StationInterval current{ordered.front(), ordered.front(), ordered.front()};
    double sum     = ordered.front();
    std::size_t n  = 1;
    for (std::size_t i = 1; i < ordered.size(); ++i) {
      if (ordered[i] - current.zMax > gap) {
        current.zMean = sum / static_cast<double>(n);
        stations.push_back(current);
        current = {ordered[i], ordered[i], ordered[i]};
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
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return hitZ[a] < hitZ[b]; });
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
    b0stub::FieldIntegralFit bendFit;
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

  EntranceFit entranceFit(const b0stub::FieldIntegralFit& bendFit,
                          const b0stub::StubFit& nonBendFit) {
    const double u = (bendFit.zReference - nonBendFit.zRef) / nonBendFit.zScale;
    return {bendFit.xReference, bendFit.txReference, bendFit.qOverP,
            nonBendFit.basisY[0] + nonBendFit.basisY[1] * u,
            nonBendFit.basisY[1] / nonBendFit.zScale};
  }

  SeedVector seedStateFromEntrance(const EntranceFit& entrance, double zFieldEntrance,
                                   double crossingAngle, bool constrainToBeamline) {
    const double ca = std::cos(crossingAngle);
    const double sa = std::sin(crossingAngle);

    double txi;
    double tyi;
    double x0;
    double y0;
    if (constrainToBeamline) {
      txi = entrance.x / zFieldEntrance;
      tyi = entrance.y / zFieldEntrance;
      x0  = 0.0;
      y0  = 0.0;
    } else {
      txi = entrance.tx;
      tyi = entrance.ty;
      x0  = entrance.x - entrance.tx * zFieldEntrance;
      y0  = entrance.y - entrance.ty * zFieldEntrance;
    }

    const b0stub::Point3 dir{txi * ca + sa, tyi, -txi * sa + ca};
    const b0stub::Point3 ref{x0 * ca, y0, -x0 * sa};
    const auto perigee = b0stub::perigeeFromRay(ref, dir, {0.0, 0.0, 0.0});

    SeedVector state;
    state << perigee.loc0, perigee.loc1, perigee.phi, perigee.theta, entrance.qOverP;
    return state;
  }

  SeedVector seedStateFromFit(const b0stub::FieldIntegralFit& bendFit,
                              const b0stub::StubFit& nonBendFit, double crossingAngle,
                              bool constrainToBeamline) {
    return seedStateFromEntrance(entranceFit(bendFit, nonBendFit), bendFit.zReference,
                                 crossingAngle, constrainToBeamline);
  }

  SeedMatrix propagateFitCovariance(const b0stub::FieldIntegralFit& bendFit,
                                    const b0stub::StubFit& nonBendFit, double crossingAngle,
                                    bool constrainToBeamline) {
    SeedMatrix result = SeedMatrix::Zero();
    if (!bendFit.covarianceValid || !nonBendFit.covarianceValid) {
      return result;
    }

    SeedMatrix fitCovariance = SeedMatrix::Zero();
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        fitCovariance(row, col) = bendFit.covariance[3 * row + col];
      }
    }
    Eigen::Matrix2d yCovariance;
    for (int row = 0; row < 2; ++row) {
      for (int col = 0; col < 2; ++col) {
        yCovariance(row, col) = nonBendFit.covarianceY[2 * row + col];
      }
    }
    const double uEntrance = (bendFit.zReference - nonBendFit.zRef) / nonBendFit.zScale;
    Eigen::Matrix2d yTransform;
    yTransform << 1.0, uEntrance, 0.0, 1.0 / nonBendFit.zScale;
    fitCovariance.block<2, 2>(3, 3) = yTransform * yCovariance * yTransform.transpose();

    const EntranceFit nominal = entranceFit(bendFit, nonBendFit);
    const std::array<double, 5> values{nominal.x, nominal.tx, nominal.qOverP, nominal.y,
                                       nominal.ty};
    SeedMatrix jacobian = SeedMatrix::Zero();
    for (int column = 0; column < 5; ++column) {
      const double sigma = std::sqrt(std::max(0.0, fitCovariance(column, column)));
      const double step  = std::max({1.0e-10, sigma * 1.0e-3, std::abs(values[column]) * 1.0e-7});
      auto plus          = nominal;
      auto minus         = nominal;
      double* plusValue  = nullptr;
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

  std::array<double, 25> seedCovarianceFromFit(const FieldIntegralFit& bendFit,
                                               const StubFit& nonBendFit, double crossingAngle,
                                               bool constrainToBeamline) {
    const SeedMatrix covariance =
        propagateFitCovariance(bendFit, nonBendFit, crossingAngle, constrainToBeamline);
    std::array<double, 25> result{};
    for (int row = 0; row < 5; ++row) {
      for (int col = 0; col < 5; ++col) {
        result[5 * row + col] = covariance(row, col);
      }
    }
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
    const bool measurable = std::isfinite(qOverP) && qOverP != 0.0 &&
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
  m_crossing_angle = detector->constant<double>("CrossingAngle") / dd4hep::rad;
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
      m_z_field_entrance = xLab * sa + zLab * ca;
    } catch (const std::exception& error) {
      throw std::runtime_error(
          std::string("B0TrackerStubSeeder: cannot derive zFieldEntrance from B0PF "
                      "geometry; set B0TrackerStubSeeder:zFieldEntrance or provide "
                      "B0PF_CenterPosition and B0PF_Length (") +
          error.what() + ")");
    }
  }
  if (!std::isfinite(m_z_field_entrance)) {
    throw std::runtime_error("B0TrackerStubSeeder: B0pf field entrance z is not finite");
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

  if (hits->empty()) {
    return;
  }

  const double ca = std::cos(m_crossing_angle);
  const double sa = std::sin(m_crossing_angle);

  // ------------------------------------------------------------------
  // Rotate hit positions into the ion frame and group them by station.
  // ------------------------------------------------------------------
  std::vector<b0stub::Point3> ion;
  ion.reserve(hits->size());
  std::vector<edm4eic::TrackerHit> ionHits;
  ionHits.reserve(hits->size());

  for (const auto& hit : *hits) {
    const auto& p          = hit.getPosition();
    const auto& covariance = hit.getPositionError();
    // TrackerHit position covariance is expressed in the sensor surface frame.
    // B0's CartesianGridXY sensors rotate with the assembly, so their local
    // x/y axes are precisely the ion-frame bend/non-bend axes used by this fit.
    const b0stub::Point3 q{p.x * ca - p.z * sa, p.y, p.x * sa + p.z * ca, covariance.xx,
                           covariance.yy};
    if (!(std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z))) {
      debug("Skipping B0 hit with non-finite position");
      continue;
    }
    ion.push_back(q);
    ionHits.push_back(hit);
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
  const auto* converter =
      m_volume_to_station.empty() ? nullptr : algorithms::GeoSvc::instance().cellIDPositionConverter();
  if (!m_stations.empty()) {
    for (std::size_t i = 0; i < ion.size(); ++i) {
      int station = -1;
      if (converter != nullptr) {
        try {
          const auto* volumeContext = converter->findContext(ionHits[i].getCellID());
          if (volumeContext != nullptr) {
            const auto found = m_volume_to_station.find(volumeContext->identifier);
            if (found != m_volume_to_station.end()) {
              station = static_cast<int>(found->second);
            }
          }
        } catch (const std::exception&) {
          station = -1;
        }
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
    std::sort(indices.begin(), indices.end(), [&](std::size_t a, std::size_t b) {
      return std::tie(ion[a].x, ion[a].y, ion[a].z) < std::tie(ion[b].x, ion[b].y, ion[b].z);
    });
  }

  // The auxiliary and field-integral fits need >= 3 points, so 3 stations is the hard floor
  // regardless of configuration (a 2-station mode would need a line fit with
  // a momentum prior instead).
  const unsigned int minStations = std::max(m_cfg.minStations, 3u);
  if (byStation.size() < minStations) {
    return;
  }

  // Sample the same ACTS/DD4hep field provider used by CKFTracking. The B0
  // field has a dipole component and a gradient, so an analytic nominal B is
  // not portable across beam-energy configurations or hit positions.
  const auto fieldProvider = m_acts_context->getFieldProvider();
  auto fieldCache = fieldProvider->makeCache(m_acts_context->getActsMagneticFieldContext());

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
      std::max(1u, m_cfg.maxCombinations / static_cast<unsigned int>(subsets.size()));

  auto makeCompatible = [&](StubCandidate& cand) -> bool {
    if (!cand.fit.valid || cand.fit.rmsX > m_cfg.maxXResidual ||
        cand.fit.rmsY > m_cfg.maxYResidual) {
      return false;
    }

    const auto zMinMax =
        std::minmax_element(cand.hitIndices.begin(), cand.hitIndices.end(),
                            [&](std::size_t a, std::size_t b) { return ion[a].z < ion[b].z; });
    const double zFirst = ion[*zMinMax.first].z;
    const double zLast  = ion[*zMinMax.second].z;
    if (std::abs(cand.fit.b1) > m_cfg.maxAbsTransverseSlope || !(m_z_field_entrance < zFirst)) {
      return false;
    }

    std::vector<b0stub::Point3> points;
    points.reserve(cand.hitIndices.size());
    std::vector<double> meshZ{m_z_field_entrance};
    for (const std::size_t idx : cand.hitIndices) {
      points.push_back(ion[idx]);
      meshZ.push_back(ion[idx].z);
    }
    const unsigned int requestedSamples = std::max(2u, m_cfg.fieldSamples);
    for (unsigned int i = 1; i + 1 < requestedSamples; ++i) {
      const double fraction = static_cast<double>(i) / (requestedSamples - 1);
      meshZ.push_back(m_z_field_entrance + fraction * (zLast - m_z_field_entrance));
    }
    std::sort(meshZ.begin(), meshZ.end());
    meshZ.erase(std::unique(meshZ.begin(), meshZ.end()), meshZ.end());

    const auto meshIndex = [&](double z) -> std::size_t {
      const auto found = std::lower_bound(meshZ.begin(), meshZ.end(), z);
      if (found == meshZ.end() || std::abs(*found - z) > 1.0e-8) {
        return meshZ.size();
      }
      return static_cast<std::size_t>(std::distance(meshZ.begin(), found));
    };

    b0stub::FieldIntegralFit previousFit;
    std::vector<b0stub::FieldIntegral> previousMoments;
    const unsigned int iterations = std::max(1u, m_cfg.fieldFitIterations);
    for (unsigned int iteration = 0; iteration < iterations; ++iteration) {
      std::vector<b0stub::FieldSample> samples;
      samples.reserve(meshZ.size());
      double maxAbsFieldY = 0.0;
      for (std::size_t i = 0; i < meshZ.size(); ++i) {
        const double z = meshZ[i];
        double x       = cand.fit.x(z);
        if (iteration > 0) {
          constexpr double kBend = 2.998e-4;
          x = previousFit.xReference + previousFit.txReference * (z - m_z_field_entrance) -
              kBend * previousFit.qOverP * previousMoments[i].second;
        }
        const double y = cand.fit.y(z);
        const Acts::Vector3 global{x * ca + z * sa, y, -x * sa + z * ca};
        const auto field = fieldProvider->getField(global, fieldCache);
        if (!field.ok() || !field.value().allFinite()) {
          return false;
        }
        const double fieldY = field.value().y() / Acts::UnitConstants::T;
        maxAbsFieldY        = std::max(maxAbsFieldY, std::abs(fieldY));
        samples.push_back({z, fieldY});
      }
      if (maxAbsFieldY < m_cfg.minAbsFieldY) {
        return false;
      }

      const auto moments = b0stub::integrateFieldSamples(samples);
      if (moments.size() != meshZ.size()) {
        return false;
      }
      std::vector<double> hitIntegrals;
      hitIntegrals.reserve(points.size());
      for (const auto& point : points) {
        const std::size_t index = meshIndex(point.z);
        if (index == meshZ.size()) {
          return false;
        }
        hitIntegrals.push_back(moments[index].second);
      }
      auto fieldFit = b0stub::fitFieldIntegral(points, hitIntegrals, m_z_field_entrance);
      if (!fieldFit.valid) {
        return false;
      }
      previousFit     = fieldFit;
      previousMoments = moments;
    }

    cand.bendFit  = previousFit;
    cand.fit.rmsX = cand.bendFit.rmsX;
    if (cand.fit.rmsX > m_cfg.maxXResidual) {
      return false;
    }
    const std::size_t firstIndex = meshIndex(zFirst);
    if (firstIndex == meshZ.size()) {
      return false;
    }
    constexpr double kBend = 2.998e-4;
    cand.txFirst =
        cand.bendFit.txReference - kBend * cand.bendFit.qOverP * previousMoments[firstIndex].first;
    if (std::abs(cand.txFirst) > m_cfg.maxAbsTransverseSlope) {
      return false;
    }

    cand.qOverP = cand.bendFit.qOverP;
    if (!std::isfinite(cand.qOverP) ||
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
    // The dipole does not bend in y, so the hits of a single track lie on a
    // straight line in y(z). Anchor on the outermost two stations and keep
    // only the intermediate hits that sit inside a road around that chord.
    // This removes most cross-track combinations before any fit, so the
    // budget below binds only in genuinely dense events rather than
    // truncating the enumeration in lexicographic order.
    const std::size_t nStations = stationHits.size();
    const auto& firstStation    = stationHits.front();
    const auto& lastStation     = stationHits.back();

    struct EndpointPair {
      std::size_t first{};
      std::size_t last{};
      b0stub::EndpointCompatibility compatibility;
    };
    std::vector<EndpointPair> compatibleEndpoints;
    compatibleEndpoints.reserve(firstStation.size() * lastStation.size());
    for (const std::size_t hFirst : firstStation) {
      for (const std::size_t hLast : lastStation) {
        ++endpointPairs;
        const auto compatibility =
            b0stub::endpointCompatibility(ion[hFirst], ion[hLast], m_cfg.maxAbsTransverseSlope,
                                          m_cfg.maxYBeamlineResidual, m_cfg.constrainToBeamline);
        if (!compatibility.valid) {
          ++endpointRejected;
          continue;
        }
        compatibleEndpoints.push_back({hFirst, hLast, compatibility});
      }
    }
    std::sort(compatibleEndpoints.begin(), compatibleEndpoints.end(),
              [&](const EndpointPair& a, const EndpointPair& b) {
                return std::tie(a.compatibility.score, ion[a.first].x, ion[a.first].y,
                                ion[a.last].x, ion[a.last].y) <
                       std::tie(b.compatibility.score, ion[b.first].x, ion[b.first].y,
                                ion[b.last].x, ion[b.last].y);
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
      const double zB          = ion[hLast].z;
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
          if (residual <= m_cfg.yRoadWidth) {
            ranked.emplace_back(residual, idx);
          }
        }
        std::sort(ranked.begin(), ranked.end(), [&](const auto& a, const auto& b) {
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
  std::sort(candidates.begin(), candidates.end(),
            [](const StubCandidate& a, const StubCandidate& b) {
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

      const std::array<double, 5> fallbackVariances{m_cfg.locaVariance, m_cfg.locbVariance,
                                                    m_cfg.phiVariance, m_cfg.thetaVariance,
                                                    m_cfg.qOverPVariance};
      for (int i = 0; i < 5; ++i) {
        if (!std::isfinite(covariance(i, i)) || covariance(i, i) <= 0.0) {
          covariance.row(i).setZero();
          covariance.col(i).setZero();
          covariance(i, i) = fallbackVariances[i];
        }
      }
      const auto modelVariances = b0stub::covarianceModelAdditions(
          emittedQOverP, m_cfg.phiModelVariance, m_cfg.phiQOverPScale, m_cfg.thetaModelVariance,
          m_cfg.thetaQOverPScale, m_cfg.qOverPModelVariance, m_cfg.qOverPRelativeUncertainty,
          m_cfg.fieldRelativeUncertainty);
      covariance(2, 2) += modelVariances[0];
      covariance(3, 3) += modelVariances[1];
      covariance(4, 4) += modelVariances[2];

      auto trackparam = track_params_output->create();
      trackparam.setType(-1); // seed
      trackparam.setLoc({static_cast<float>(state(0)), static_cast<float>(state(1))});
      trackparam.setPhi(static_cast<float>(state(2)));
      trackparam.setTheta(static_cast<float>(state(3)));
      trackparam.setQOverP(static_cast<float>(emittedQOverP));
      trackparam.setTime(10);

      edm4eic::Cov6f cov;
      for (int row = 0; row < 5; ++row) {
        for (int col = row; col < 5; ++col) {
          cov(row, col) = static_cast<float>(covariance(row, col));
        }
      }
      cov(5, 5) = m_cfg.timeVariance;
      trackparam.setCovariance(cov);

      auto seed = seeds->create();
      seed.setPerigee({0.f, 0.f, 0.f});
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
