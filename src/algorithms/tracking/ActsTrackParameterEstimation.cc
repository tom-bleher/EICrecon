// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Extracted from ACTS v47.7.0, commit
// 2790f1b05c0c94cb3b569cbb18262b24867c7855:
// Core/src/Seeding/EstimateTrackParamsFromSeed.cpp and
// Core/include/Acts/Seeding/detail/CircleFit.hpp.
// Upstream source SHA256:
// 1cf302f132a0d9ce7837e4649c8bc46c22c07b3690b3feb4fdc4f0c450265d16
// Upstream CircleFit.hpp SHA256:
// 88841acf4d6d95bc86dbf46a24b0a48dda7c1a1d0f8cfd398d91d4672d22fa2d
// Changes: private namespace, optional result instead of ACTS error enum,
// and native-API dispatch. Fit arithmetic is unchanged.
#include "ActsTrackParameterEstimation.h"
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Utilities/MathHelpers.hpp>
#include <Eigen/Eigenvalues>
#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#if EICRECON_HAS_ACTS_MULTIPOINT
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#else
namespace {
using namespace Acts;
namespace detail {

  /// A circle in a plane, as fitted to a set of points.
  struct CircleFit {
    /// Circle center.
    Vector2 center = Vector2::Zero();
    /// Circle radius.
    double radius = 0;
  };

  /// Algebraic circle fit (Taubin) over the transverse `(x, y)` projection.
  ///
  /// Fits `A(x^2+y^2) + B*x + C*y + D = 0` with `A` free to reach zero, so
  /// near-collinear input degrades to a line, reported as `std::nullopt`,
  /// instead of becoming ill-conditioned. Weights are relative factors on the
  /// residuals; an empty span means uniform.
  ///
  /// @param points the points whose transverse projection is fitted
  /// @param weights optional per-point weights (empty span = uniform)
  /// @return the fitted circle, or `std::nullopt` if no finite circle is defined
  inline std::optional<CircleFit> fitCircleTaubin(std::span<const Vector3> points,
                                                  std::span<const double> weights = {}) {
    const std::size_t n = points.size();
    if (n < 3) {
      return std::nullopt;
    }
    assert((weights.empty() || weights.size() == n) &&
           "weights must be empty or match the number of points");

    const auto weightAt = [&](std::size_t i) { return weights.empty() ? 1 : weights[i]; };

    // Weighted centroid of the transverse projection.
    double sumW  = 0;
    double meanX = 0;
    double meanY = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const double w = weightAt(i);
      sumW += w;
      meanX += w * points[i].x();
      meanY += w * points[i].y();
    }
    if (sumW <= 0) {
      return std::nullopt;
    }
    const double invW = 1 / sumW;
    meanX *= invW;
    meanY *= invW;

    // Centroid-relative weighted moments, z = x^2 + y^2.
    double mxx = 0;
    double myy = 0;
    double mxy = 0;
    double mxz = 0;
    double myz = 0;
    double mzz = 0;
    double mz  = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const double w = weightAt(i);
      const double x = points[i].x() - meanX;
      const double y = points[i].y() - meanY;
      const double z = x * x + y * y;
      mxx += w * x * x;
      myy += w * y * y;
      mxy += w * x * y;
      mxz += w * x * z;
      myz += w * y * z;
      mzz += w * z * z;
      mz += w * z;
    }
    mxx *= invW;
    myy *= invW;
    mxy *= invW;
    mxz *= invW;
    myz *= invW;
    mzz *= invW;
    mz *= invW;
    if (mz <= 0) {
      // All points coincide.
      return std::nullopt;
    }

    // Taubin eigenproblem `M a = lambda N a` for a = (A, B, C), with D = -A*mz.
    SquareMatrix3 m;
    m << mzz - mz * mz, mxz, myz, //
        mxz, mxx, mxy,            //
        myz, mxy, myy;
    SquareMatrix3 nMat = SquareMatrix3::Zero();
    nMat(0, 0)         = 4 * mz;
    nMat(1, 1)         = 1;
    nMat(2, 2)         = 1;

    Eigen::GeneralizedSelfAdjointEigenSolver<SquareMatrix3> solver(m, nMat);
    if (solver.info() != Eigen::Success) {
      return std::nullopt;
    }
    // Ascending eigenvalues, so the first eigenvector minimizes the Taubin cost.
    const Vector3 a     = solver.eigenvectors().col(0);
    const double coeffA = a(0);
    const double coeffB = a(1);
    const double coeffC = a(2);
    if (coeffA == 0) {
      // Perfectly straight: no finite circle.
      return std::nullopt;
    }

    const double coeffD = -coeffA * mz;
    const double r2 =
        (coeffB * coeffB + coeffC * coeffC - 4 * coeffA * coeffD) / (4 * coeffA * coeffA);
    if (!std::isfinite(r2) || r2 <= 0) {
      return std::nullopt;
    }
    const double radius = std::sqrt(r2);

    // A radius far beyond the point spread sqrt(mz) is a line.
    if (constexpr double maxRadiusToSpread = 1e6; radius > maxRadiusToSpread * std::sqrt(mz)) {
      return std::nullopt;
    }

    const Vector2 centerRel(-coeffB / (2 * coeffA), -coeffC / (2 * coeffA));
    return CircleFit{.center = centerRel + Vector2(meanX, meanY), .radius = radius};
  }

  /// Refine a circle fit by minimizing the radial residuals
  /// `sum_i (|p_i - center| - R)^2` with Gauss-Newton steps, on the transverse
  /// `(x, y)` projection. Weights as in @ref fitCircleTaubin.
  ///
  /// @param fit the circle fit to refine, e.g. from @ref fitCircleTaubin
  /// @param points the points whose transverse projection is fitted
  /// @param iterations the maximum number of Gauss-Newton iterations
  /// @param weights optional per-point weights (empty span = uniform)
  /// @return the refined circle, or `std::nullopt` if a step drives the radius
  ///         non-positive, i.e. the refinement collapses the circle
  inline std::optional<CircleFit> refineCircleGeometric(CircleFit fit,
                                                        std::span<const Vector3> points,
                                                        const std::size_t iterations,
                                                        std::span<const double> weights = {}) {
    assert((weights.empty() || weights.size() == points.size()) &&
           "weights must be empty or match the number of points");

    const auto weightAt = [&](std::size_t i) { return weights.empty() ? 1 : weights[i]; };

    for (std::size_t it = 0; it < iterations; ++it) {
      SquareMatrix3 jtj = SquareMatrix3::Zero();
      Vector3 jtr       = Vector3::Zero();
      for (std::size_t i = 0; i < points.size(); ++i) {
        const Vector2 d   = points[i].head<2>() - fit.center;
        const double dist = d.norm();
        if (dist < std::numeric_limits<double>::epsilon()) {
          continue;
        }
        const double w        = weightAt(i);
        const double residual = dist - fit.radius;
        // Jacobian of the residual w.r.t. (cx, cy, R).
        const Vector3 j(-d.x() / dist, -d.y() / dist, -1);
        jtj += w * j * j.transpose();
        jtr += w * j * residual;
      }
      const Vector3 delta = jtj.ldlt().solve(-jtr);
      if (!delta.allFinite()) {
        // Singular normal equations: keep the fit as it is.
        break;
      }
      fit.center.x() += delta.x();
      fit.center.y() += delta.y();
      fit.radius += delta.z();
      if (fit.radius <= 0) {
        return std::nullopt;
      }
      if (delta.norm() < 1e-9) {
        break;
      }
    }
    return fit;
  }

} // namespace detail

Transform3 estimationFrameLocalToGlobal(const Vector3& sp0, const Vector3& sp1,
                                        const Vector3& bField) {
  // Define a new coordinate frame with its origin at the bottom space point, z
  // axis long the magnetic field direction and y axis perpendicular to vector
  // from the bottom to middle space point. Hence, the projection of the middle
  // space point on the transverse plane will be located at the x axis of the
  // new frame.
  //
  // A vanishing field falls back to global z, a reference along the field to
  // any orthogonal y.
  const Vector3 relVec   = sp1 - sp0;
  const double bMag      = bField.norm();
  const Vector3 newZAxis = (bMag > std::numeric_limits<double>::epsilon())
                               ? Vector3(bField / bMag)
                               : Vector3(Vector3::UnitZ());
  Vector3 newYAxis       = newZAxis.cross(relVec);
  if (newYAxis.norm() < std::numeric_limits<double>::epsilon()) {
    newYAxis = newZAxis.unitOrthogonal();
  }
  newYAxis.normalize();
  const Vector3 newXAxis = newYAxis.cross(newZAxis);
  RotationMatrix3 rotation;
  rotation.col(0) = newXAxis;
  rotation.col(1) = newYAxis;
  rotation.col(2) = newZAxis;
  // The center of the new frame is at the bottom space point
  const Translation3 translation(sp0);
  // The transform which constructs the new frame
  return translation * rotation;
}

std::optional<Acts::FreeVector> fallbackEstimate(std::span<const Vector3> spacePoints,
                                                 const Vector3& bField, double t0,
                                                 std::size_t geometricRefineIterations,
                                                 std::span<const double> weights,
                                                 std::size_t referenceIndex) {
  if (spacePoints.size() < 3) {
    return std::nullopt;
  }
  assert((weights.empty() || weights.size() == spacePoints.size()) &&
         "weights must be empty or match the space points");
  assert(referenceIndex < spacePoints.size() && "reference index must point into the space points");

  const auto w = [&](std::size_t i) { return weights.empty() ? 1 : weights[i]; };

  const Vector3& reference = spacePoints[referenceIndex];

  // Estimation frame: field along local +z, fixed by the first two points. It
  // is anchored to the seed rather than to the reference point, so that the
  // fitted helix does not depend on where the parameters are reported.
  const Transform3 transform =
      estimationFrameLocalToGlobal(spacePoints.front(), spacePoints[1], bField);
  const Transform3 toLocal = transform.inverse();
  std::vector<Vector3> local;
  local.reserve(spacePoints.size());
  for (const Vector3& sp : spacePoints) {
    local.push_back(toLocal * sp);
  }

  FreeVector params            = FreeVector::Zero();
  params.segment<3>(eFreePos0) = reference;
  params[eFreeTime]            = t0;

  const double bMag = bField.norm();

  std::optional<detail::CircleFit> circle = detail::fitCircleTaubin(local, weights);
  if (circle.has_value() && geometricRefineIterations > 0) {
    circle = detail::refineCircleGeometric(*circle, local, geometricRefineIterations, weights);
  }

  if (!circle.has_value()) {
    // Straight-line limit: the direction is the principal axis, q/p stays zero.
    double sumW  = 0;
    Vector3 mean = Vector3::Zero();
    for (std::size_t i = 0; i < local.size(); ++i) {
      sumW += w(i);
      mean += w(i) * local[i];
    }
    if (sumW <= 0) {
      return std::nullopt;
    }
    mean /= sumW;
    SquareMatrix3 cov = SquareMatrix3::Zero();
    for (std::size_t i = 0; i < local.size(); ++i) {
      const Vector3 d = local[i] - mean;
      cov += w(i) * d * d.transpose();
    }
    Eigen::SelfAdjointEigenSolver<SquareMatrix3> solver(cov);
    // A non-positive largest eigenvalue means the points coincide.
    if (solver.info() != Eigen::Success || (solver.eigenvalues()[2] <= 0)) {
      return std::nullopt;
    }
    Vector3 localDir = solver.eigenvectors().col(2);
    if (localDir.dot(local.back() - local.front()) < 0) {
      localDir = -localDir;
    }
    localDir.normalize();
    params.segment<3>(eFreeDir0) = transform.linear() * localDir;
    return params;
  }

  const Vector2 center = circle->center;
  const double radius  = circle->radius;

  const Vector2 firstXy = local.front().head<2>();
  const Vector2 refXy   = local[referenceIndex].head<2>();

  // Sense of travel around the circle, +1 for counterclockwise. Taken from the
  // first two points so that it does not depend on the reference point, which
  // may be the last one and so have no successor.
  const Vector2 firstRadial = firstXy - center;
  const Vector2 firstCcwTangent(-firstRadial.y(), firstRadial.x());
  const double rotSense = (firstCcwTangent.dot(local[1].head<2>() - firstXy) >= 0) ? 1 : -1;

  // Tangent at the reference: its radial vector rotated by +90 degrees,
  // oriented along travel.
  const Vector2 radial  = refXy - center;
  const Vector2 tangent = (rotSense * Vector2(-radial.y(), radial.x())).normalized();

  // q v x B points to the center; with B along local +z, v x B = (v_y, -v_x).
  const double toCenterProj =
      tangent.y() * (center.x() - refXy.x()) - tangent.x() * (center.y() - refXy.y());
  const double qSign = (toCenterProj >= 0) ? 1 : -1;

  // z is linear in the arc length s = radius * turning angle. Consecutive
  // points are assumed less than half a turn apart. The slope is invariant
  // under a shift of s, so s is measured from the first point rather than from
  // the reference, which keeps the unwrapping sequential for any reference.
  const double phiFirst = std::atan2(firstXy.y() - center.y(), firstXy.x() - center.x());
  double phiPrev        = phiFirst;
  double sumW           = 0;
  double sumS           = 0;
  double sumZ           = 0;
  double sumSS          = 0;
  double sumSZ          = 0;
  for (std::size_t i = 0; i < local.size(); ++i) {
    const Vector3& p = local[i];
    double phi       = std::atan2(p.y() - center.y(), p.x() - center.x());
    while (phi - phiPrev > std::numbers::pi) {
      phi -= 2 * std::numbers::pi;
    }
    while (phi - phiPrev < -std::numbers::pi) {
      phi += 2 * std::numbers::pi;
    }
    phiPrev         = phi;
    const double wi = w(i);
    const double s  = rotSense * radius * (phi - phiFirst);
    const double z  = p.z();
    sumW += wi;
    sumS += wi * s;
    sumZ += wi * z;
    sumSS += wi * s * s;
    sumSZ += wi * s * z;
  }
  const double denom  = sumW * sumSS - sumS * sumS;
  const double lambda = (std::abs(denom) > std::numeric_limits<double>::epsilon())
                            ? (sumW * sumSZ - sumS * sumZ) / denom
                            : 0;

  // s is the transverse arc length, so the local direction is (t_x, t_y,
  // lambda).
  Vector3 localDir(tangent.x(), tangent.y(), lambda);
  localDir.normalize();
  params.segment<3>(eFreeDir0) = transform.linear() * localDir;

  // p_T = |B| * R, so q/p_T = qSign / (R * |B|). Without a field q/p stays
  // zero.
  if (bMag > std::numeric_limits<double>::epsilon()) {
    const double qOverPt = qSign / (radius * bMag);
    params[eFreeQOverP]  = qOverPt / fastHypot(1, lambda);
  }

  return params;
}

} // namespace
#endif

std::optional<Acts::FreeVector> eicrecon::acts_compat::estimateTrackParamsFromSpacePoints(
    std::span<const Acts::Vector3> points, const Acts::Vector3& field, double time,
    std::size_t refinementIterations, std::span<const double> weights, std::size_t referenceIndex) {
  if (points.size() < 3 || referenceIndex >= points.size() ||
      (!weights.empty() && weights.size() != points.size()) || !field.allFinite() ||
      !std::isfinite(time)) {
    return std::nullopt;
  }
  for (const auto& point : points) {
    if (!point.allFinite()) {
      return std::nullopt;
    }
  }
  for (const auto weight : weights) {
    if (weight < 0. || !std::isfinite(weight)) {
      return std::nullopt;
    }
  }
#if EICRECON_HAS_ACTS_MULTIPOINT
  const auto result = Acts::estimateTrackParamsFromSpacePoints(
      points, field, time, refinementIterations, weights, referenceIndex);
  if (!result.ok() || !result->allFinite()) {
    return std::nullopt;
  }
  return *result;
#else
  const auto result =
      fallbackEstimate(points, field, time, refinementIterations, weights, referenceIndex);
  if (!result || !result->allFinite()) {
    return std::nullopt;
  }
  return result;
#endif
}
