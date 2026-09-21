// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/UnitVectors.hpp>
#include <Eigen/Core>
#include <cmath>

namespace eicrecon {

/// Convert a charged track's bound covariance to [x,y,z,px,py,pz], in mm and GeV.
/// Requires finite, nonzero q/p. The surface Jacobian includes the dependence
/// of perigee position on direction when the impact parameter is nonzero.
inline Eigen::Matrix<double, 6, 6>
boundToCartesianCovariance(const Acts::Surface& surface, const Acts::GeometryContext& gctx,
                           const Acts::BoundVector& parameters,
                           const Eigen::Matrix<double, 6, 6>& covariance,
                           double absoluteCharge = Acts::UnitConstants::e) {
  const auto direction =
      Acts::makeDirectionFromPhiTheta(parameters[Acts::eBoundPhi], parameters[Acts::eBoundTheta]);
  const auto position    = surface.localToGlobal(gctx, parameters.head<2>(), direction);
  const auto boundToFree = surface.boundToFreeJacobian(gctx, position, direction);
  const double qOverP    = parameters[Acts::eBoundQOverP];
  const double momentum  = std::abs(absoluteCharge / qOverP);

  Eigen::Matrix<double, 6, 6> jacobian;
  jacobian.topRows<3>()    = boundToFree.block<3, 6>(Acts::eFreePos0, 0) / Acts::UnitConstants::mm;
  jacobian.bottomRows<3>() = momentum * boundToFree.block<3, 6>(Acts::eFreeDir0, 0);
  jacobian.block<3, 1>(3, Acts::eBoundQOverP) -= momentum / qOverP * direction;
  jacobian.bottomRows<3>() /= Acts::UnitConstants::GeV;
  return jacobian * covariance * jacobian.transpose();
}

} // namespace eicrecon
