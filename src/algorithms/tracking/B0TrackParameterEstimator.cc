// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#include "B0TelescopeActs.h"

#if EICRECON_HAS_B0_TELESCOPE
#include "B0TrackParameterEstimator.h"
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/GenericBoundTrackParameters.hpp>
#endif
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/PropagatorOptions.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>

namespace eicrecon::b0telescope {

std::optional<PerigeeEstimate> estimateTrackParameters(
    const std::array<Acts::Vector3, 3>& positions, double firstHitTime,
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    std::shared_ptr<const Acts::MagneticFieldProvider> field,
    const B0TelescopeSeedingConfig& cfg) {
  if (!field || !std::isfinite(firstHitTime) ||
      !std::all_of(positions.begin(), positions.end(), [](const auto& p) { return p.allFinite(); })) {
    return std::nullopt;
  }
  auto cache = field->makeCache(mctx);
  const auto sampled = field->getField(positions[1], cache);
  if (!sampled.ok() || !sampled->allFinite() || sampled->norm() < 1.e-9 * Acts::UnitConstants::T) {
    return std::nullopt;
  }
  // Avoid degeneracies in the conformal estimator before it normalizes vectors.
  const auto chord = (positions[2] - positions[0]).eval();
  if (chord.norm() < Acts::UnitConstants::um ||
      (positions[1] - positions[0]).norm() < Acts::UnitConstants::um ||
      (positions[2] - positions[1]).norm() < Acts::UnitConstants::um ||
      chord.cross(*sampled).norm() < 1.e-12 * chord.norm() * sampled->norm()) {
    return std::nullopt;
  }
  auto transform = Acts::Transform3::Identity();
  transform.translation() = positions[0];
  auto firstPlane = Acts::Surface::makeShared<Acts::PlaneSurface>(transform);
  const auto local = Acts::estimateTrackParamsFromSeed(
      gctx, *firstPlane, positions[0], firstHitTime, positions[1], positions[2], *sampled);
  if (!local.ok() || !local->allFinite() ||
      std::abs((*local)[Acts::eBoundQOverP]) < 1.e-12 / Acts::UnitConstants::GeV) {
    return std::nullopt; // Unresolved curvature: no invented momentum or charge sign.
  }

  constexpr double mm = Acts::UnitConstants::mm / edm4eic::unit::mm;
  constexpr double ns = Acts::UnitConstants::ns / edm4eic::unit::ns;
  constexpr double inverseGeV = edm4eic::unit::GeV / Acts::UnitConstants::GeV;
  const std::array<double, 6> sigmas{
      cfg.positionSigma * mm, cfg.positionSigma * mm,
      cfg.angleSigma * Acts::UnitConstants::rad, cfg.angleSigma * Acts::UnitConstants::rad,
      std::max(std::abs((*local)[Acts::eBoundQOverP]) * cfg.relativeMomentumSigma,
               cfg.qOverPSigmaFloor * inverseGeV), cfg.timeSigma * ns};
  auto covariance = Acts::BoundMatrix::Zero().eval();
  for (std::size_t i = 0; i < sigmas.size(); ++i) {
    covariance(i, i) = sigmas[i] * sigmas[i];
  }
  // Match the existing CKF's pion mass hypothesis. There is no PID decision here.
  const Acts::BoundTrackParameters initial(firstPlane, *local, covariance,
                                           Acts::ParticleHypothesis::pion());
  const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  // Field-only seed transport. Material effects belong to the downstream CKF.
  // In particular, never relabel the local hit tangent as a direction at the IP.
  Acts::Propagator<Acts::EigenStepper<>> propagator{Acts::EigenStepper<>{field}};
  decltype(propagator)::Options<> options(gctx, mctx);
  options.direction = Acts::Direction::Backward();
  options.pathLimit = cfg.propagationPathLimit * mm;
  const auto propagated = propagator.propagate(initial, *perigee, options);
  if (!propagated.ok() || !propagated->endParameters ||
      !propagated->endParameters->covariance()) {
    return std::nullopt;
  }
  const auto& end = *propagated->endParameters;
  PerigeeEstimate result{end.parameters(), *end.covariance()};
  if (!result.parameters.allFinite() || !result.covariance.allFinite() ||
      Eigen::LLT<Acts::BoundMatrix>(result.covariance).info() != Eigen::Success) {
    return std::nullopt;
  }
  return result;
}

} // namespace eicrecon::b0telescope
#endif
