// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0SeedTransport.h"

#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <cmath>
#include <tuple>

namespace eicrecon {

std::optional<Acts::BoundTrackParameters>
transportB0Seed(const Acts::BoundTrackParameters& start, double distance,
                const std::shared_ptr<const Acts::MagneticFieldProvider>& field,
                const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx) {
  if (!std::isfinite(distance) || distance < 0. || !start.parameters().allFinite() ||
      !start.covariance() || !start.covariance()->allFinite()) {
    return std::nullopt;
  }
  if (distance == 0.) {
    return start;
  }

  // A path-length target avoids the poorly conditioned closest-approach search
  // for a nearly axial track on a preselected upstream perigee line.
  using Stepper    = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper>;
  Propagator propagator{Stepper{field}};
  Propagator::Options<> options(gctx, mctx);
  options.direction = Acts::Direction::Backward();
  options.pathLimit = distance;
  auto result       = propagator.propagate(start, options);
  if (!result.ok() || !result->endParameters ||
      std::abs(result->pathLength + distance) > options.surfaceTolerance) {
    return std::nullopt;
  }

  const auto& endpoint = *result->endParameters;
  const auto perigee   = Acts::Surface::makeShared<Acts::PerigeeSurface>(endpoint.position(gctx));
  Stepper stepper{field};
  auto state = stepper.makeState(options.stepping);
  stepper.initialize(state, endpoint);
  // The curvilinear-to-perigee covariance conversion needs path derivatives,
  // including magnetic curvature, even though no further step is taken.
  if (!stepper.prepareCurvilinearState(state)) {
    return std::nullopt;
  }
  auto bound = stepper.boundState(state, *perigee);
  if (!bound.ok()) {
    return std::nullopt;
  }
  auto parameters = std::get<0>(*bound);
  if (!parameters.parameters().allFinite() || !parameters.covariance() ||
      !parameters.covariance()->allFinite()) {
    return std::nullopt;
  }
  return parameters;
}

} // namespace eicrecon
