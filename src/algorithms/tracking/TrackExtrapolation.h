// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Direction.hpp>
#include <Acts/Definitions/Units.hpp>
#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/TrackParameters.hpp>
#endif
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>
#include <cmath>

namespace eicrecon {

/// Whether the reference point of a track lies before its first measurement along the momentum,
/// as for a track expressed where it was produced
template <typename track_proxy_t>
bool isReferenceBeforeFirstMeasurement(const track_proxy_t& track,
                                       const Acts::GeometryContext& gctx) {
  auto first = Acts::findFirstMeasurementState(track);
  if (!first.ok()) {
    return false;
  }
  const auto measured = track.createParametersFromState(*first);
  return (track.createParametersAtReference().position(gctx) - measured.position(gctx))
             .dot(measured.direction()) < 0;
}

/// Extrapolates a track from its first measurement against the momentum to a surface, e.g. a
/// perigee, as a track starts before its first measurement. A positive `approach`
/// first steps that far back without a target: inside a dipole the tangent of a track bending
/// towards the beam line meets it downstream, and a propagation that targets the beam line from
/// there ends downstream as well.
template <typename track_proxy_t, typename propagator_t, typename propagator_options_t>
Acts::Result<void> extrapolateTrackBackward(track_proxy_t& track, const Acts::Surface& surface,
                                            const propagator_t& propagator,
                                            propagator_options_t options, double approach = 0.) {
  auto first = Acts::findFirstMeasurementState(track);
  if (!first.ok()) {
    return first.error();
  }
  options.direction = Acts::Direction::Backward();
  auto start        = track.createParametersFromState(*first);
  if (approach > 0.) {
    auto stepOptions      = options;
    stepOptions.pathLimit = approach;
    auto step             = propagator.propagate(start, stepOptions);
    if (!step.ok()) {
      return step.error();
    }
    // A step that stops short has left the tracking world, which then does not reach the surface
    if (std::abs(step->pathLength) < approach - 1. * Acts::UnitConstants::mm) {
      return Acts::Result<void>::failure(
          Acts::TrackExtrapolationError::ReferenceSurfaceUnreachable);
    }
    start = step->endParameters.value();
  }
#if Acts_VERSION_MAJOR >= 46
  auto propagation =
      propagator.template propagate<propagator_options_t, Acts::ForcedSurfaceReached>(
          start, surface, options);
#else
  auto propagation =
      propagator.template propagate<Acts::BoundTrackParameters, propagator_options_t,
                                    Acts::ForcedSurfaceReached>(start, surface, options);
#endif
  if (!propagation.ok()) {
    return propagation.error();
  }
  track.setReferenceSurface(surface.getSharedPtr());
  track.parameters() = propagation->endParameters.value().parameters();
  track.covariance() = propagation->endParameters.value().covariance().value();
  return Acts::Result<void>::success();
}

/// Extrapolates a track to the perigee surface with Acts::extrapolateTrackToReferenceSurface and
/// accepts the result only if it lies before the first measurement. The helper propagates towards
/// where the straight tangent at the first or last state meets the perigee line; for a track bending
/// towards the beam line in a dipole that is downstream, where the propagation either leaves the
/// world or ends after the measurements. The extrapolation is then retried backward, first half-way
/// towards the perigee. A displaced track that passes the beam line nowhere near its origin is kept
/// on the fallback surface, e.g. the perigee of its seed, unless that is the perigee itself.
///
/// If every attempt fails, an error is returned and the track parameters are not usable.
template <typename track_proxy_t, typename propagator_t, typename propagator_options_t>
Acts::Result<void>
extrapolateTrackToPerigee(track_proxy_t& track, const Acts::Surface& perigee,
                          const Acts::Surface& fallback, const propagator_t& propagator,
                          const propagator_options_t& options, const Acts::Logger& logger) {
  const auto& gctx = options.geoContext;
  // A reference point after the first measurement is as wrong as a failed extrapolation
  auto check = [&](Acts::Result<void> result, const char* attempt) {
    if (result.ok() && !isReferenceBeforeFirstMeasurement(track, gctx)) {
      result =
          Acts::Result<void>::failure(Acts::TrackExtrapolationError::ReferenceSurfaceUnreachable);
    }
    if (!result.ok()) {
      ACTS_DEBUG(attempt << " extrapolation failed with error " << result.error().message());
    }
    return result;
  };

  auto result = check(Acts::extrapolateTrackToReferenceSurface(
                          track, perigee, propagator, options,
                          Acts::TrackExtrapolationStrategy::firstOrLast, logger),
                      "Nearest-state");
  if (result.ok()) {
    return result;
  }

  auto first = Acts::findFirstMeasurementState(track);
  if (!first.ok()) {
    return first.error();
  }
  const double approach =
      0.5 * (track.createParametersFromState(*first).position(gctx) - perigee.center(gctx)).norm();
  result =
      check(extrapolateTrackBackward(track, perigee, propagator, options, approach), "Backward");
  if (result.ok() || fallback.center(gctx) == perigee.center(gctx)) {
    return result;
  }
  return check(extrapolateTrackBackward(track, fallback, propagator, options), "Fallback");
}

} // namespace eicrecon
