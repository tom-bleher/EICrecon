// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Direction.hpp>
#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/TrackParameters.hpp>
#endif
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>

namespace eicrecon {

/// Extrapolates a track from its first measurement against the momentum to a surface,
/// as a track starts before its first measurement
template <typename track_proxy_t, typename propagator_t, typename propagator_options_t>
Acts::Result<void> extrapolateTrackBackward(track_proxy_t& track, const Acts::Surface& surface,
                                            const propagator_t& propagator,
                                            propagator_options_t options) {
  auto first = Acts::findFirstMeasurementState(track);
  if (!first.ok()) {
    return first.error();
  }
  options.direction = Acts::Direction::Backward();
  auto propagation = propagator.template propagate<Acts::BoundTrackParameters, propagator_options_t,
                                                   Acts::ForcedSurfaceReached>(
      track.createParametersFromState(*first), surface, options);
  if (!propagation.ok()) {
    return propagation.error();
  }
  track.setReferenceSurface(surface.getSharedPtr());
  track.parameters() = propagation->endParameters.value().parameters();
  track.covariance() = propagation->endParameters.value().covariance().value();
  return Acts::Result<void>::success();
}

/// Extrapolates a track to the perigee surface with Acts::extrapolateTrackToReferenceSurface.
/// The helper propagates towards where the straight tangent at the first or last state meets the
/// perigee line. The tangent of a track bending towards the beam line in a dipole meets it
/// downstream, and the propagation leaves the world. The extrapolation is then retried backward;
/// a displaced track that passes the beam line nowhere near its origin is kept on the fallback
/// surface, e.g. the perigee of its seed.
template <typename track_proxy_t, typename propagator_t, typename propagator_options_t>
Acts::Result<void>
extrapolateTrackToPerigee(track_proxy_t& track, const Acts::Surface& perigee,
                          const Acts::Surface& fallback, const propagator_t& propagator,
                          const propagator_options_t& options, const Acts::Logger& logger) {
  auto result = Acts::extrapolateTrackToReferenceSurface(
      track, perigee, propagator, options, Acts::TrackExtrapolationStrategy::firstOrLast, logger);
  if (result.ok()) {
    return result;
  }
  ACTS_DEBUG("Extrapolation failed with error " << result.error().message()
                                                << ", retrying backward");
  result = extrapolateTrackBackward(track, perigee, propagator, options);
  if (result.ok()) {
    return result;
  }
  return extrapolateTrackBackward(track, fallback, propagator, options);
}

} // namespace eicrecon
