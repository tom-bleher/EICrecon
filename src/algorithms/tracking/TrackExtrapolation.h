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
#include <Acts/Utilities/Result.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>

namespace eicrecon {

/// Where extrapolateTrackInTwoStages left a track
enum class TwoStageExtrapolation {
  /// On the target surface, e.g. the perigee at the interaction point
  Target,
  /// On the reference plane, because the target could not be reached from there
  ReferencePlane,
};

/// Extrapolates a track in two stages. Stage 1 propagates from the first measurement against the
/// momentum, with the navigated propagator and so with the material of the tracking geometry, to
/// a plane across the beam just upstream of the measurements. Stage 2 propagates from that plane
/// to the target surface without navigation or material.
///
/// This is for a spectrometer far from the interaction point, such as the B0 tracker inside its
/// dipole. Acts::extrapolateTrackToReferenceSurface takes the direction from where the straight
/// tangent at a measurement meets the perigee line, which for a track bending towards the beam
/// line is downstream, so the track is lost or expressed after its measurements. Starting from
/// the first measurement and stepping backward to a plane avoids that. Without navigation, stage 2
/// does not pick up material between the plane and the target that does not belong on the track,
/// such as surfaces it crosses at a grazing angle. The field is still followed, since the solenoid
/// bends the track on its way to the interaction point.
///
/// A result on the target is accepted only if it lies upstream of the plane. Otherwise, as for a
/// decay product that passes the beam line nowhere near its origin, the track is left on the plane.
/// An error is returned, and the track is unchanged, if the plane cannot be reached.
template <typename track_proxy_t, typename navigated_t, typename navigated_options_t,
          typename free_t, typename free_options_t>
Acts::Result<TwoStageExtrapolation>
extrapolateTrackInTwoStages(track_proxy_t& track, const Acts::Surface& plane,
                            const Acts::Surface& target, const navigated_t& navigated,
                            navigated_options_t navigatedOptions, const free_t& free,
                            free_options_t freeOptions) {
  const auto& gctx = navigatedOptions.geoContext;
  auto first       = Acts::findFirstMeasurementState(track);
  if (!first.ok()) {
    return first.error();
  }
  const auto start = track.createParametersFromState(*first);
  // The plane must lie upstream of the first measurement, along the momentum
  const double planeOffset = (plane.center(gctx) - start.position(gctx)).dot(start.direction());
  if (!(planeOffset < 0.)) {
    return Acts::Result<TwoStageExtrapolation>::failure(
        Acts::TrackExtrapolationError::ReferenceSurfaceUnreachable);
  }

  navigatedOptions.direction = Acts::Direction::Backward();
#if Acts_VERSION_MAJOR >= 46
  auto toPlane = navigated.template propagate<navigated_options_t, Acts::SurfaceReached>(
      start, plane, navigatedOptions);
#else
  auto toPlane = navigated.template propagate<Acts::BoundTrackParameters, navigated_options_t,
                                              Acts::SurfaceReached>(start, plane, navigatedOptions);
#endif
  if (!toPlane.ok()) {
    return toPlane.error();
  }
  if (!toPlane->endParameters.has_value()) {
    return Acts::Result<TwoStageExtrapolation>::failure(
        Acts::TrackExtrapolationError::ReferenceSurfaceUnreachable);
  }
  const auto& onPlane = toPlane->endParameters.value();

  freeOptions.direction = Acts::Direction::Backward();
#if Acts_VERSION_MAJOR >= 46
  auto toTarget = free.template propagate<free_options_t, Acts::ForcedSurfaceReached>(
      onPlane, target, freeOptions);
#else
  auto toTarget = free.template propagate<Acts::BoundTrackParameters, free_options_t,
                                          Acts::ForcedSurfaceReached>(onPlane, target, freeOptions);
#endif
  const bool reached =
      toTarget.ok() && toTarget->endParameters.has_value() &&
      (toTarget->endParameters->position(gctx) - onPlane.position(gctx)).dot(onPlane.direction()) <
          0.;
  const auto& end = reached ? toTarget->endParameters.value() : onPlane;
  track.setReferenceSurface(end.referenceSurface().getSharedPtr());
  track.parameters() = end.parameters();
  track.covariance() = end.covariance().value();
  return reached ? TwoStageExtrapolation::Target : TwoStageExtrapolation::ReferencePlane;
}

} // namespace eicrecon
