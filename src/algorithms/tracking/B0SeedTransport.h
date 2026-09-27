// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/TrackParameters.hpp>
#endif
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <memory>
#include <optional>

namespace eicrecon {

/// Transport a seed upstream by a fixed path length in the magnetic field and
/// express it, including covariance, on a perigee centered at the endpoint.
/// This local reference change intentionally has no material interaction.
std::optional<Acts::BoundTrackParameters>
transportB0Seed(const Acts::BoundTrackParameters& start, double distance,
                const std::shared_ptr<const Acts::MagneticFieldProvider>& field,
                const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx);

} // namespace eicrecon
