// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <array>
#include <memory>
#include <optional>
#include "B0TelescopeSeedingConfig.h"

namespace eicrecon::b0telescope {

struct PerigeeEstimate {
  Acts::BoundVector parameters;
  Acts::BoundMatrix covariance;
};

// Input positions/time and returned parameters/covariance use ACTS units.
// No beam-spot constraint: the origin perigee is a reference surface, not a vertex.
std::optional<PerigeeEstimate> estimateTrackParameters(
    const std::array<Acts::Vector3, 3>& positions, double firstHitTime,
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    std::shared_ptr<const Acts::MagneticFieldProvider> field,
    const B0TelescopeSeedingConfig& config);

} // namespace eicrecon::b0telescope
