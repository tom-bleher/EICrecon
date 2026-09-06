// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <cstddef>
#include <optional>
#include <span>

namespace eicrecon::acts_compat {
// Native ACTS multipoint estimator when available, otherwise the MPL-2.0
// ACTS 47.7 implementation. Positions/field use ACTS units; points must be in
// track order. The homogeneous-field estimate is at points[referenceIndex].
// Invalid inputs or nonfinite/degenerate estimates return nullopt.
std::optional<Acts::FreeVector> estimateTrackParamsFromSpacePoints(
    std::span<const Acts::Vector3> points, const Acts::Vector3& field, double time = 0.,
    std::size_t refinementIterations = 0, std::span<const double> weights = {},
    std::size_t referenceIndex = 0);
} // namespace eicrecon::acts_compat
