// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#pragma once

#include "B0TelescopeActs.h"
#include "B0TelescopeSeedingConfig.h"

#if EICRECON_HAS_B0_TELESCOPE
#include <array>
#include <cstdint>

namespace eicrecon::b0telescope {

// Points are stored contiguously by physical station. Offsets delimit four ranges.
// All coordinates in the ACTS container are laboratory coordinates, in ACTS units.
void findTelescopeSeeds(const SpacePoints& points, const std::array<std::uint32_t, 5>& offsets,
                        double detectorRotation, const B0TelescopeSeedingConfig& config,
                        Seeds& output);

} // namespace eicrecon::b0telescope
#endif
