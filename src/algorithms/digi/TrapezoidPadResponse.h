// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
//
// Pad response of AC-coupled LGAD pixels along one axis: a deposit on the metal electrode is
// collected by that pad alone, and a deposit in the inter-pad gap is shared linearly between
// the two pads whose metal edges bound the gap (Tornago et al., NIM A 1003 (2021) 165319;
// Dutta et al., NIM A (2025) 170224). The 2D response is the product of the x and y responses.

#pragma once

#include <algorithm>
#include <cmath>

namespace eicrecon {

/// Fraction of a deposit at `offset` from a pad centre collected by that pad.
/// Equals 1 on the metal (|offset| <= electrode/2), falls linearly across the gap and reaches 0
/// at the metal edge of the neighbouring pad. An electrode as wide as the pitch means no sharing.
inline double trapezoidPadShare(double offset, double pitch, double electrode) {
  if (electrode >= pitch) {
    return std::abs(offset) < 0.5 * pitch ? 1.0 : 0.0;
  }
  return std::clamp((pitch - 0.5 * electrode - std::abs(offset)) / (pitch - electrode), 0.0, 1.0);
}

/// Inverse of trapezoidPadShare for a deposit shared by two adjacent pads: distance of the
/// deposit from the leading pad centre, towards the neighbour, given the neighbour's fraction
/// f = A_neighbour / (A_lead + A_neighbour) of the pair amplitude.
inline double trapezoidSharedOffset(double f, double pitch, double electrode) {
  const double metal = std::min(electrode, pitch);
  return 0.5 * metal + std::clamp(f, 0.0, 0.5) * (pitch - metal);
}

} // namespace eicrecon
