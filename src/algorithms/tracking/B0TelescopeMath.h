// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace eicrecon::b0telescope {

// Numerical positions in mm. This small geometry kernel has no ACTS/JANA dependency.
struct Point {
  double x{}, y{}, z{};
};

inline bool finite(const Point& p) {
  return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

// The realistic readout numbers front/back layers 1..8; the ideal one uses 1..4.
// This is an explicit readout contract, not an inference from the hits in an event.
inline std::optional<unsigned> stationFromLayer(int layer, unsigned layersPerStation) {
  if ((layersPerStation != 1 && layersPerStation != 2) || layer < 1 ||
      layer > static_cast<int>(4 * layersPerStation)) {
    return std::nullopt;
  }
  return static_cast<unsigned>(layer - 1) / layersPerStation;
}

// Inverse of the detector's rotation about laboratory y; no IP constraint.
inline Point toTelescopeFrame(const Point& p, double rotation) {
  const double c = std::cos(rotation), s = std::sin(rotation);
  return {c * p.x - s * p.z, p.y, s * p.x + c * p.z};
}

inline constexpr std::array<std::array<unsigned, 3>, 4> stationTriples{{
    {{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}},
}};

inline bool compatibleDoublet(const Point& first, const Point& last,
                              double minDeltaZ, double maxSlope) {
  if (!finite(first) || !finite(last)) {
    return false;
  }
  const double dz = last.z - first.z;
  return dz >= minDeltaZ && dz > 0.0 &&
         std::abs(last.x - first.x) <= maxSlope * dz &&
         std::abs(last.y - first.y) <= maxSlope * dz;
}

// Rectangular sagitta window. X is the bending plane, so its window is wider.
// This is a compatibility score, NOT a fitted chi-square or an IP distance.
inline std::optional<double> tripletScore(const Point& a, const Point& b, const Point& c,
                                         double maxResidualX, double maxResidualY) {
  if (!finite(a) || !finite(b) || !finite(c) || !(a.z < b.z && b.z < c.z) ||
      !(maxResidualX > 0.0 && maxResidualY > 0.0)) {
    return std::nullopt;
  }
  const double f = (b.z - a.z) / (c.z - a.z);
  const double dx = (b.x - std::lerp(a.x, c.x, f)) / maxResidualX;
  const double dy = (b.y - std::lerp(a.y, c.y, f)) / maxResidualY;
  if (!std::isfinite(dx) || !std::isfinite(dy) || std::abs(dx) > 1.0 || std::abs(dy) > 1.0) {
    return std::nullopt;
  }
  return dx * dx + dy * dy;
}

// Bound the *whole event* before seeding. Do not keep an input-order-biased prefix
// when a busy event exceeds this starter implementation's combinatorial budget.
inline bool withinCombinationBudget(const std::array<std::size_t, 4>& counts,
                                    std::size_t budget) {
  std::size_t total = 0;
  for (const auto& triple : stationTriples) {
    if (counts[triple[0]] == 0 || counts[triple[1]] == 0 || counts[triple[2]] == 0) {
      continue;
    }
    std::size_t product = 1;
    for (const auto station : triple) {
      if (counts[station] > (budget - total) / product) {
        return false;
      }
      product *= counts[station];
    }
    total += product;
  }
  return true;
}

} // namespace eicrecon::b0telescope
