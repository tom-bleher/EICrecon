// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
//
// Tabulated pad response for pixel sensors with signal sharing: the fraction of a deposit that
// each pad of a (2 radius + 1)^2 neighbourhood collects, on a grid of deposit positions over
// the pad containing it. Tables can come from measured maps or detailed simulation.
//
// Text format (lengths in mm, '#' starts a comment line):
//   pitch_x pitch_y nx ny radius
//   then nx * ny rows, x index slowest, each with (2 radius + 1)^2 fractions ordered by pad
//   offset (di, dj), di slowest, for deposits at the centres of an nx by ny grid over the pad.

#pragma once

#include <algorithm>
#include <cmath>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace eicrecon {

class PadResponseTable {
public:
  explicit PadResponseTable(std::istream& in) {
    std::stringstream data;
    for (std::string line; std::getline(in, line);) {
      if (line.find_first_not_of(" \t") != std::string::npos && line.front() != '#') {
        data << line << '\n';
      }
    }
    if (!(data >> m_pitch_x >> m_pitch_y >> m_nx >> m_ny >> m_radius) || m_nx < 2 || m_ny < 2 ||
        m_radius < 0) {
      throw std::runtime_error("PadResponseTable: bad header");
    }
    m_values.resize(static_cast<std::size_t>(m_nx * m_ny * pads()));
    for (auto& v : m_values) {
      if (!(data >> v)) {
        throw std::runtime_error("PadResponseTable: too few values");
      }
    }
  }

  double pitchX() const { return m_pitch_x; }
  double pitchY() const { return m_pitch_y; }

  /// Fraction collected by a pad for a deposit at offset (dx, dy) from its centre (mm)
  double fraction(double dx, double dy) const {
    // The deposit sits in the pad k pitches away from this one, at offset (ux, uy) from its
    // centre; this pad is then at offset -k from the pad containing the deposit
    const auto kx = static_cast<int>(std::lround(dx / m_pitch_x));
    const auto ky = static_cast<int>(std::lround(dy / m_pitch_y));
    if (std::abs(kx) > m_radius || std::abs(ky) > m_radius) {
      return 0.0;
    }
    const double ux = dx - kx * m_pitch_x;
    const double uy = dy - ky * m_pitch_y;
    const int pad   = (m_radius - kx) * (2 * m_radius + 1) + (m_radius - ky);
    // Bilinear interpolation between the grid points at bin centres, clamped at the pad edges
    const double fx = std::clamp((ux / m_pitch_x + 0.5) * m_nx - 0.5, 0.0, m_nx - 1.0);
    const double fy = std::clamp((uy / m_pitch_y + 0.5) * m_ny - 0.5, 0.0, m_ny - 1.0);
    const int ix    = std::min(static_cast<int>(fx), m_nx - 2);
    const int iy    = std::min(static_cast<int>(fy), m_ny - 2);
    const double tx = fx - ix;
    const double ty = fy - iy;
    return (1 - tx) * (1 - ty) * at(ix, iy, pad) + tx * (1 - ty) * at(ix + 1, iy, pad) +
           (1 - tx) * ty * at(ix, iy + 1, pad) + tx * ty * at(ix + 1, iy + 1, pad);
  }

private:
  int pads() const { return (2 * m_radius + 1) * (2 * m_radius + 1); }
  double at(int ix, int iy, int pad) const {
    return m_values[static_cast<std::size_t>((ix * m_ny + iy) * pads() + pad)];
  }

  double m_pitch_x{0};
  double m_pitch_y{0};
  int m_nx{0};
  int m_ny{0};
  int m_radius{0};
  std::vector<double> m_values;
};

} // namespace eicrecon
