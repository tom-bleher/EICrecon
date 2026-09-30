// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
//
// Calibrated two-pad position curve: the distance of a deposit from the leading pad centre,
// towards the neighbour, as a function of the neighbour's fraction f = A_n / (A_lead + A_n).
// Derived from simulation truth or test-beam tracks for a given sensor and readout.
//
// Text format ('#' starts a comment line): rows "f offset_mm" with increasing f in [0, 0.5].

#pragma once

#include <algorithm>
#include <istream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace eicrecon {

class SharedOffsetCurve {
public:
  explicit SharedOffsetCurve(std::istream& in) {
    for (std::string line; std::getline(in, line);) {
      if (line.find_first_not_of(" \t") == std::string::npos || line.front() == '#') {
        continue;
      }
      double f = 0, offset = 0;
      if (!(std::istringstream(line) >> f >> offset) || (!m_f.empty() && f <= m_f.back())) {
        throw std::runtime_error("SharedOffsetCurve: bad or non-increasing row: " + line);
      }
      m_f.push_back(f);
      m_offset.push_back(offset);
    }
    if (m_f.size() < 2) {
      throw std::runtime_error("SharedOffsetCurve: need at least two rows");
    }
  }

  /// Offset (mm) for fraction f, interpolated linearly and clamped to the tabulated range
  double offset(double f) const {
    if (f <= m_f.front()) {
      return m_offset.front();
    }
    if (f >= m_f.back()) {
      return m_offset.back();
    }
    const auto i = static_cast<std::size_t>(
        std::distance(m_f.begin(), std::upper_bound(m_f.begin(), m_f.end(), f)));
    const double t = (f - m_f[i - 1]) / (m_f[i] - m_f[i - 1]);
    return m_offset[i - 1] + t * (m_offset[i] - m_offset[i - 1]);
  }

private:
  std::vector<double> m_f;
  std::vector<double> m_offset;
};

} // namespace eicrecon
