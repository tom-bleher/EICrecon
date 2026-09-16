// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace eicrecon::b0 {

/// Run-scoped physical-station identity, independent of ACTS layer numbering.
/// Inputs must be sensor-surface centres transformed into the ion frame [mm].
class StationMap {
public:
  struct Surface {
    std::uint64_t id{};
    double z{};
  };
  struct Station {
    double zMin{};
    double zMax{};
    double zMean{};
    std::vector<std::uint64_t> surfaces;
  };

  explicit StationMap(std::vector<Surface> surfaces, double gap = 50.0) : m_gap(gap) {
    if (!std::isfinite(gap) || gap <= 0.0) {
      throw std::invalid_argument("B0 station gap must be finite and positive");
    }
    std::map<std::uint64_t, double> unique;
    for (const auto& surface : surfaces) {
      if (!std::isfinite(surface.z)) {
        throw std::invalid_argument("B0 surface centre must be finite");
      }
      const auto [it, inserted] = unique.emplace(surface.id, surface.z);
      if (!inserted && it->second != surface.z) {
        throw std::invalid_argument("A B0 surface ID has conflicting positions");
      }
    }
    surfaces.clear();
    for (const auto& [id, z] : unique) {
      surfaces.push_back({id, z});
    }
    std::sort(surfaces.begin(), surfaces.end(), [](const Surface& a, const Surface& b) {
      return a.z != b.z ? a.z < b.z : a.id < b.id;
    });
    for (const auto& surface : surfaces) {
      if (m_stations.empty() || surface.z - m_stations.back().zMax > m_gap) {
        m_stations.push_back({surface.z, surface.z, surface.z, {surface.id}});
      } else {
        auto& station = m_stations.back();
        station.zMax = surface.z;
        station.zMean += (surface.z - station.zMean) / (station.surfaces.size() + 1.0);
        station.surfaces.push_back(surface.id);
      }
      m_surfaceStation.emplace(surface.id, m_stations.size() - 1);
    }
  }

  [[nodiscard]] std::optional<std::size_t> stationForSurface(std::uint64_t id) const {
    const auto found = m_surfaceStation.find(id);
    return found == m_surfaceStation.end() ? std::nullopt
                                          : std::optional<std::size_t>{found->second};
  }

  /// Diagnostic position fallback only. Production consumers should use IDs.
  /// A tie is deliberately unresolved rather than depending on container order.
  [[nodiscard]] std::optional<std::size_t> stationForZ(double z) const {
    if (!std::isfinite(z)) {
      return std::nullopt;
    }
    std::optional<std::size_t> result;
    double distance = std::numeric_limits<double>::infinity();
    bool tied = false;
    for (std::size_t i = 0; i < m_stations.size(); ++i) {
      const auto& station = m_stations[i];
      if (z < station.zMin - m_gap || z > station.zMax + m_gap) {
        continue;
      }
      const double candidate = std::abs(z - station.zMean);
      if (candidate < distance) {
        distance = candidate;
        result = i;
        tied = false;
      } else if (candidate == distance) {
        tied = true;
      }
    }
    return tied ? std::nullopt : result;
  }

  [[nodiscard]] const std::vector<Station>& stations() const { return m_stations; }
  [[nodiscard]] double gap() const { return m_gap; }
  [[nodiscard]] std::size_t surfaceCount() const { return m_surfaceStation.size(); }

private:
  double m_gap;
  std::vector<Station> m_stations;
  std::map<std::uint64_t, std::size_t> m_surfaceStation;
};

inline double ionFrameZ(double x, double z, double crossingAngle) {
  return x * std::sin(crossingAngle) + z * std::cos(crossingAngle);
}

} // namespace eicrecon::b0
