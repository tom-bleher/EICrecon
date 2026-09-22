// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <Acts/EventData/TrackStateType.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace eicrecon {

using B0SurfaceStationMap = std::unordered_map<std::uint64_t, unsigned int>;

/// Use physical station positions, independent of ACTS volume/layer numbering.
/// Positions are sensor-surface centres in the ion frame [mm]. Consecutive
/// positions separated by at most gap belong to the same physical station.
inline B0SurfaceStationMap
makeB0SurfaceStationMap(const std::vector<std::pair<std::uint64_t, double>>& surfaceZ, double gap) {
  if (!std::isfinite(gap) || gap <= 0.0) {
    throw std::invalid_argument("B0 station gap must be finite and positive");
  }

  // Repeated identical surfaces are harmless; conflicting identities are not.
  std::map<std::uint64_t, double> uniqueSurfaces;
  for (const auto& [id, z] : surfaceZ) {
    if (!std::isfinite(z)) {
      throw std::invalid_argument("B0 surface centre must be finite");
    }
    const auto [it, inserted] = uniqueSurfaces.emplace(id, z);
    if (!inserted && it->second != z) {
      throw std::invalid_argument("A B0 surface ID has conflicting positions");
    }
  }

  std::vector<std::pair<std::uint64_t, double>> surfaces(uniqueSurfaces.begin(),
                                                         uniqueSurfaces.end());
  std::sort(surfaces.begin(), surfaces.end(), [](const auto& a, const auto& b) {
    return a.second != b.second ? a.second < b.second : a.first < b.first;
  });

  B0SurfaceStationMap result;
  result.reserve(surfaces.size());
  unsigned int station = 0;
  double previousZ     = 0.0;
  for (const auto& [id, z] : surfaces) {
    if (!result.empty() && z - previousZ > gap) {
      ++station;
    }
    // Retain membership from clustering rather than reassigning by position.
    result.emplace(id, station);
    previousZ = z;
  }
  return result;
}

struct B0TrackStationCounts {
  std::size_t stations{};
  std::size_t unmappedMeasurements{};
};

/// Count only the fitted measurements, excluding holes, outliers and material states.
template <typename Track>
B0TrackStationCounts countB0TrackStations(const Track& track,
                                          const B0SurfaceStationMap& stationMap) {
  std::set<unsigned int> stations;
  std::size_t unmapped = 0;
  for (const auto& state : track.trackStatesReversed()) {
    const auto flags = state.typeFlags();
#if Acts_VERSION_MAJOR >= 45
    if (!flags.isMeasurement() || flags.isHole() || flags.isOutlier()) {
#else
    if (!flags.test(Acts::TrackStateFlag::MeasurementFlag) ||
        flags.test(Acts::TrackStateFlag::HoleFlag) ||
        flags.test(Acts::TrackStateFlag::OutlierFlag)) {
#endif
      continue;
    }
    const auto found = stationMap.find(state.referenceSurface().geometryId().value());
    if (found == stationMap.end()) {
      ++unmapped;
    } else {
      stations.insert(found->second);
    }
  }
  return {stations.size(), unmapped};
}

} // namespace eicrecon
