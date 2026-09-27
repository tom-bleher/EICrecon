// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <Acts/EventData/TrackStateType.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace eicrecon {

using SurfaceStationMap = std::unordered_map<std::uint64_t, std::size_t>;

/// Group the complete sensitive geometry of a telescope by gaps in global z.
/// Front/back faces and staggered modules share a station; ACTS layer IDs need
/// not be consecutive. Build this once, independently of which hits an event has.
inline SurfaceStationMap
makeSurfaceStationMap(std::vector<std::pair<std::uint64_t, double>> surfaceZ, double gap) {
  if (!std::isfinite(gap) || gap <= 0.) {
    throw std::invalid_argument("StationZGap must be finite and positive");
  }
  for (const auto& [id, z] : surfaceZ) {
    if (!std::isfinite(z)) {
      throw std::invalid_argument("Station surface position must be finite");
    }
  }
  std::ranges::sort(surfaceZ, [](const auto& a, const auto& b) {
    return a.second < b.second || (a.second == b.second && a.first < b.first);
  });
  SurfaceStationMap result;
  std::size_t station = 0;
  double lastZ        = 0.;
  for (const auto& [id, z] : surfaceZ) {
    if (!result.empty() && z - lastZ > gap) {
      ++station;
    }
    if (!result.emplace(id, station).second) {
      throw std::invalid_argument("Duplicate station surface identifier");
    }
    lastZ = z;
  }
  return result;
}

/// Require accepted measurements in distinct physical stations. Unknown
/// measurement surfaces fail the selection rather than inventing coverage.
template <typename Track>
bool hasTrackStationCoverage(const Track& track, const SurfaceStationMap& stationMap,
                             std::size_t minimum) {
  if (minimum == 0) {
    return true;
  }
  std::set<std::size_t> stations;
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
      return false;
    }
    stations.insert(found->second);
  }
  return stations.size() >= minimum;
}

} // namespace eicrecon
