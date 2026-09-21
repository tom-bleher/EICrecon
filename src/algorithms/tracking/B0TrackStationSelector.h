// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <Acts/EventData/TrackStateType.hpp>
#include <cstddef>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "B0TrackerStubSeeder.h"

namespace eicrecon {

using B0SurfaceStationMap = std::unordered_map<std::uint64_t, unsigned int>;

/// Use physical station positions, independent of ACTS volume/layer numbering.
inline B0SurfaceStationMap
makeB0SurfaceStationMap(const std::vector<std::pair<std::uint64_t, double>>& surfaceZ, double gap) {
  std::vector<double> positions;
  for (const auto& [id, z] : surfaceZ) {
    positions.push_back(z);
  }
  const auto stations = b0stub::clusterStations(positions, gap);
  B0SurfaceStationMap result;
  for (const auto& [id, z] : surfaceZ) {
    const int station = b0stub::assignStation(z, stations, gap);
    if (station >= 0) {
      result.emplace(id, static_cast<unsigned int>(station));
    }
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
