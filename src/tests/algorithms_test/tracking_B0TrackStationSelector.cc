// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <catch2/catch_test_macros.hpp>
#include <Acts/EventData/TrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <algorithm>
#include <memory>
#include <vector>

#include "algorithms/tracking/B0TrackStationSelector.h"

TEST_CASE("B0 final track selection counts physical stations, not sensors",
          "[tracking][b0stations]") {
  // Deliberately scrambled IDs and input order: front/back stacks and staggered
  // chips share a station, regardless of ACTS sensitive/layer/volume numbering.
  std::vector<std::pair<std::uint64_t, double>> surfaces{{91, 6850.},   {74, 6096.}, {12, 6104.},
                                                         {33, 6097.65}, {58, 6376.}, {9, 6625.}};
  const auto stationMap = eicrecon::makeB0SurfaceStationMap(surfaces, 50.);
  std::ranges::reverse(surfaces);
  CHECK(stationMap == eicrecon::makeB0SurfaceStationMap(surfaces, 50.));
  REQUIRE(stationMap.size() == surfaces.size());
  CHECK(stationMap.at(74) == stationMap.at(12));
  CHECK(stationMap.at(74) == stationMap.at(33));
  CHECK(stationMap.at(74) != stationMap.at(58));

  ActsExamples::TrackContainer tracks{std::make_shared<Acts::VectorTrackContainer>(),
                                      std::make_shared<Acts::VectorMultiTrajectory>()};
  auto track     = tracks.makeTrack();
  const auto add = [&](std::uint64_t id, bool outlier = false, bool hole = false) {
    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Transform3::Identity());
    surface->assignGeometryId(Acts::GeometryIdentifier{id});
    auto state = track.appendTrackState(Acts::TrackStatePropMask::None);
    state.setReferenceSurface(surface);
#if Acts_VERSION_MAJOR >= 45
    if (hole) {
      state.typeFlags().setIsHole();
    } else if (outlier) {
      state.typeFlags().setIsOutlier();
    } else {
      state.typeFlags().setIsMeasurement();
    }
#else
    state.typeFlags().set(hole      ? Acts::TrackStateFlag::HoleFlag
                          : outlier ? Acts::TrackStateFlag::OutlierFlag
                                    : Acts::TrackStateFlag::MeasurementFlag);
#endif
  };

  add(74);
  add(12);
  add(33);
  CHECK(eicrecon::countB0TrackStations(track, stationMap).stations == 1);
  add(58);
  add(9, true);
  add(91, false, true);
  // Four measurements plus an outlier and hole still span only two stations.
  CHECK(eicrecon::countB0TrackStations(track, stationMap).stations == 2);
  add(999); // An unknown surface cannot manufacture a third B0 station.
  const auto missing = eicrecon::countB0TrackStations(track, stationMap);
  CHECK(missing.stations == 2);
  CHECK(missing.unmappedMeasurements == 1);
  add(9);
  CHECK(eicrecon::countB0TrackStations(track, stationMap).stations == 3);
  add(91);
  CHECK(eicrecon::countB0TrackStations(track, stationMap).stations == 4);
}
