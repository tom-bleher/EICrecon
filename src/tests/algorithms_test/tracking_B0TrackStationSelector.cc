// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <catch2/catch_test_macros.hpp>
#include <Acts/EventData/TrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
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

TEST_CASE("B0 station mapping preserves four-station layouts", "[tracking][b0stations]") {
  using eicrecon::makeB0SurfaceStationMap;
  const auto disks = makeB0SurfaceStationMap({{9, 5902.}, {8, 6172.}, {7, 6442.}, {6, 6712.}}, 50.);
  REQUIRE(disks.size() == 4);
  CHECK(disks.at(9) == 0);
  CHECK(disks.at(8) == 1);
  CHECK(disks.at(7) == 2);
  CHECK(disks.at(6) == 3);

  std::vector<std::pair<std::uint64_t, double>> faces{{91, 6000.}, {23, 6007.}, {5, 6276.},
                                                      {1, 6283.},  {17, 6525.}, {18, 6532.},
                                                      {48, 6750.}, {49, 6757.}};
  const auto mapped = makeB0SurfaceStationMap(faces, 50.);
  REQUIRE(mapped.size() == 8);
  for (std::size_t i = 0; i < faces.size(); ++i) {
    CHECK(mapped.at(faces[i].first) == i / 2);
  }
  std::ranges::reverse(faces);
  CHECK(mapped == makeB0SurfaceStationMap(faces, 50.));
}

TEST_CASE("B0 station mapping validates surface identity", "[tracking][b0stations]") {
  using eicrecon::makeB0SurfaceStationMap;
  const auto expected = makeB0SurfaceStationMap({{91, 6000.}, {23, 6007.}}, 50.);
  CHECK(expected == makeB0SurfaceStationMap({{23, 6007.}, {91, 6000.}, {23, 6007.}}, 50.));
  const auto coincident = makeB0SurfaceStationMap({{91, 6000.}, {23, 6000.}}, 50.);
  REQUIRE(coincident.size() == 2);
  CHECK(coincident.at(91) == coincident.at(23));
  CHECK_THROWS_AS(makeB0SurfaceStationMap({{91, 6000.}, {91, 6007.}}, 50.), std::invalid_argument);
  CHECK_THROWS_AS(makeB0SurfaceStationMap({{91, 6007.}, {91, 6000.}}, 50.), std::invalid_argument);
}

TEST_CASE("B0 station mapping rejects invalid gaps and positions", "[tracking][b0stations]") {
  using eicrecon::makeB0SurfaceStationMap;
  const double infinity = std::numeric_limits<double>::infinity();
  const double nan      = std::numeric_limits<double>::quiet_NaN();
  for (const double gap : std::array{0.0, -1.0, infinity, -infinity, nan}) {
    CHECK_THROWS_AS(makeB0SurfaceStationMap({{1, 6000.}}, gap), std::invalid_argument);
    CHECK_THROWS_AS(makeB0SurfaceStationMap({}, gap), std::invalid_argument);
  }
  for (const double z : std::array{infinity, -infinity, nan}) {
    CHECK_THROWS_AS(makeB0SurfaceStationMap({{1, z}}, 50.), std::invalid_argument);
  }
}

TEST_CASE("B0 station mapping retains single-linkage gap semantics", "[tracking][b0stations]") {
  using eicrecon::makeB0SurfaceStationMap;
  CHECK(makeB0SurfaceStationMap({}, 50.).empty());
  const auto one = makeB0SurfaceStationMap({{91, 6000.}}, 50.);
  REQUIRE(one.size() == 1);
  CHECK(one.at(91) == 0);
  CHECK(one.find(999) == one.end());

  // Equality stays in the same station; compare consecutive positions, not
  // distance from the first surface or the station mean.
  const auto mapped = makeB0SurfaceStationMap({{4, 100.01}, {3, 50.}, {2, 0.}, {1, -50.}}, 50.);
  REQUIRE(mapped.size() == 4);
  CHECK(mapped.at(1) == 0);
  CHECK(mapped.at(2) == 0);
  CHECK(mapped.at(3) == 0);
  CHECK(mapped.at(4) == 1);
}
