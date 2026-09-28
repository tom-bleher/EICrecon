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
#include <limits>
#include <stdexcept>
#include <vector>

#include "algorithms/tracking/TrackStationSelector.h"

TEST_CASE("B0 final track selection counts physical stations, not sensors",
          "[tracking][b0stations]") {
  // Deliberately scrambled IDs and input order: front/back stacks and staggered
  // chips share a station, regardless of ACTS sensitive/layer/volume numbering.
  std::vector<std::pair<std::uint64_t, double>> surfaces{{91, 6850.},   {74, 6096.}, {12, 6104.},
                                                         {33, 6097.65}, {58, 6376.}, {9, 6625.}};
  const auto stationMap = eicrecon::makeSurfaceStationMap(surfaces, 50.);
  std::ranges::reverse(surfaces);
  CHECK(stationMap == eicrecon::makeSurfaceStationMap(surfaces, 50.));
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
  CHECK_FALSE(eicrecon::hasTrackStationCoverage(track, stationMap, 3));
  add(58);
  add(9, true);
  add(91, false, true);
  // Four measurements plus an outlier and hole still span only two stations.
  CHECK_FALSE(eicrecon::hasTrackStationCoverage(track, stationMap, 3));

  add(9);
  CHECK(eicrecon::hasTrackStationCoverage(track, stationMap, 3));
  add(91);
  CHECK(eicrecon::hasTrackStationCoverage(track, stationMap, 4));
  add(999); // Unmapped measurements invalidate even otherwise sufficient coverage.
  CHECK_FALSE(eicrecon::hasTrackStationCoverage(track, stationMap, 3));
  CHECK(eicrecon::hasTrackStationCoverage(track, {}, 0));
}

TEST_CASE("Station mapping supports single-plane and paired-face geometries",
          "[tracking][b0stations]") {
  const auto single =
      eicrecon::makeSurfaceStationMap({{1, 6100.}, {2, 6380.}, {3, 6630.}, {4, 6850.}}, 50.);
  const auto paired = eicrecon::makeSurfaceStationMap({{1, 6096.},
                                                       {2, 6104.},
                                                       {3, 6376.},
                                                       {4, 6384.},
                                                       {5, 6626.},
                                                       {6, 6634.},
                                                       {7, 6846.},
                                                       {8, 6854.}},
                                                      50.);
  for (std::uint64_t station = 1; station <= 4; ++station) {
    CHECK(single.at(station) == station - 1);
    CHECK(paired.at(2 * station - 1) == single.at(station));
    CHECK(paired.at(2 * station) == single.at(station));
  }
  CHECK(eicrecon::makeSurfaceStationMap({}, 50.).empty());
}

TEST_CASE("Station mapping rejects invalid geometry and gap settings", "[tracking][b0stations]") {
  for (double gap : {0., -1., std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS(eicrecon::makeSurfaceStationMap({{1, 6100.}}, gap), std::invalid_argument);
  }
  CHECK_THROWS_AS(
      eicrecon::makeSurfaceStationMap({{1, std::numeric_limits<double>::quiet_NaN()}}, 50.),
      std::invalid_argument);
  CHECK_THROWS_AS(eicrecon::makeSurfaceStationMap({{1, 6100.}, {1, 6300.}}, 50.),
                  std::invalid_argument);
}
