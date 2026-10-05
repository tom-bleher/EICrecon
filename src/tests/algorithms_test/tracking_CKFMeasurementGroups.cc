// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/EventData/TrackStateType.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <catch2/catch_test_macros.hpp>
#include <array>
#include <limits>
#include <memory>

#include "algorithms/tracking/CKFTracking.h"

namespace {
void addState(ActsExamples::TrackProxy& track, double z, bool outlier = false, bool hole = false) {
  auto transform              = Acts::Transform3::Identity();
  transform.translation().z() = z;
  auto state                  = track.appendTrackState(Acts::TrackStatePropMask::None);
  state.setReferenceSurface(Acts::Surface::makeShared<Acts::PlaneSurface>(transform));
#if Acts_VERSION_MAJOR >= 45
  state.typeFlags().setIsMeasurement();
  if (outlier) {
    state.typeFlags().setIsOutlier();
  }
  if (hole) {
    state.typeFlags().setIsHole();
  }
#else
  state.typeFlags().set(Acts::TrackStateFlag::MeasurementFlag);
  if (outlier) {
    state.typeFlags().set(Acts::TrackStateFlag::OutlierFlag);
  }
  if (hole) {
    state.typeFlags().set(Acts::TrackStateFlag::HoleFlag);
  }
#endif
}
} // namespace

TEST_CASE("CKF measurement groups count stations rather than their two faces", "[B0Tracking]") {
#if Acts_VERSION_MAJOR >= 45
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  const Acts::GeometryContext gctx;
#endif
  auto tracks = std::make_shared<Acts::VectorTrackContainer>();
  auto states = std::make_shared<Acts::VectorMultiTrajectory>();
  ActsExamples::TrackContainer container(tracks, states);
  auto track = container.makeTrack();
  CHECK(eicrecon::CKFTracking::countMeasurementGroups(track, gctx, 50.) == 0);
  for (const double z : {5897., 5903., 6197., 6203.}) {
    addState(track, z);
  }
  CHECK(eicrecon::CKFTracking::countMeasurementGroups(track, gctx, 50.) == 2);
  addState(track, 6497., true);
  addState(track, 6503., false, true);
  CHECK(eicrecon::CKFTracking::countMeasurementGroups(track, gctx, 50.) == 2);
  addState(track, 6501.);
  CHECK(eicrecon::CKFTracking::countMeasurementGroups(track, gctx, 50.) == 3);

  auto reversed = container.makeTrack();
  for (const double z : {6503., 6497., 6203., 6197., 5903., 5897.}) {
    addState(reversed, z);
  }
  CHECK(eicrecon::CKFTracking::countMeasurementGroups(reversed, gctx, 50.) == 3);
  auto missing = container.makeTrack();
  for (const double z : {5900., 6500., 6800.}) {
    addState(missing, z);
  }
  CHECK(eicrecon::CKFTracking::countMeasurementGroups(missing, gctx, 50.) == 3);
}

TEST_CASE("CKF measurement groups reject nonfinite surface positions", "[B0Tracking]") {
#if Acts_VERSION_MAJOR >= 45
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  const Acts::GeometryContext gctx;
#endif
  auto tracks = std::make_shared<Acts::VectorTrackContainer>();
  auto states = std::make_shared<Acts::VectorMultiTrajectory>();
  ActsExamples::TrackContainer container(tracks, states);
  auto track = container.makeTrack();
  addState(track, 5900.);
  addState(track, std::numeric_limits<double>::quiet_NaN());
  addState(track, 6500.);
  CHECK(eicrecon::CKFTracking::countMeasurementGroups(track, gctx, 50.) == 0);
}
