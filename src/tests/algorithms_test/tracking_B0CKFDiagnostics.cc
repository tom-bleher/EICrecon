// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/EventData/ProxyAccessor.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <catch2/catch_test_macros.hpp>

#include "algorithms/tracking/B0ReconstructionCounters.h"
#include "algorithms/tracking/CKFTrackingConfig.h"

namespace ckfdiag = eicrecon::b0counters::ckfdiag;

TEST_CASE("ckfdiag resolveRow keeps accepted rows only when column exists",
          "[tracking][b0ckfdiag]") {
  // No column (central tracking, older files): every row resolves, as before.
  CHECK(ckfdiag::resolveRow(false, ckfdiag::kAccepted));
  CHECK(ckfdiag::resolveRow(false, ckfdiag::kTooFewStations));
  CHECK(ckfdiag::resolveRow(false, ckfdiag::kFindFailed));
  // With column: only accepted candidates resolve.
  CHECK(ckfdiag::resolveRow(true, ckfdiag::kAccepted));
  CHECK_FALSE(ckfdiag::resolveRow(true, ckfdiag::kNoValidMeasurement));
  CHECK_FALSE(ckfdiag::resolveRow(true, ckfdiag::kTooFewHits));
  CHECK_FALSE(ckfdiag::resolveRow(true, ckfdiag::kTooFewStations));
  CHECK_FALSE(ckfdiag::resolveRow(true, ckfdiag::kSmoothingFailed));
  CHECK_FALSE(ckfdiag::resolveRow(true, ckfdiag::kExtrapolationFailed));
  CHECK_FALSE(ckfdiag::resolveRow(true, ckfdiag::kFindFailed));
  CHECK_FALSE(ckfdiag::resolveRow(true, ckfdiag::kNoCandidates));
}

TEST_CASE("ckfdiag diagnostics are opt-in and off by default",
          "[tracking][b0ckfdiag]") {
  CHECK(eicrecon::CKFTrackingConfig{}.keepB0CKFDiagnostics == false);
}

namespace {

/// Mutable container with a ckf_status column holding the given statuses.
ActsExamples::TrackContainer makeStatusContainer(const std::vector<unsigned>& statuses) {
  auto tracks = std::make_shared<Acts::VectorTrackContainer>();
  auto states = std::make_shared<Acts::VectorMultiTrajectory>();
  ActsExamples::TrackContainer container{tracks, states};
  container.addColumn<unsigned int>(ckfdiag::kStatusColumn);
  Acts::ProxyAccessor<unsigned int> writer(ckfdiag::kStatusColumn);
  for (const unsigned status : statuses) {
    auto track = container.makeTrack();
    writer(track) = status;
  }
  return container;
}

} // namespace

TEST_CASE("ckfdiag column detection works on real containers, including empty",
          "[tracking][b0ckfdiag]") {
  ActsExamples::TrackContainer plain{
      std::make_shared<Acts::VectorTrackContainer>(),
      std::make_shared<Acts::VectorMultiTrajectory>()};
  CHECK_FALSE(plain.hasColumn(ckfdiag::kStatusColumn));

  auto mixed = makeStatusContainer(
      {ckfdiag::kAccepted, ckfdiag::kTooFewStations, ckfdiag::kAccepted,
       ckfdiag::kFindFailed, ckfdiag::kTooFewHits, ckfdiag::kAccepted});
  CHECK(mixed.hasColumn(ckfdiag::kStatusColumn));

  ActsExamples::TrackContainer empty{
      std::make_shared<Acts::VectorTrackContainer>(),
      std::make_shared<Acts::VectorMultiTrajectory>()};
  empty.addColumn<unsigned int>(ckfdiag::kStatusColumn);
  CHECK(empty.hasColumn(ckfdiag::kStatusColumn));
  CHECK(empty.size() == 0);
}

TEST_CASE("mixed accepted/rejected/marker rows resolve exactly the accepted set",
          "[tracking][b0ckfdiag]") {
  auto mixed = makeStatusContainer(
      {ckfdiag::kAccepted, ckfdiag::kTooFewStations, ckfdiag::kAccepted,
       ckfdiag::kFindFailed, ckfdiag::kNoCandidates, ckfdiag::kAccepted});
  Acts::ConstProxyAccessor<unsigned int> reader(ckfdiag::kStatusColumn);
  std::vector<std::size_t> selected;
  std::size_t index = 0;
  for (const auto& track : mixed) {
    if (ckfdiag::resolveRow(mixed.hasColumn(ckfdiag::kStatusColumn),
                            reader(track))) {
      selected.push_back(index);
    }
    ++index;
  }
  CHECK(selected == std::vector<std::size_t>{0, 2, 5});

  // An accepted-only container resolves every row, matching the pre-diagnostic
  // behavior for the same physics content.
  auto acceptedOnly = makeStatusContainer(
      {ckfdiag::kAccepted, ckfdiag::kAccepted, ckfdiag::kAccepted});
  std::size_t kept = 0;
  for (const auto& track : acceptedOnly) {
    if (ckfdiag::resolveRow(acceptedOnly.hasColumn(ckfdiag::kStatusColumn),
                            reader(track))) {
      ++kept;
    }
  }
  CHECK(kept == acceptedOnly.size());
  CHECK(kept == selected.size());
}
