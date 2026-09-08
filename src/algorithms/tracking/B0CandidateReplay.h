// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <edm4eic/TrackerHitCollection.h>
#include <podio/ObjectID.h>

namespace eicrecon {

/// One reconstructed hit in an externally supplied B0 candidate.
/// `index` is the RecHit collection order / ObjectID.index in the job that
/// produced the graph. It is rebound to the live collectionID at replay.
struct B0ReplayHit {
  int index{-1};
  std::optional<std::uint64_t> cellID;
};

struct B0ReplayEvent {
  std::int64_t run{};
  std::int64_t event{};
  std::vector<std::vector<B0ReplayHit>> candidates;
};

struct B0ReplayFile {
  int format_version{1};
  std::map<std::pair<std::int64_t, std::int64_t>, B0ReplayEvent> events;
};

B0ReplayFile loadB0ReplayFile(const std::string& path);

/// Bind file indices onto the live RecHits collection. Verifies ObjectID.index
/// and, when present, cellID. Returns ObjectIDs for B0TrackerStubSeeder replay.
std::vector<std::vector<podio::ObjectID>>
bindReplayCandidates(const edm4eic::TrackerHitCollection& hits,
                     const std::vector<std::vector<B0ReplayHit>>& candidates);

} // namespace eicrecon
