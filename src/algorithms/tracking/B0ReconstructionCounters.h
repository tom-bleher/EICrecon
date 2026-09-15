// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace eicrecon::b0counters {

/// Run-level B0 reconstruction failure accounting.
///
/// Increments are atomic so they are safe across JANA worker threads. They
/// must not change reconstruction output. Dump a JSON snapshot with
/// `jsonSummary()`, or set `B0_TRACKING_COUNTERS_FILE` and call `flush()`
/// (also invoked at process exit).
enum class SeederStat {
  events = 0,
  insufficientStations,
  endpointRejected,
  xResidualRejected,
  yResidualRejected,
  momentumRejected,
  combinationBudgetExceeded,
  overlapRejectedSubset,
  overlapRejectedUnrelated,
  otherRejected,
  _count
};

enum class CkfStat {
  seeds = 0,
  ckfNoTrack,
  ckfTooFewHits,
  ckfTooFewStations,
  smoothingFailure,
  extrapolationFailure,
  tracksAccepted,
  _count
};

enum class AmbiguityStat { tracksIn = 0, tracksOut, ambiguityRejected, _count };

/// "stub" for the data-driven chain, "truth" for the truth-seeded chain.
std::string_view chainFromAlgorithmName(std::string_view name);
bool isB0AlgorithmName(std::string_view name);

void incrementSeeder(SeederStat stat, std::uint64_t n = 1);
void incrementCkf(std::string_view chain, CkfStat stat, std::uint64_t n = 1);
void incrementAmbiguity(std::string_view chain, AmbiguityStat stat, std::uint64_t n = 1);

void reset();
std::uint64_t seederValue(SeederStat stat);
std::uint64_t ckfValue(std::string_view chain, CkfStat stat);
std::uint64_t ambiguityValue(std::string_view chain, AmbiguityStat stat);

/// Compact JSON object with seeder / ckf / ambiguity counts.
std::string jsonSummary();

/// Write JSON to `$B0_TRACKING_COUNTERS_FILE` when set, and emit one info line
/// if any counter is non-zero. Safe to call more than once.
void flush();

} // namespace eicrecon::b0counters
