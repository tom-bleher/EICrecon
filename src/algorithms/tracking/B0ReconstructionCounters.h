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
  // findTracks outcomes (appended; existing ordinals are stable).
  ckfFindFailed,
  ckfFindEmpty,
  ckfFindErrCkf,
  ckfFindErrPropagation,
  _count
};

enum class AmbiguityStat { tracksIn = 0, tracksOut, ambiguityRejected, _count };

/// Per-seed CKF failure diagnostics carried inside the UNFILTERED ACTS track
/// containers for B0 chains (see CKFTracking::process).
///
/// Every processed seed owns at least one row: accepted and CKF-rejected
/// candidates keep their fitted states, while seeds with zero candidates
/// leave a stateless marker row recording the findTracks outcome. The status
/// column lets downstream stages (AmbiguitySolver, ActsToTracks, B0Trackers)
/// separate finding failures from disambiguation without changing physics
/// output. Columns exist only when CKFTracking runs with numB0StationsMin > 0,
/// so central tracking containers are unaffected; consumers must probe for the
/// column and keep every row when it is absent (older files, unit tests).
namespace ckfdiag {

// Row outcome in the "ckf_status" dynamic column.
inline constexpr unsigned kAccepted = 0;
inline constexpr unsigned kNoValidMeasurement = 1;
inline constexpr unsigned kTooFewHits = 2;
inline constexpr unsigned kTooFewStations = 3;
inline constexpr unsigned kSmoothingFailed = 4;
inline constexpr unsigned kExtrapolationFailed = 5;
// Stateless markers (seeds with zero candidates).
inline constexpr unsigned kFindFailed = 6;
inline constexpr unsigned kNoCandidates = 7;

// findTracks error classes in the "ckf_find_err" column (class * 1000 + value).
// kFindErrCkf covers errors raised by the CKF actor itself (update, measurement
// selection, max steps); kFindErrPropagation covers everything propagated from
// navigation/propagation (stepper, navigator, material).
inline constexpr unsigned kFindErrNone = 0;
inline constexpr unsigned kFindErrCkf = 1;
inline constexpr unsigned kFindErrPropagation = 2;

inline constexpr const char kStatusColumn[] = "ckf_status";
inline constexpr const char kFindErrColumn[] = "ckf_find_err";

inline unsigned packFindError(unsigned errorClass, unsigned errorValue) {
  return errorClass * 1000u + (errorValue % 1000u);
}

inline unsigned findErrorClass(unsigned packed) { return packed / 1000u; }

inline unsigned findErrorValue(unsigned packed) { return packed % 1000u; }

// Whether a row enters ambiguity resolution / EDM conversion: every row when
// the container carries no ckf_status column (central tracking, older files),
// otherwise only accepted candidates. Unit-tested in tracking_B0CKFDiagnostics.
inline bool resolveRow(bool haveStatusColumn, unsigned status) {
  return !haveStatusColumn || status == kAccepted;
}

} // namespace ckfdiag

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
