// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0ReconstructionCounters.h"

#include <spdlog/spdlog.h>
#include <array>
#include <atomic>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>

namespace eicrecon::b0counters {
namespace {

  constexpr std::size_t kSeederN    = static_cast<std::size_t>(SeederStat::_count);
  constexpr std::size_t kCkfN       = static_cast<std::size_t>(CkfStat::_count);
  constexpr std::size_t kAmbiguityN = static_cast<std::size_t>(AmbiguityStat::_count);

  struct ChainCounts {
    std::array<std::atomic<std::uint64_t>, kCkfN> ckf{};
    std::array<std::atomic<std::uint64_t>, kAmbiguityN> ambiguity{};
  };

  std::array<std::atomic<std::uint64_t>, kSeederN> g_seeder{};
  ChainCounts g_stub;
  ChainCounts g_truth;

  ChainCounts& chainCounts(std::string_view chain) { return chain == "truth" ? g_truth : g_stub; }

  void resetArray(auto& values) {
    for (auto& value : values) {
      value.store(0, std::memory_order_relaxed);
    }
  }

  std::uint64_t load(const auto& values, std::size_t index) {
    return values[index].load(std::memory_order_relaxed);
  }

  void add(auto& values, std::size_t index, std::uint64_t n) {
    values[index].fetch_add(n, std::memory_order_relaxed);
  }

  bool anyNonZero() {
    for (std::size_t i = 0; i < kSeederN; ++i) {
      if (load(g_seeder, i) != 0) {
        return true;
      }
    }
    const auto checkChain = [](const ChainCounts& chain) {
      for (std::size_t i = 0; i < kCkfN; ++i) {
        if (load(chain.ckf, i) != 0) {
          return true;
        }
      }
      for (std::size_t i = 0; i < kAmbiguityN; ++i) {
        if (load(chain.ambiguity, i) != 0) {
          return true;
        }
      }
      return false;
    };
    return checkChain(g_stub) || checkChain(g_truth);
  }

  void appendObject(std::ostringstream& out, const auto& names, const auto& values) {
    out << '{';
    for (std::size_t i = 0; i < names.size(); ++i) {
      if (i > 0) {
        out << ',';
      }
      out << '"' << names[i] << '"' << ':' << load(values, i);
    }
    out << '}';
  }

  struct FlushAtExit {
    ~FlushAtExit() { flush(); }
  };

  // Register process-exit dump. Reconstruction never reads these counters.
  const FlushAtExit kFlushAtExit;

} // namespace

std::string_view chainFromAlgorithmName(std::string_view name) {
  return name.find("TruthSeeded") != std::string_view::npos ? "truth" : "stub";
}

bool isB0AlgorithmName(std::string_view name) {
  return name.find("B0Tracker") != std::string_view::npos;
}

void incrementSeeder(SeederStat stat, std::uint64_t n) {
  add(g_seeder, static_cast<std::size_t>(stat), n);
}

void incrementCkf(std::string_view chain, CkfStat stat, std::uint64_t n) {
  add(chainCounts(chain).ckf, static_cast<std::size_t>(stat), n);
}

void incrementAmbiguity(std::string_view chain, AmbiguityStat stat, std::uint64_t n) {
  add(chainCounts(chain).ambiguity, static_cast<std::size_t>(stat), n);
}

void reset() {
  resetArray(g_seeder);
  resetArray(g_stub.ckf);
  resetArray(g_stub.ambiguity);
  resetArray(g_truth.ckf);
  resetArray(g_truth.ambiguity);
}

std::uint64_t seederValue(SeederStat stat) {
  return load(g_seeder, static_cast<std::size_t>(stat));
}

std::uint64_t ckfValue(std::string_view chain, CkfStat stat) {
  return load(chainCounts(chain).ckf, static_cast<std::size_t>(stat));
}

std::uint64_t ambiguityValue(std::string_view chain, AmbiguityStat stat) {
  return load(chainCounts(chain).ambiguity, static_cast<std::size_t>(stat));
}

std::string jsonSummary() {
  static constexpr std::array<std::string_view, kSeederN> kSeederNames{
      "events",
      "insufficientStations",
      "endpointRejected",
      "xResidualRejected",
      "yResidualRejected",
      "momentumRejected",
      "combinationBudgetExceeded",
      "overlapRejectedSubset",
      "overlapRejectedUnrelated",
      "otherRejected",
  };
  static constexpr std::array<std::string_view, kCkfN> kCkfNames{
      "seeds",
      "ckfNoTrack",
      "ckfTooFewHits",
      "ckfTooFewStations",
      "smoothingFailure",
      "extrapolationFailure",
      "tracksAccepted",
  };
  static constexpr std::array<std::string_view, kAmbiguityN> kAmbiguityNames{
      "tracksIn",
      "tracksOut",
      "ambiguityRejected",
  };

  std::ostringstream out;
  out << "{\"seeder\":";
  appendObject(out, kSeederNames, g_seeder);
  const auto appendChain = [&](std::string_view label, const ChainCounts& chain) {
    out << '"' << label << "\":{\"ckf\":";
    appendObject(out, kCkfNames, chain.ckf);
    out << ",\"ambiguity\":";
    appendObject(out, kAmbiguityNames, chain.ambiguity);
    out << '}';
  };
  out << ",\"overlapRejected\":"
      << (seederValue(SeederStat::overlapRejectedSubset) +
          seederValue(SeederStat::overlapRejectedUnrelated));
  out << ",\"chains\":{";
  appendChain("stub", g_stub);
  out << ',';
  appendChain("truth", g_truth);
  out << "}}";
  return out.str();
}

void flush() {
  if (!anyNonZero()) {
    return;
  }
  const std::string json = jsonSummary();
  if (const char* path = std::getenv("B0_TRACKING_COUNTERS_FILE");
      path != nullptr && *path != '\0') {
    std::ofstream file(path);
    if (file) {
      file << json << '\n';
    }
  }
  if (auto logger = spdlog::default_logger()) {
    logger->info("B0TrackingCounters {}", json);
  }
}

} // namespace eicrecon::b0counters
