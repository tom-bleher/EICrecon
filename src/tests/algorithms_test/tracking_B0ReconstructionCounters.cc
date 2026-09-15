// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "algorithms/tracking/B0ReconstructionCounters.h"

using eicrecon::b0counters::AmbiguityStat;
using eicrecon::b0counters::chainFromAlgorithmName;
using eicrecon::b0counters::CkfStat;
using eicrecon::b0counters::isB0AlgorithmName;
using eicrecon::b0counters::jsonSummary;
using eicrecon::b0counters::reset;
using eicrecon::b0counters::SeederStat;

TEST_CASE("B0 counter names split stub and truth chains", "[B0ReconstructionCounters]") {
  CHECK(chainFromAlgorithmName("B0TrackerCKFTrajectories") == "stub");
  CHECK(chainFromAlgorithmName("B0TrackerCKFTruthSeededTrajectories") == "truth");
  CHECK(chainFromAlgorithmName("B0TrackerAmbiguityResolutionSolver") == "stub");
  CHECK(chainFromAlgorithmName("B0TrackerTruthSeededAmbiguityResolutionSolver") == "truth");
  CHECK(isB0AlgorithmName("B0TrackerCKFTrajectories"));
  CHECK_FALSE(isB0AlgorithmName("CentralCKFTrajectories"));
}

TEST_CASE("B0 counters accumulate without affecting reconstruction state",
          "[B0ReconstructionCounters]") {
  reset();
  eicrecon::b0counters::incrementSeeder(SeederStat::xResidualRejected);
  eicrecon::b0counters::incrementSeeder(SeederStat::overlapRejectedSubset, 2);
  eicrecon::b0counters::incrementCkf("stub", CkfStat::ckfNoTrack);
  eicrecon::b0counters::incrementCkf("truth", CkfStat::ckfTooFewStations, 4);
  eicrecon::b0counters::incrementAmbiguity("stub", AmbiguityStat::ambiguityRejected, 3);

  CHECK(eicrecon::b0counters::seederValue(SeederStat::xResidualRejected) == 1);
  CHECK(eicrecon::b0counters::seederValue(SeederStat::overlapRejectedSubset) == 2);
  CHECK(eicrecon::b0counters::ckfValue("stub", CkfStat::ckfNoTrack) == 1);
  CHECK(eicrecon::b0counters::ckfValue("truth", CkfStat::ckfTooFewStations) == 4);
  CHECK(eicrecon::b0counters::ambiguityValue("stub", AmbiguityStat::ambiguityRejected) == 3);
  CHECK(eicrecon::b0counters::ckfValue("truth", CkfStat::ckfNoTrack) == 0);

  const auto json = jsonSummary();
  CHECK(json.find("\"xResidualRejected\":1") != std::string::npos);
  CHECK(json.find("\"overlapRejected\":2") != std::string::npos);
  CHECK(json.find("\"ckfNoTrack\":1") != std::string::npos);
  CHECK(json.find("\"ckfTooFewStations\":4") != std::string::npos);
  CHECK(json.find("\"ambiguityRejected\":3") != std::string::npos);

  reset();
  CHECK(eicrecon::b0counters::seederValue(SeederStat::xResidualRejected) == 0);
  CHECK(jsonSummary().find("\"xResidualRejected\":0") != std::string::npos);
}

TEST_CASE("B0 counters flush to B0_TRACKING_COUNTERS_FILE", "[B0ReconstructionCounters]") {
  reset();
  eicrecon::b0counters::incrementSeeder(SeederStat::momentumRejected, 7);

  const auto path = std::filesystem::temp_directory_path() / "b0_tracking_counters_test.json";
  std::filesystem::remove(path);
  setenv("B0_TRACKING_COUNTERS_FILE", path.c_str(), 1);
  eicrecon::b0counters::flush();

  std::ifstream in(path);
  REQUIRE(in);
  std::string body((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  CHECK(body.find("\"momentumRejected\":7") != std::string::npos);

  unsetenv("B0_TRACKING_COUNTERS_FILE");
  std::filesystem::remove(path);
  reset();
}