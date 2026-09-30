// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <DD4hep/DD4hepUnits.h>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <initializer_list>

#include "algorithms/digi/SiliconTrackerDigi.h"
#include "algorithms/digi/SiliconTrackerDigiConfig.h"

namespace {

// Digitize deposits (keV) that all fall into one cell
edm4eic::RawTrackerHitCollection digitize(const eicrecon::SiliconTrackerDigiConfig& cfg,
                                          std::initializer_list<double> edeps_keV) {
  eicrecon::SiliconTrackerDigi algo("SiliconTrackerDigi");
  algo.applyConfig(cfg);
  algo.init();

  edm4hep::EventHeaderCollection headers;
  headers.create(1, 0);
  edm4hep::MCParticleCollection particles;
  auto particle = particles.create();
  edm4hep::SimTrackerHitCollection sim_hits;
  for (const double edep : edeps_keV) {
    auto hit = sim_hits.create();
    hit.setCellID(42);
    hit.setEDep(edep * dd4hep::keV);
    hit.setParticle(particle);
  }

  edm4eic::RawTrackerHitCollection raw_hits;
  edm4eic::MCRecoTrackerHitLinkCollection links;
  edm4eic::MCRecoTrackerHitAssociationCollection associations;
  algo.process({&headers, &sim_hits}, {&raw_hits, &links, &associations});
  return raw_hits;
}

} // namespace

TEST_CASE("threshold per SimTrackerHit by default", "[SiliconTrackerDigi]") {
  const eicrecon::SiliconTrackerDigiConfig cfg{.threshold = 10 * dd4hep::keV};
  const auto raw_hits = digitize(cfg, {12, 9});
  REQUIRE(raw_hits.size() == 1);
  CHECK(raw_hits[0].getCharge() == 12);
}

TEST_CASE("threshold on the summed cell energy", "[SiliconTrackerDigi]") {
  const eicrecon::SiliconTrackerDigiConfig cfg{.threshold          = 2 * dd4hep::keV,
                                               .thresholdOnCellSum = true};

  SECTION("two small shares of one cell pass together") {
    const auto raw_hits = digitize(cfg, {1.2, 1.4});
    REQUIRE(raw_hits.size() == 1);
    CHECK(raw_hits[0].getCharge() == 3);
  }

  SECTION("a cell below threshold is dropped") { CHECK(digitize(cfg, {0.8, 0.9}).empty()); }
}

TEST_CASE("finer charge counts", "[SiliconTrackerDigi]") {
  const eicrecon::SiliconTrackerDigiConfig cfg{
      .threshold = 2 * dd4hep::keV, .thresholdOnCellSum = true, .countsPerKeV = 10};
  const auto raw_hits = digitize(cfg, {1.24, 1.37});
  REQUIRE(raw_hits.size() == 1);
  CHECK(raw_hits[0].getCharge() == 26);
}
