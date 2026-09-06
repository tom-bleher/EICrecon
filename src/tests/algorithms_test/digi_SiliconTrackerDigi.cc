// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <catch2/catch_test_macros.hpp>
#include <edm4hep/MCParticleCollection.h>
#include <set>

#include "algorithms/digi/SiliconTrackerDigi.h"

TEST_CASE("Silicon digitizer truth links contain only contributing deposits", "[silicondigi]") {
  eicrecon::SiliconTrackerDigi algorithm("test_silicon_contributors");
  eicrecon::SiliconTrackerDigiConfig config;
  constexpr float threshold = 1.0e-5F;
  config.threshold          = threshold;
  config.timeResolution     = 0.;
  algorithm.applyConfig(config);
  algorithm.init();

  edm4hep::EventHeaderCollection headers;
  headers.create(0, 0);
  edm4hep::MCParticleCollection particles;
  auto particle = particles.create();
  particle.setPDG(2212);
  edm4hep::SimTrackerHitCollection simulated;
  const auto add = [&](uint64_t cell, float energy) {
    auto hit = simulated.create();
    hit.setCellID(cell);
    hit.setEDep(energy);
    hit.setParticle(particle);
  };
  add(1, 2 * threshold);
  add(1, threshold / 2); // Same cell, but contributes neither charge nor truth.
  add(1, 3 * threshold);
  add(2, threshold / 2); // No raw hit should exist for this cell.
  add(3, threshold);     // Exactly at threshold is retained.

  edm4eic::RawTrackerHitCollection raw;
  edm4eic::MCRecoTrackerHitLinkCollection links;
  edm4eic::MCRecoTrackerHitAssociationCollection associations;
  algorithm.process({&headers, &simulated}, {&raw, &links, &associations});

  REQUIRE(raw.size() == 2);
  for (const auto hit : raw) {
    REQUIRE((hit.getCellID() == 1 || hit.getCellID() == 3));
    CHECK(hit.getCharge() == (hit.getCellID() == 1 ? 50 : 10));
  }
  REQUIRE(links.size() == 3);
  REQUIRE(associations.size() == 3);
  std::set<int> linked, associated;
  for (const auto link : links) {
    linked.insert(link.getTo().getObjectID().index);
    CHECK(link.getFrom().getCellID() == link.getTo().getCellID());
    CHECK(link.getTo().getEDep() >= threshold);
  }
  for (const auto association : associations) {
    associated.insert(association.getSimHit().getObjectID().index);
    CHECK(association.getRawHit().getCellID() == association.getSimHit().getCellID());
    CHECK(association.getSimHit().getEDep() >= threshold);
  }
  CHECK(linked == std::set<int>{0, 2, 4});
  CHECK(associated == linked);
}
