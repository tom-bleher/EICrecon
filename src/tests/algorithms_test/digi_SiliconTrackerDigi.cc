// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Evaluator/DD4hepUnits.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <cstdint>

#include "algorithms/digi/SiliconTrackerDigi.h"
#include "algorithms/digi/SiliconTrackerDigiConfig.h"

using Catch::Approx;
using eicrecon::SiliconTrackerDigi;
using eicrecon::SiliconTrackerDigiConfig;

namespace {

void addSimHit(edm4hep::SimTrackerHitCollection& sim_hits,
               edm4hep::MCParticleCollection& particles, std::uint64_t cell_id, double e_dep,
               double time) {
  auto particle = particles.create();
  particle.setPDG(2212);
  particle.setCharge(1.0F);
  particle.setMass(0.938272F);
  particle.setTime(0.0F);

  auto hit = sim_hits.create();
  hit.setCellID(cell_id);
  hit.setEDep(e_dep);
  hit.setTime(time);
  hit.setParticle(particle);
}

struct DigiOutput {
  edm4eic::RawTrackerHitCollection raw_hits;
  edm4eic::MCRecoTrackerHitLinkCollection links;
  edm4eic::MCRecoTrackerHitAssociationCollection associations;
};

DigiOutput runDigi(const edm4hep::EventHeaderCollection& headers,
                   const edm4hep::SimTrackerHitCollection& sim_hits, double threshold) {
  SiliconTrackerDigi algo("silicon_tracker_digi_test");
  SiliconTrackerDigiConfig cfg;
  cfg.threshold      = threshold;
  cfg.timeResolution = 0.0;
  algo.applyConfig(cfg);
  algo.init();

  DigiOutput output;
  algo.process({&headers, &sim_hits}, {&output.raw_hits, &output.links, &output.associations});
  return output;
}

} // namespace

TEST_CASE("SiliconTrackerDigi thresholds the summed channel signal", "[SiliconTrackerDigi]") {
  edm4hep::EventHeaderCollection headers;
  headers.create(1, 0);
  edm4hep::SimTrackerHitCollection sim_hits;
  edm4hep::MCParticleCollection particles;

  constexpr std::uint64_t cell_id = 0x1234;
  addSimHit(sim_hits, particles, cell_id, 6.0 * dd4hep::keV, 10.0);
  addSimHit(sim_hits, particles, cell_id, 6.0 * dd4hep::keV, 10.1);

  auto output = runDigi(headers, sim_hits, 10.0 * dd4hep::keV);

  REQUIRE(output.raw_hits.size() == 1);
  CHECK(output.raw_hits[0].getCellID() == cell_id);
  CHECK(output.raw_hits[0].getCharge() == 12);
  REQUIRE(output.associations.size() == 2);
  CHECK(output.associations[0].getWeight() == Approx(0.5));
  CHECK(output.associations[1].getWeight() == Approx(0.5));
}

TEST_CASE("SiliconTrackerDigi keeps truth weights consistent with accumulated charge",
          "[SiliconTrackerDigi]") {
  edm4hep::EventHeaderCollection headers;
  headers.create(2, 0);
  edm4hep::SimTrackerHitCollection sim_hits;
  edm4hep::MCParticleCollection particles;

  constexpr std::uint64_t cell_id = 0x5678;
  addSimHit(sim_hits, particles, cell_id, 12.0 * dd4hep::keV, 7.0);
  addSimHit(sim_hits, particles, cell_id, 6.0 * dd4hep::keV, 7.1);

  auto output = runDigi(headers, sim_hits, 10.0 * dd4hep::keV);

  REQUIRE(output.raw_hits.size() == 1);
  CHECK(output.raw_hits[0].getCharge() == 18);
  REQUIRE(output.links.size() == 2);
  REQUIRE(output.associations.size() == 2);

  const double weight_sum = output.associations[0].getWeight() + output.associations[1].getWeight();
  CHECK(weight_sum == Approx(1.0));
  CHECK(output.associations[0].getWeight() == Approx(2.0 / 3.0));
  CHECK(output.associations[1].getWeight() == Approx(1.0 / 3.0));
}

TEST_CASE("SiliconTrackerDigi drops a channel whose summed signal is below threshold",
          "[SiliconTrackerDigi]") {
  edm4hep::EventHeaderCollection headers;
  headers.create(3, 0);
  edm4hep::SimTrackerHitCollection sim_hits;
  edm4hep::MCParticleCollection particles;

  constexpr std::uint64_t cell_id = 0x9abc;
  addSimHit(sim_hits, particles, cell_id, 4.0 * dd4hep::keV, 3.0);
  addSimHit(sim_hits, particles, cell_id, 5.0 * dd4hep::keV, 3.1);

  auto output = runDigi(headers, sim_hits, 10.0 * dd4hep::keV);

  CHECK(output.raw_hits.empty());
  CHECK(output.links.empty());
  CHECK(output.associations.empty());
}

TEST_CASE("SiliconTrackerDigi preserves single-hit behavior", "[SiliconTrackerDigi]") {
  edm4hep::EventHeaderCollection headers;
  headers.create(4, 0);
  edm4hep::SimTrackerHitCollection sim_hits;
  edm4hep::MCParticleCollection particles;

  constexpr std::uint64_t cell_id = 0xdef0;
  addSimHit(sim_hits, particles, cell_id, 12.0 * dd4hep::keV, 5.0);

  auto output = runDigi(headers, sim_hits, 10.0 * dd4hep::keV);

  REQUIRE(output.raw_hits.size() == 1);
  CHECK(output.raw_hits[0].getCharge() == 12);
  REQUIRE(output.associations.size() == 1);
  CHECK(output.associations[0].getWeight() == Approx(1.0));
}
