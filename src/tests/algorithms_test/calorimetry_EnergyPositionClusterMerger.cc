// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/MCRecoClusterParticleAssociationCollection.h>
#include <edm4eic/MCRecoClusterParticleLinkCollection.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/Vector3f.h>
#include <memory>

#include "algorithms/calorimetry/EnergyPositionClusterMerger.h"
#include "algorithms/calorimetry/EnergyPositionClusterMergerConfig.h"

using eicrecon::EnergyPositionClusterMerger;
using eicrecon::EnergyPositionClusterMergerConfig;

TEST_CASE("unmatched energy clusters are dropped or passed through",
          "[EnergyPositionClusterMerger]") {
  const float EPSILON = 1e-5;

  edm4hep::MCParticleCollection mcparts;
  auto mc1 = mcparts.create();
  auto mc2 = mcparts.create();

  // energy clusters: one matching the position cluster in eta/phi/energy, one far from it
  edm4eic::ClusterCollection energy_clus;
  auto ec_matched = energy_clus.create();
  ec_matched.setEnergy(5.0 * edm4eic::unit::GeV);
  ec_matched.setPosition(edm4hep::Vector3f{800, 0, 0});
  auto ec_alone = energy_clus.create();
  ec_alone.setEnergy(3.0 * edm4eic::unit::GeV);
  ec_alone.setPosition(edm4hep::Vector3f{0, 800, -1500});
  ec_alone.setNhits(7);

  edm4eic::MCRecoClusterParticleAssociationCollection energy_assoc;
  auto ea1 = energy_assoc.create();
  ea1.setRec(ec_matched);
  ea1.setSim(mc1);
  ea1.setWeight(1.0);
  auto ea2 = energy_assoc.create();
  ea2.setRec(ec_alone);
  ea2.setSim(mc2);
  ea2.setWeight(0.8);

  edm4eic::ClusterCollection pos_clus;
  auto pc = pos_clus.create();
  pc.setEnergy(4.5 * edm4eic::unit::GeV);
  pc.setPosition(edm4hep::Vector3f{810, 5, 0});
  edm4eic::MCRecoClusterParticleAssociationCollection pos_assoc;

  for (bool keep : {false, true}) {
    EnergyPositionClusterMerger algo("EnergyPositionClusterMerger");
    EnergyPositionClusterMergerConfig cfg;
    cfg.keepUnmatchedEnergyClusters = keep;
    algo.applyConfig(cfg);
    algo.init();

    auto merged       = std::make_unique<edm4eic::ClusterCollection>();
    auto merged_links = std::make_unique<edm4eic::MCRecoClusterParticleLinkCollection>();
    auto merged_assoc = std::make_unique<edm4eic::MCRecoClusterParticleAssociationCollection>();
    algo.process({&energy_clus, &energy_assoc, &pos_clus, &pos_assoc},
                 {merged.get(), merged_links.get(), merged_assoc.get()});

    // the matched pair: energy from the energy cluster, position from the position cluster
    REQUIRE(merged->size() == (keep ? 2 : 1));
    REQUIRE_THAT((*merged)[0].getEnergy(), Catch::Matchers::WithinAbs(5.0, EPSILON));
    REQUIRE_THAT((*merged)[0].getPosition().x, Catch::Matchers::WithinAbs(810, EPSILON));
    if (keep) {
      // the unmatched energy cluster, with its own energy, position, hit count and association
      const auto& out = (*merged)[1];
      REQUIRE_THAT(out.getEnergy(), Catch::Matchers::WithinAbs(3.0, EPSILON));
      REQUIRE_THAT(out.getPosition().z, Catch::Matchers::WithinAbs(-1500, EPSILON));
      REQUIRE(out.getNhits() == 7);
      REQUIRE(out.getClusters().size() == 1);
      REQUIRE(merged_assoc->size() == 2);
      REQUIRE((*merged_assoc)[1].getRec() == out);
      REQUIRE((*merged_assoc)[1].getSim() == mc2);
      REQUIRE_THAT((*merged_assoc)[1].getWeight(), Catch::Matchers::WithinAbs(0.8, EPSILON));
    } else {
      REQUIRE(merged_assoc->size() == 1);
    }
  }
}
