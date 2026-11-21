// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/Vector3f.h>
#include <podio/RelationRange.h>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <vector>

#include "algorithms/tracking/TruthTrackSeeding.h"
#include "algorithms/tracking/TruthTrackSeedingConfig.h"

namespace {

const double kFieldStrength = 1.3 * Acts::UnitConstants::T;
const Acts::Vector3 kFieldDirection(std::sin(0.5), std::cos(0.5), 0.);
const std::vector<double> kStations{5850., 6272., 6601., 6880.}; // mm

/// An event: reconstructed hits, their raw hits, the sim hits and particles behind them, and the
/// raw hit links that the digitization makes
struct Event {
  edm4hep::MCParticleCollection particles;
  edm4hep::SimTrackerHitCollection simHits;
  edm4eic::RawTrackerHitCollection rawHits;
  edm4eic::TrackerHitCollection hits;
  edm4eic::MCRecoTrackerHitLinkCollection links;

  /// Point of a helix in the field, starting at z = 5800 mm along z and 10 mrad along the field
  static Acts::Vector3 helixAt(double q, double p, double z) {
    const Acts::Vector3 across = Acts::Vector3::UnitZ().cross(kFieldDirection);
    const double radius        = p / kFieldStrength;
    const double dz            = z - 5800.;
    return q * (radius - std::sqrt(radius * radius - dz * dz)) * across +
           0.01 * dz * kFieldDirection + z * Acts::Vector3::UnitZ();
  }

  /// A reconstructed hit at a position, with its raw hit
  edm4eic::TrackerHit addHit(const Acts::Vector3& position, float time) {
    auto raw = rawHits.create();
    auto hit = hits.create(0ULL,
                           edm4hep::Vector3f{static_cast<float>(position.x()),
                                             static_cast<float>(position.y()),
                                             static_cast<float>(position.z())},
                           edm4eic::CovDiag3f{4.e-4F, 4.e-4F, 0.F}, time, 0.F, 0.F, 0.F);
    hit.setRawHit(raw);
    return hit;
  }

  /// Links a reconstructed hit to a particle through a new sim hit
  void link(const edm4eic::TrackerHit& hit, const edm4hep::MCParticle& particle, float time,
            bool bySecondary = false) {
    auto sim = simHits.create();
    sim.setParticle(particle);
    sim.setTime(time);
    sim.setProducedBySecondary(bySecondary);
    auto lnk = links.create();
    lnk.setFrom(hit.getRawHit());
    lnk.setTo(sim);
  }

  /// A particle of charge q and momentum p with hits on the given stations (indices of kStations)
  edm4hep::MCParticle addTrack(double q, double p, std::initializer_list<std::size_t> stations) {
    auto particle = particles.create();
    for (const auto s : stations) {
      const auto time = static_cast<float>(s + 1);
      link(addHit(helixAt(q, p, kStations[s]), time), particle, time);
    }
    return particle;
  }
};

/// The seeds of an event
struct Seeds {
  edm4eic::TrackSeedCollection seeds;
  edm4eic::TrackParametersCollection parameters;
};

Seeds findSeeds(const Event& event, const eicrecon::TruthTrackSeedingConfig& cfg = {
                                        .minHitDistance = 150. * Acts::UnitConstants::mm}) {
  const Acts::ConstantBField field(kFieldStrength * kFieldDirection);
  eicrecon::TruthTrackSeeding algo("TruthTrackSeeding");
  algo.applyConfig(cfg);
  algo.init();
  Seeds out;
  algo.findSeeds(event.hits, event.links, field, Acts::MagneticFieldContext{}, out.seeds,
                 out.parameters);
  return out;
}

} // namespace

// The test environment has no Acts magnetic field service, so the tests call findSeeds() with an
// Acts::ConstantBField instead of process().

TEST_CASE("TruthTrackSeeding seeds every triplet of a particle's hits", "[TruthTrackSeeding]") {
  for (const double q : {+1., -1.}) {
    for (const double p : {5., 40., 100.}) {
      Event event;
      event.addTrack(q, p, {0, 1, 2, 3});
      const auto out = findSeeds(event);

      // Stations 1-2-3, 1-2-4, 1-3-4 and 2-3-4
      REQUIRE(out.seeds.size() == 4);
      for (const auto seed : out.seeds) {
        CHECK(seed.getHits().size() == 3);
        CHECK(seed.getParams().getQOverP() == Catch::Approx(q / p).epsilon(0.02));
      }
    }
  }
}

TEST_CASE("TruthTrackSeeding needs three hits apart", "[TruthTrackSeeding]") {
  SECTION("fewer than three stations") {
    Event event;
    event.addTrack(1., 40., {0, 3});
    CHECK(findSeeds(event).seeds.empty());
  }

  SECTION("a second hit on a station is skipped") {
    Event event;
    const auto particle = event.addTrack(1., 40., {0, 1, 2});
    // An overlapping sensor, 1 mm downstream of the station 2 hit and after it in time
    event.link(event.addHit(Event::helixAt(1., 40., kStations[1] + 1.), 2.5F), particle, 2.5F);
    const auto out = findSeeds(event);
    REQUIRE(out.seeds.size() == 1);
    for (const auto& hit : out.seeds[0].getHits()) {
      CHECK(hit.getTime() != Catch::Approx(2.5F));
    }
  }

  SECTION("hits beyond maxHitsPerParticle are not used") {
    Event event;
    event.addTrack(1., 40., {0, 1, 2, 3});
    CHECK(findSeeds(event, {.minHitDistance = 150., .maxHitsPerParticle = 3}).seeds.size() == 1);
  }
}

TEST_CASE("TruthTrackSeeding takes only a particle's own hits", "[TruthTrackSeeding]") {
  SECTION("a hit made by a secondary does not count") {
    Event event;
    const auto particle = event.addTrack(1., 40., {0, 1});
    event.link(event.addHit(Event::helixAt(1., 40., kStations[2]), 3.F), particle, 3.F, true);
    CHECK(findSeeds(event).seeds.empty());
  }

  SECTION("a hit shared by two particles seeds both") {
    Event event;
    event.addTrack(1., 40., {0, 1, 2});
    // A second particle shares the three hits of the first and has its own on station 4
    auto second = event.particles.create();
    for (const auto hit : event.hits) {
      event.link(hit, second, hit.getTime());
    }
    event.link(event.addHit(Event::helixAt(1., 40., kStations[3]), 4.F), second, 4.F);
    // The first particle: stations 1-2-3; the second: the four triplets of stations 1 to 4
    CHECK(findSeeds(event).seeds.size() == 5);
  }

  SECTION("hits without a particle are ignored") {
    Event event;
    for (const auto s : {0, 1, 2}) {
      auto hit = event.addHit(Event::helixAt(1., 40., kStations[s]), 1.F);
      auto sim = event.simHits.create();
      auto lnk = event.links.create();
      lnk.setFrom(hit.getRawHit());
      lnk.setTo(sim);
    }
    CHECK(findSeeds(event).seeds.empty());
  }
}
