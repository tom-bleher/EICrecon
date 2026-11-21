// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "TruthTrackSeeding.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Result.hpp>
#include <edm4eic/unit_system.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Vector3f.h>
#include <podio/ObjectID.h>
#include <algorithm>
#include <cstdint>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "algorithms/tracking/HelixSeed.h"

namespace eicrecon {

namespace {

  /// Orders ObjectIDs by collection and index, which older podio does not define
  struct ObjectIDLess {
    bool operator()(const podio::ObjectID& a, const podio::ObjectID& b) const {
      return std::tie(a.collectionID, a.index) < std::tie(b.collectionID, b.index);
    }
  };

} // namespace

void TruthTrackSeeding::process(const Input& input, const Output& output) const {
  const auto [hits, rawHitLinks] = input;
  auto [seeds, parameters]       = output;
  findSeeds(*hits, *rawHitLinks, *m_geoSvc->getFieldProvider(),
            m_geoSvc->getActsMagneticFieldContext(), *seeds, *parameters);
  debug("Created {} truth track seeds from {} hits", seeds->size(), hits->size());
}

void TruthTrackSeeding::findSeeds(const edm4eic::TrackerHitCollection& hits,
                                  const edm4eic::MCRecoTrackerHitLinkCollection& rawHitLinks,
                                  const Acts::MagneticFieldProvider& fieldProvider,
                                  const Acts::MagneticFieldContext& fieldContext,
                                  edm4eic::TrackSeedCollection& seeds,
                                  edm4eic::TrackParametersCollection& parameters) const {
  // The particles of each raw hit, with the time of their sim hit there
  std::map<podio::ObjectID, std::vector<std::pair<float, podio::ObjectID>>, ObjectIDLess>
      particlesOfRawHit;
  for (const auto& link : rawHitLinks) {
    const auto rawHit = link.getFrom();
    const auto simHit = link.getTo();
    if (!rawHit.isAvailable() || !simHit.isAvailable() || simHit.isProducedBySecondary()) {
      continue;
    }
    const auto particle = simHit.getParticle();
    if (!particle.isAvailable()) {
      continue;
    }
    particlesOfRawHit[rawHit.getObjectID()].emplace_back(simHit.getTime(), particle.getObjectID());
  }

  // The hits of each particle, as (time, hit index)
  const double toActs = Acts::UnitConstants::mm / edm4eic::unit::mm;
  std::vector<Acts::Vector3> positions;
  positions.reserve(hits.size());
  std::map<podio::ObjectID, std::vector<std::pair<float, std::uint32_t>>, ObjectIDLess>
      hitsOfParticle;
  for (std::uint32_t i = 0; i < hits.size(); ++i) {
    const auto pos = hits[i].getPosition();
    positions.emplace_back(toActs * Acts::Vector3(pos.x, pos.y, pos.z));
    const auto rawHit = hits[i].getRawHit();
    if (!rawHit.isAvailable()) {
      continue;
    }
    const auto it = particlesOfRawHit.find(rawHit.getObjectID());
    if (it == particlesOfRawHit.end()) {
      continue;
    }
    for (const auto& [time, particle] : it->second) {
      hitsOfParticle[particle].emplace_back(time, i);
    }
  }

  auto fieldCache = fieldProvider.makeCache(fieldContext);
  for (auto& [particle, particleHits] : hitsOfParticle) {
    std::sort(particleHits.begin(), particleHits.end());
    std::vector<std::uint32_t> kept;
    for (const auto& [time, i] : particleHits) {
      if (kept.size() == m_cfg.maxHitsPerParticle) {
        break;
      }
      if (kept.empty() || (positions[i] - positions[kept.back()]).norm() >= m_cfg.minHitDistance) {
        kept.push_back(i);
      }
    }

    for (std::size_t a = 0; a < kept.size(); ++a) {
      for (std::size_t b = a + 1; b < kept.size(); ++b) {
        const auto field = fieldProvider.getField(positions[kept[b]], fieldCache);
        // No field outside the field map
        if (!field.ok() || !(field.value().norm() > 0.)) {
          continue;
        }
        for (std::size_t c = b + 1; c < kept.size(); ++c) {
          addHelixSeed(m_cfg.prior, {kept[a], kept[b], kept[c]}, 0.F, hits, positions,
                       field.value(), seeds, parameters);
        }
      }
    }
  }
}

} // namespace eicrecon
