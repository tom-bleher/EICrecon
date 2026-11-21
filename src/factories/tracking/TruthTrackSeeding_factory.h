// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "algorithms/tracking/TruthTrackSeeding.h"
#include "algorithms/tracking/TruthTrackSeedingConfig.h"
#include "extensions/jana/JOmniFactory.h"
#include "services/geometry/acts/ACTSGeo_service.h"

namespace eicrecon {

class TruthTrackSeeding_factory
    : public JOmniFactory<TruthTrackSeeding_factory, TruthTrackSeedingConfig> {

private:
  using AlgoT = eicrecon::TruthTrackSeeding;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::TrackerHit> m_hits_input{this};
  PodioInput<edm4eic::MCRecoTrackerHitLink> m_raw_hit_links_input{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds_output{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};

  ParameterRef<double> m_minHitDistance{
      this, "minHitDistance", config().minHitDistance,
      "min distance between consecutive hits of a particle that seed it, mm"};
  ParameterRef<std::size_t> m_maxHitsPerParticle{
      this, "maxHitsPerParticle", config().maxHitsPerParticle,
      "max hits per particle, in time order; every triplet of them is a seed"};

  Service<ACTSGeo_service> m_ACTSGeoSvc{this};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void Process(int32_t /* run_number */, uint64_t /* event_number */) {
    m_algo->process({m_hits_input(), m_raw_hit_links_input()},
                    {m_seeds_output().get(), m_trackparams_output().get()});
  }
};

} // namespace eicrecon
