// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "algorithms/tracking/B0TripletSeeding.h"
#include "algorithms/tracking/B0TripletSeedingConfig.h"
#include "extensions/jana/JOmniFactory.h"
#include "services/geometry/acts/ACTSGeo_service.h"

namespace eicrecon {

class B0TripletSeeding_factory
    : public JOmniFactory<B0TripletSeeding_factory, B0TripletSeedingConfig> {

private:
  using AlgoT = eicrecon::B0TripletSeeding;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::TrackerHit> m_hits_input{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds_output{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};

  Service<ACTSGeo_service> m_ACTSGeoSvc{this};

  ParameterRef<double> m_stationGap{this, "stationGap", config().stationGap,
                                    "min distance in z between stations"};
  ParameterRef<double> m_maxResidual{
      this, "maxResidual", config().maxResidual,
      "max distance of the middle hit from the chord along the magnetic field"};
  ParameterRef<double> m_minMomentum{this, "minMomentum", config().minMomentum,
                                     "min momentum of the seed"};
  ParameterRef<double> m_anchorDistance{
      this, "anchorDistance", config().anchorDistance,
      "distance upstream of the first hit at which seeds are expressed"};
  ParameterRef<std::size_t> m_maxSeeds{this, "maxSeeds", config().maxSeeds,
                                       "max number of seeds per event"};
  ParameterRef<double> m_positionError{this, "position_Error", config().positionError,
                                       "error on the position at the first hit"};
  ParameterRef<double> m_angleError{this, "angle_Error", config().angleError,
                                    "error on the direction at the first hit"};
  ParameterRef<double> m_qOverPRelativeError{this, "qOverP_RelativeError",
                                             config().qOverPRelativeError, "relative error on q/p"};
  ParameterRef<double> m_timeError{this, "time_Error", config().timeError, "error on time"};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void Process(int32_t /* run_number */, uint64_t /* event_number */) {
    m_algo->process({m_hits_input()}, {m_seeds_output().get(), m_trackparams_output().get()});
  }
};

} // namespace eicrecon
