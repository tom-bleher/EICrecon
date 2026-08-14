// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <JANA/JEvent.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <memory>

#include "algorithms/tracking/B0TrackerStubSeeder.h"
#include "algorithms/tracking/B0TrackerStubSeederConfig.h"
#include "extensions/jana/JOmniFactory.h"

namespace eicrecon {

class B0TrackerStubSeeder_factory
    : public JOmniFactory<B0TrackerStubSeeder_factory, B0TrackerStubSeederConfig> {

private:
  using AlgoT = eicrecon::B0TrackerStubSeeder;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::Measurement2D> m_measurements_input{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds_output{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};

  ParameterRef<std::string> m_readout{this, "readout", config().readout,
                                      "DD4hep readout used to decode B0 layers"};
  ParameterRef<float> m_crossingAngle{this, "crossingAngle", config().crossingAngle,
                                      "ion-beam rotation about y [rad]"};
  ParameterRef<float> m_bFieldY{this, "bFieldY", config().bFieldY, "B0pf dipole field along y [T]"};
  ParameterRef<float> m_zFieldEntrance{this, "zFieldEntrance", config().zFieldEntrance,
                                       "ion-frame z of B0pf field entrance [mm]"};
  ParameterRef<float> m_stationZGap{this, "stationZGap", config().stationZGap,
                                    "ion-frame z gap [mm] that starts a new B0 station"};
  ParameterRef<unsigned int> m_minStations{this, "minStations", config().minStations,
                                           "minimum distinct stations per seed"};
  ParameterRef<unsigned int> m_maxCombinations{this, "maxCombinations", config().maxCombinations,
                                               "cap on hit combinations per event"};
  ParameterRef<unsigned int> m_maxSeeds{this, "maxSeeds", config().maxSeeds,
                                        "maximum seeds per event"};
  ParameterRef<unsigned int> m_maxSharedHits{this, "maxSharedHits", config().maxSharedHits,
                                             "max shared hits between accepted seeds"};
  ParameterRef<int> m_charge{this, "charge", config().charge, "assumed charge"};
  ParameterRef<bool> m_constrainToBeamline{this, "constrainToBeamline",
                                           config().constrainToBeamline,
                                           "use origin->field-entrance chord for the direction"};
  ParameterRef<float> m_momentumPrior{this, "momentumPrior", config().momentumPrior,
                                      "fallback momentum [GeV]"};
  ParameterRef<float> m_pMin{this, "pMin", config().pMin, "min accepted fitted momentum [GeV]"};
  ParameterRef<float> m_pMax{this, "pMax", config().pMax, "max accepted fitted momentum [GeV]"};
  ParameterRef<float> m_locaError{this, "locaError", config().locaError, "seed loc0 error"};
  ParameterRef<float> m_locbError{this, "locbError", config().locbError, "seed loc1 error"};
  ParameterRef<float> m_phiError{this, "phiError", config().phiError, "seed phi error"};
  ParameterRef<float> m_thetaError{this, "thetaError", config().thetaError, "seed theta error"};
  ParameterRef<float> m_qOverPError{this, "qOverPError", config().qOverPError, "seed q/p error"};
  ParameterRef<float> m_timeError{this, "timeError", config().timeError, "seed time error"};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void ChangeRun(int32_t /* run_number */) {}

  void Process(int32_t /* run_number */, uint64_t /* event_number */) {
    m_algo->process({m_measurements_input()},
                    {m_seeds_output().get(), m_trackparams_output().get()});
  }
};

} // namespace eicrecon
