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

  // Not read here: the algorithm takes the ACTS geometry and field from algorithms::ActsSvc. The
  // member declares that dependency, as in the other ACTS factories, and makes JANA acquire the
  // service before Configure() constructs the algorithm. Do not remove it as unused.
  Service<ACTSGeo_service> m_ACTSGeoSvc{this};

  ParameterRef<double> m_stationGap{
      this, "stationGap", config().stationGap,
      "z gap between consecutive hits that starts a new station, in Acts units (mm)"};
  ParameterRef<double> m_maxResidual{
      this, "maxResidual", config().maxResidual,
      "max distance of the middle hit from the chord along the magnetic field, in Acts units (mm)"};
  ParameterRef<double> m_minMomentum{this, "minMomentum", config().minMomentum,
                                     "min momentum of the seed, in Acts units (GeV)"};
  ParameterRef<double> m_anchorDistance{
      this, "anchorDistance", config().anchorDistance,
      "distance upstream of the first hit at which seeds are expressed, in Acts units (mm)"};
  ParameterRef<std::size_t> m_maxSeeds{this, "maxSeeds", config().maxSeeds,
                                       "max number of seeds per event"};
  ParameterRef<double> m_timeWindowSigmas{
      this, "timeWindowSigmas", config().timeWindowSigmas,
      "hit-time window width in units of the combined hit time error, 0 turns it off"};
  ParameterRef<double> m_timeWindowTolerance{this, "timeWindowTolerance",
                                             config().timeWindowTolerance,
                                             "extra width of the hit-time window, in ns"};
  ParameterRef<double> m_timeErrorFloor{this, "timeErrorFloor", config().timeErrorFloor,
                                        "min hit time error in the time window, in ns"};
  ParameterRef<double> m_slowParticleMass{
      this, "slowParticleMass", config().slowParticleMass,
      "mass per unit charge of the slowest particles to seed, in Acts units (GeV)"};
  ParameterRef<std::size_t> m_fieldSamples{
      this, "fieldSamples", config().fieldSamples,
      "field samples along a triplet for its effective field, 0 for the field at the middle hit"};
  ParameterRef<double> m_fieldTermAllowance{
      this, "fieldTermAllowance", config().fieldTermAllowance,
      "extra raw residual admitted before the field term is subtracted, in Acts units (mm)"};
  ParameterRef<double> m_positionError{this, "positionError", config().positionError,
                                       "error on the seed position, in Acts units (mm)"};
  ParameterRef<double> m_angleError{this, "angleError", config().angleError,
                                    "error on the seed direction, in Acts units (rad)"};
  ParameterRef<double> m_qOverPRelativeError{this, "qOverPRelativeError",
                                             config().qOverPRelativeError, "relative error on q/p"};
  ParameterRef<double> m_timeError{this, "timeError", config().timeError,
                                   "error on the seed time, in ns"};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void Process(int32_t /* run_number */, uint64_t /* event_number */) override {
    m_algo->process({m_hits_input()}, {m_seeds_output().get(), m_trackparams_output().get()});
  }
};

} // namespace eicrecon
