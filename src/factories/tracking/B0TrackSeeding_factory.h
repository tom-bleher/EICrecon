// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "algorithms/tracking/B0TrackSeeding.h"
#include "algorithms/tracking/B0TrackSeedingConfig.h"
#include "extensions/jana/JOmniFactory.h"
#include "services/geometry/acts/ACTSGeo_service.h"

namespace eicrecon {

class B0TrackSeeding_factory : public JOmniFactory<B0TrackSeeding_factory, B0TrackSeedingConfig> {

private:
  using AlgoT = eicrecon::B0TrackSeeding;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::TrackerHit> m_hits_input{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds_output{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};

  // The algorithm takes the ACTS geometry and field from algorithms::ActsSvc; this declares the
  // dependency, so that JANA acquires the service before Configure()
  Service<ACTSGeo_service> m_ACTSGeoSvc{this};

  ParameterRef<double> m_rMin{this, "rMin", config().rMin,
                              "min distance of a hit along the ion beam, mm"};
  ParameterRef<double> m_rMax{this, "rMax", config().rMax,
                              "max distance of a hit along the ion beam, mm"};
  ParameterRef<double> m_zMax{this, "zMax", config().zMax,
                              "max distance of a hit from the beam along the field, mm"};
  ParameterRef<double> m_rMinMiddle{this, "rMinMiddle", config().rMinMiddle,
                                    "min distance of a middle hit along the ion beam, mm"};
  ParameterRef<double> m_rMaxMiddle{this, "rMaxMiddle", config().rMaxMiddle,
                                    "max distance of a middle hit along the ion beam, mm"};
  ParameterRef<double> m_deltaRMin{this, "deltaRMin", config().deltaRMin,
                                   "min distance along the beam of a doublet, mm"};
  ParameterRef<double> m_deltaRMax{this, "deltaRMax", config().deltaRMax,
                                   "max distance along the beam of a doublet, mm"};
  ParameterRef<double> m_cotThetaMax{this, "cotThetaMax", config().cotThetaMax,
                                     "max slope of a doublet along the field"};
  ParameterRef<double> m_minPt{this, "minPt", config().minPt,
                               "min momentum across the field of the Acts finders, GeV"};
  ParameterRef<double> m_sigmaScattering{this, "sigmaScattering", config().sigmaScattering,
                                         "tolerance of the Acts scattering cut, in sigmas"};
  ParameterRef<double> m_radLengthPerSeed{this, "radLengthPerSeed", config().radLengthPerSeed,
                                          "material per station of the scattering cut, in X0"};
  ParameterRef<double> m_hitErrorScale{
      this, "hitErrorScale", config().hitErrorScale,
      "factor on the hit errors given to the Acts finders, the minimum sigmas of their slope cut"};
  ParameterRef<double> m_fieldBinWidth{this, "fieldBinWidth", config().fieldBinWidth,
                                       "|B| bin width of the Acts finders, tesla"};
  ParameterRef<std::size_t> m_maxSeedsPerMiddle{
      this, "maxSeedsPerMiddle", config().maxSeedsPerMiddle, "max seeds per middle hit"};
  ParameterRef<double> m_timeWindowSigmas{this, "timeWindowSigmas", config().timeWindowSigmas,
                                          "hit-time window in units of the hit time errors, 0 off"};
  ParameterRef<double> m_timeErrorFloor{this, "timeErrorFloor", config().timeErrorFloor,
                                        "min hit time error in the time window, ns"};
  ParameterRef<double> m_slowParticleMass{
      this, "slowParticleMass", config().slowParticleMass,
      "mass per unit charge of the slowest seeded particle, GeV"};
  ParameterRef<double> m_minMomentum{this, "minMomentum", config().minMomentum,
                                     "min seed momentum, GeV"};
  ParameterRef<double> m_anchorDistance{
      this, "anchorDistance", config().anchorDistance,
      "distance of the seed perigee upstream of its first hit, mm"};
  ParameterRef<std::size_t> m_maxSeeds{this, "maxSeeds", config().maxSeeds, "max seeds per event"};
  ParameterRef<double> m_positionError{this, "positionError", config().positionError,
                                       "seed position error, mm"};
  ParameterRef<double> m_angleError{this, "angleError", config().angleError,
                                    "seed direction error, rad"};
  ParameterRef<double> m_qOverPRelativeError{
      this, "qOverPRelativeError", config().qOverPRelativeError, "relative seed q/p error"};
  ParameterRef<double> m_timeError{this, "timeError", config().timeError, "seed time error, ns"};

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
