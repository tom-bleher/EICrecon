// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#pragma once

#include "algorithms/tracking/B0TelescopeSeeding.h"
#if EICRECON_HAS_B0_TELESCOPE
#include "extensions/jana/JOmniFactory.h"
#include "services/geometry/acts/ACTSGeo_service.h"
#include <memory>

namespace eicrecon {
class B0TelescopeSeeding_factory
    : public JOmniFactory<B0TelescopeSeeding_factory, B0TelescopeSeedingConfig> {
private:
  std::unique_ptr<B0TelescopeSeeding> m_algo;
  PodioInput<edm4eic::TrackerHit> m_hits{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds{this};
  PodioOutput<edm4eic::TrackParameters> m_parameters{this};
  Service<ACTSGeo_service> m_geometry{this};

  ParameterRef<unsigned> m_layers{this, "layersPerStation", config().layersPerStation,
                                 "Readout layers per physical station: 2 realistic, 1 ideal"};
  ParameterRef<double> m_dz{this, "minDeltaZ", config().minDeltaZ, "Minimum station separation [mm]"};
  ParameterRef<double> m_slope{this, "maxSlope", config().maxSlope, "Maximum absolute telescope x/z or y/z slope"};
  ParameterRef<double> m_rx{this, "maxResidualX", config().maxResidualX, "Bending-plane sagitta window [mm]"};
  ParameterRef<double> m_ry{this, "maxResidualY", config().maxResidualY, "Non-bending sagitta window [mm]"};
  ParameterRef<std::size_t> m_combinations{this, "maxCombinations", config().maxCombinations, "Reject whole event beyond this triplet budget"};
  ParameterRef<std::size_t> m_perMiddle{this, "maxSeedsPerMiddle", config().maxSeedsPerMiddle, "Maximum seeds per middle point per station triple"};
  ParameterRef<std::size_t> m_maxSeeds{this, "maxSeeds", config().maxSeeds, "Maximum ranked seeds per event"};
  ParameterRef<double> m_position{this, "positionSigma", config().positionSigma, "Initial position uncertainty at first-hit plane [mm]"};
  ParameterRef<double> m_angle{this, "angleSigma", config().angleSigma, "Initial angular uncertainty [rad]"};
  ParameterRef<double> m_relativeP{this, "relativeMomentumSigma", config().relativeMomentumSigma, "Initial relative momentum uncertainty"};
  ParameterRef<double> m_qop{this, "qOverPSigmaFloor", config().qOverPSigmaFloor, "Minimum initial q/p uncertainty [1/GeV]"};
  ParameterRef<double> m_time{this, "timeSigma", config().timeSigma, "Initial time uncertainty [ns]"};
  ParameterRef<double> m_path{this, "propagationPathLimit", config().propagationPathLimit, "Maximum backward seed transport path [mm]"};

public:
  void Configure() {
    m_algo = std::make_unique<B0TelescopeSeeding>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }
  void Process(int32_t, uint64_t) {
    m_algo->process({m_hits()}, {m_seeds().get(), m_parameters().get()});
  }
};
} // namespace eicrecon
#endif
