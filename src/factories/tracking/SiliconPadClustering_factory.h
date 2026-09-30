// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include "algorithms/tracking/SiliconPadClustering.h"
#include "extensions/jana/JOmniFactory.h"
#include "services/geometry/dd4hep/DD4hep_service.h"

namespace eicrecon {

class SiliconPadClustering_factory
    : public JOmniFactory<SiliconPadClustering_factory, SiliconPadClusteringConfig> {

public:
  using AlgoT = eicrecon::SiliconPadClustering;

private:
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::TrackerHit> m_pad_hits_input{this};
  PodioOutput<edm4eic::TrackerHit> m_cluster_hits_output{this};

  ParameterRef<std::string> m_readout{this, "readout", config().readout};
  ParameterRef<double> m_deltaT{this, "deltaT", config().deltaT};
  ParameterRef<double> m_electrode_x{this, "electrodeX", config().electrode_x};
  ParameterRef<double> m_electrode_y{this, "electrodeY", config().electrode_y};
  ParameterRef<double> m_single_pad_resolution{this, "singlePadResolution",
                                               config().single_pad_resolution};
  ParameterRef<double> m_shared_pad_resolution{this, "sharedPadResolution",
                                               config().shared_pad_resolution};

  Service<DD4hep_service> m_geoSvc{this};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void Process(int32_t /* run_number */, uint64_t /* event_number */) {
    m_algo->process({m_pad_hits_input()}, {m_cluster_hits_output().get()});
  }
};

} // namespace eicrecon
