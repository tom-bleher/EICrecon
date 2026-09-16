// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "algorithms/digi/B0ACLGADCalibration.h"
#include "algorithms/interfaces/UniqueIDGenSvc.h"
#include "algorithms/tracking/ActsGeometryProvider.h"
#include "extensions/jana/JOmniFactory.h"
#include <DD4hep/VolumeManager.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <Eigen/Core>
#include <map>
#include <memory>
#include <string>

namespace eicrecon {
struct B0ACLGADConfig { std::string calibrationFile; bool allowExperimental = false; };
class B0ACLGAD_factory : public JOmniFactory<B0ACLGAD_factory, B0ACLGADConfig> {
  struct Binding {
    const dd4hep::VolumeManagerContext* context;
    const Acts::Surface* surface;
    double halfX, halfY, halfZ;
    Eigen::Matrix2d jacobian;
    Acts::Vector2 origin;
  };
  PodioInput<edm4hep::EventHeader> m_headers{this};
  PodioInput<edm4hep::SimTrackerHit> m_simHits{this};
  PodioOutput<edm4eic::RawTrackerHit> m_raw{this};
  PodioOutput<edm4eic::MCRecoTrackerHitLink> m_links{this};
  PodioOutput<edm4eic::MCRecoTrackerHitAssociation> m_associations{this};
  PodioOutput<edm4eic::TrackerHit> m_channelHits{this};
  PodioOutput<edm4eic::Measurement2D> m_measurements{this};
  ParameterRef<std::string> m_file{this, "calibrationFile", config().calibrationFile,
                                    "Required physical-pitch B0 response calibration (.cfg)"};
  ParameterRef<bool> m_experimental{this, "allowExperimental", config().allowExperimental,
                                    "Explicitly permit uncalibrated response studies (not production)"};
  b0::aclgad::Calibration m_calibration;
  std::shared_ptr<const ActsGeometryProvider> m_provider;
  std::map<std::uint64_t, Binding> m_bindings;
  std::vector<b0::aclgad::Sensor> m_sensors;
  const algorithms::UniqueIDGenSvc& m_uid = algorithms::UniqueIDGenSvc::instance();
public:
  void Configure();
  void ChangeRun(int32_t) {}
  void Process(int32_t, uint64_t);
};
} // namespace eicrecon
