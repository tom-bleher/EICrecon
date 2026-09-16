// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "B0TelescopeSeeding.h"
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <mutex>

namespace eicrecon {
struct B0TelescopeTruthSeedingConfig {
  std::size_t minStations{3};
  double stationZGap{50.0};
  double minMomentumGeV{0.1};
  double relativeMomentumSmear{0.1};
  double positionSigmaMm{1.0};
  double slopeSigma{0.01};
  double qOverPSigma{0.1};
  double timeSigmaNs{1.0};
  double measurementPurityMin{0.5}; // strict inequality: ties are not truth matches
  double maxProjectionMm{0.1}; // bound the thin-sensor tangent approximation
  std::string referenceFile;
};
using B0TelescopeTruthSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4hep::EventHeaderCollection, edm4hep::SimTrackerHitCollection,
                      edm4eic::Measurement2DCollection, edm4eic::MCRecoTrackerHitAssociationCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;
class B0TelescopeTruthSeeding : public B0TelescopeTruthSeedingAlgorithm,
                               public WithPodConfig<B0TelescopeTruthSeedingConfig> {
public:
  explicit B0TelescopeTruthSeeding(std::string_view name)
      : B0TelescopeTruthSeedingAlgorithm{name, {"header", "simHits", "measurements", "rawAssociations"},
                                         {"seeds", "parameters"}, "Truth reference from actual B0 crossings"} {}
  void init() final;
  void process(const Input&, const Output&) const final;
private:
  std::shared_ptr<const ActsGeometryProvider> m_provider;
  std::unique_ptr<b0::TelescopeGeometry> m_geometry;
  mutable std::mutex m_fileMutex;
};
} // namespace eicrecon
