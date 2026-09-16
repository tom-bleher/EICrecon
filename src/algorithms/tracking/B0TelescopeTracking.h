// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "B0TelescopeGeometry.h"
#include "CKFTracking.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include <edm4hep/EventHeaderCollection.h>
#include <array>
#include <memory>
#include <mutex>
#include <string>

namespace eicrecon {
struct B0TelescopeTrackingConfig {
  bool doFinalWeakPriorRefit{false}; // explicit; never inferred from station count
  bool fitSeedMeasurementsDirectly{false};
  int particleHypothesisPdg{2212};
  bool useSeedParticleHypothesis{false}; // truth-reference studies only
  std::size_t minStations{3};
  double stationZGap{50.0};
  double chi2Cut{50.0};
  std::size_t maxBranchesPerSurface{5};
  std::array<double, 6> weakPriorVariances{1e6, 1e6, 1.0, 1.0, 1.0, 1e6};
  double weakPriorScale{1.0};
  double initializationStepMm{1.0};
  // Prompt mode is a POST-FIT selection, NOT a hit-derived beam-spot prior.
  // Displaced mode (false) never requires extrapolation to the interaction point.
  bool promptSelection{false};
  double promptMaxD0Mm{5.0};
  double promptMaxZ0Mm{200.0};
  std::string diagnosticsFile; // optional per-chain JSONL, including innovations
};
using B0TelescopeTrackingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4hep::EventHeaderCollection, edm4eic::TrackSeedCollection,
                      edm4eic::Measurement2DCollection>,
    algorithms::Output<Acts::ConstVectorMultiTrajectory*, Acts::ConstVectorTrackContainer*>>;
class B0TelescopeTracking : public B0TelescopeTrackingAlgorithm,
                           public WithPodConfig<B0TelescopeTrackingConfig> {
public:
  explicit B0TelescopeTracking(std::string_view name)
      : B0TelescopeTrackingAlgorithm{name, {"header", "seeds", "measurements"},
                                     {"states", "tracks"}, "B0 finding and independent fitting"} {}
  void init() final;
  void process(const Input&, const Output&) const final;
private:
  std::shared_ptr<const ActsGeometryProvider> m_provider;
  std::unique_ptr<b0::TelescopeGeometry> m_geometry;
  std::shared_ptr<CKFTracking::CKFTrackingFunction> m_finder;
  std::unique_ptr<const Acts::Logger> m_logger;
  mutable std::mutex m_diagnosticsMutex;
};
} // namespace eicrecon
