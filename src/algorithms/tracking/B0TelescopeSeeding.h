// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "B0TelescopeGeometry.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include <algorithms/algorithm.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <memory>
#include <mutex>
#include <edm4hep/EventHeaderCollection.h>
#include <string>
#include <string_view>

namespace eicrecon {
struct B0TelescopeSeedingConfig {
  double stationZGap{50.0};
  b0::SearchOptions search;
  double curvatureSignificance{3.0};
  double proposalMaxMomentum{10000.0};
  int forceCharge{0};
  bool bothCharges{false};
  // Search-covariance floors, not detector resolutions or final-fit priors.
  double seedPositionSigma{0.2}; // mm
  double seedAngularSigma{0.001}; // radians
  double seedQOverPSigma{0.02}; // GeV^-1
  std::string stationMapFile;
  std::string diagnosticsFile;
};
using B0TelescopeSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4hep::EventHeaderCollection, edm4eic::Measurement2DCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;
class B0TelescopeSeeding : public B0TelescopeSeedingAlgorithm,
                          public WithPodConfig<B0TelescopeSeedingConfig> {
public:
  explicit B0TelescopeSeeding(std::string_view name)
      : B0TelescopeSeedingAlgorithm{name, {"header", "measurements"}, {"seeds", "parameters"},
                                    "Local, station-aware B0 telescope seeding"} {}
  void init() final;
  void process(const Input&, const Output&) const final;
private:
  std::shared_ptr<const ActsGeometryProvider> m_provider;
  std::unique_ptr<b0::TelescopeGeometry> m_geometry;
  mutable std::mutex m_fileMutex;
};
namespace b0 {
/// Bind a local ion-frame state to the first physical sensor plane. Used by
/// both realistic and truth seeds. The fit state is NOT an origin perigee.
void writeTelescopeSeed(const TelescopeFit& fit, const Acts::Surface& surface,
                        const Acts::GeometryContext& context, const Eigen::Matrix3d& labToIon,
                        double time, double timeVariance,
                        edm4eic::MutableTrackParameters parameters);
}
} // namespace eicrecon
