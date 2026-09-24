// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <memory>
#include <string>
#include <string_view>

#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"
#include "algorithms/tracking/B0TripletSeedingConfig.h"

namespace eicrecon {

using B0TripletSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Track seeds for the B0 tracker, a telescope of planar stations inside the
/// B0 magnet. Every triplet of hits on three stations that is compatible with
/// a helix in the local field becomes a seed, without assuming a vertex or a
/// charge, so tracks from displaced decays are seeded as well.
class B0TripletSeeding : public B0TripletSeedingAlgorithm,
                         public WithPodConfig<B0TripletSeedingConfig> {
public:
  B0TripletSeeding(std::string_view name)
      : B0TripletSeedingAlgorithm{name,
                                  {"inputTrackerHits"},
                                  {"outputTrackSeeds", "outputTrackParameters"},
                                  "create track seeds from B0 tracker hit triplets"} {}

  void init() final {};
  void process(const Input& input, const Output& output) const final;

private:
  const algorithms::ActsSvc& m_actsSvc{algorithms::ActsSvc::instance()};
  std::shared_ptr<const ActsGeometryProvider> m_geoSvc{m_actsSvc.acts_geometry_provider()};
};

} // namespace eicrecon
