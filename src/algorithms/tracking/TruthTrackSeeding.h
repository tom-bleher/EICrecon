// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <algorithms/algorithm.h>
#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <memory>
#include <string>
#include <string_view>

#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"
#include "algorithms/tracking/TruthTrackSeedingConfig.h"

namespace eicrecon {

using TruthTrackSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection, edm4eic::MCRecoTrackerHitLinkCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Track seeds from truth pattern recognition: each particle is seeded from triplets of its own
/// hits, with the helix estimate of the hit seeders (addHelixSeed).
///
/// A hit belongs to every particle that left one of its sim hits, except through a secondary
/// (sim hits flagged produced-by-secondary). A particle's hits are taken in time order, skipping
/// hits closer than minHitDistance to the last one kept, up to maxHitsPerParticle, and every
/// triplet of them becomes a seed, in the field at its middle hit. Particles with fewer than three
/// such hits get no seed. As in ActsExamples' TruthSeedingAlgorithm, this is the seeding of a
/// perfect pattern recognition, so it differs from the hit seeds only in which triplets it takes.
class TruthTrackSeeding : public TruthTrackSeedingAlgorithm,
                          public WithPodConfig<TruthTrackSeedingConfig> {
public:
  TruthTrackSeeding(std::string_view name)
      : TruthTrackSeedingAlgorithm{name,
                                   {"inputTrackerHits", "inputRawHitLinks"},
                                   {"outputTrackSeeds", "outputTrackParameters"},
                                   "create track seeds from triplets of each particle's own hits"} {
  }

  void init() final {}
  void process(const Input& input, const Output& output) const final;

  /// The truth seeds of the hits in a field
  void findSeeds(const edm4eic::TrackerHitCollection& hits,
                 const edm4eic::MCRecoTrackerHitLinkCollection& rawHitLinks,
                 const Acts::MagneticFieldProvider& fieldProvider,
                 const Acts::MagneticFieldContext& fieldContext,
                 edm4eic::TrackSeedCollection& seeds,
                 edm4eic::TrackParametersCollection& parameters) const;

private:
  const std::shared_ptr<const ActsGeometryProvider> m_geoSvc{
      algorithms::ActsSvc::instance().acts_geometry_provider()};
};

} // namespace eicrecon
