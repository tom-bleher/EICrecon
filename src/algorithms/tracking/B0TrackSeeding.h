// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Seeding2/TripletSeeder.hpp>
#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"
#include "algorithms/tracking/B0TrackSeedingConfig.h"

namespace eicrecon {

using B0TrackSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// B0 track seeds from the Acts triplet seeder.
///
/// The Acts finders assume a uniform field along z. The field of the B0 dipole changes direction
/// off its midplane, so each middle hit gets its own frame, a rotation that takes the field at that
/// hit to z. The finders work on differences from the middle hit, so that is all they need.
class B0TrackSeeding : public B0TrackSeedingAlgorithm, public WithPodConfig<B0TrackSeedingConfig> {
public:
  B0TrackSeeding(std::string_view name)
      : B0TrackSeedingAlgorithm{name,
                                {"inputTrackerHits"},
                                {"outputTrackSeeds", "outputTrackParameters"},
                                "create B0 track seeds with the Acts triplet seeder"} {}

  void process(const Input& input, const Output& output) const final;

  /// Adds the seed of three hits on a helix in the field at the middle one, unless they have no
  /// finite non-zero q/p. Positions are in Acts units, the time in EDM4eic units.
  static std::optional<edm4eic::MutableTrackSeed>
  addSeed(const B0TrackSeedingConfig& cfg, const Acts::Vector3& first, const Acts::Vector3& middle,
          const Acts::Vector3& last, const Acts::Vector3& field, float time,
          edm4eic::TrackSeedCollection& seeds, edm4eic::TrackParametersCollection& parameters);

private:
  const std::shared_ptr<const ActsGeometryProvider> m_geoSvc{
      algorithms::ActsSvc::instance().acts_geometry_provider()};
  Acts::TripletSeeder m_seeder;
};

} // namespace eicrecon
