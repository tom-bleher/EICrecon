// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Seeding2/BroadTripletSeedFilter.hpp>
#include <Acts/Seeding2/DoubletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeeder.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <Eigen/Core>
#include <cstddef>
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
/// The Acts finders assume a uniform field along z, layers in radius and tracks from a beam line.
/// The B0 tracker is four planar stations in a dipole whose field direction varies off the
/// midplane, and its tracks include decay products from upstream. So each middle hit gets its own
/// frame, z' along the local field and x' along the ion beam, in which the stations are layers at
/// r' ~ distance along the beam; the beam-line cuts are opened.
class B0TrackSeeding : public B0TrackSeedingAlgorithm, public WithPodConfig<B0TrackSeedingConfig> {
public:
  B0TrackSeeding(std::string_view name)
      : B0TrackSeedingAlgorithm{name,
                                {"inputTrackerHits"},
                                {"outputTrackSeeds", "outputTrackParameters"},
                                "create B0 track seeds with the Acts triplet seeder"} {}

  void init() final;
  void process(const Input& input, const Output& output) const final;

  /// A hit in Acts units, with the field at its position
  struct Hit {
    Acts::Vector3 position;
    Acts::Vector3 field;
    double time;
    double variance;
    std::size_t index;
  };

  /// Seed parameters on the perigee at their anchor, and their diagonal prior
  struct Seed {
    Acts::Vector3 anchor;
    Acts::BoundVector parameters;
    Acts::BoundVector errors;
  };

  /// Global to seeding frame; rows x' (beam, across the field), y', z' (field)
  static Acts::RotationMatrix3 frameRotation(const Acts::Vector3& field, const Acts::Vector3& beam);

  /// Seed parameters of a hit triplet, none below minMomentum
  static std::optional<Seed> makeSeed(const B0TrackSeedingConfig& cfg, const Hit& a, const Hit& b,
                                      const Hit& c);

private:
  const algorithms::ActsSvc& m_actsSvc{algorithms::ActsSvc::instance()};
  const std::shared_ptr<const ActsGeometryProvider> m_geoSvc{m_actsSvc.acts_geometry_provider()};
  Acts::Vector3 m_beam{0., 0., 1.};

  // The Acts seeding setup, built in init(); the finders are created for the field at each middle
  // hit
  Acts::DoubletSeedFinder::Config m_bottomConfig;
  Acts::DoubletSeedFinder::Config m_topConfig;
  Acts::TripletSeedFinder::Config m_tripletConfig;
  Acts::BroadTripletSeedFilter::Config m_filterConfig;
  std::unique_ptr<const Acts::Logger> m_filterLogger;
  std::unique_ptr<const Acts::TripletSeeder> m_seeder;
};

} // namespace eicrecon
