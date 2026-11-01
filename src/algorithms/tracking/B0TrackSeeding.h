// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <cstddef>
#include <memory>
#include <optional>
#include <string_view>

#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"
#include "algorithms/tracking/B0TrackSeedingConfig.h"

namespace eicrecon {

using B0TrackSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Track seeds for the B0 tracker with the Acts triplet seeder (DoubletSeedFinder,
/// TripletSeedFinder, BroadTripletSeedFilter, TripletSeeder) and
/// Acts::estimateTrackParamsFromSeed.
///
/// The Acts finders assume a uniform field along their z axis, layers ordered in radius and
/// tracks from a beam line. The B0 tracker is four planar stations inside a dipole with a
/// quadrupole gradient, whose field turns off the midplane, and its tracks include decay products
/// from anywhere upstream. So for each middle hit the hits are rotated into a frame whose z' axis
/// is along the field at that hit and whose x' axis is along the ion beam: the stations become
/// layers at r' = distance along the beam, and the field is locally a field along z' of the
/// magnitude there. The cuts and seed weights that assume a beam line are opened up, and a hit-time
/// window rejects hits from other bunch crossings. Each seed is expressed on a perigee just
/// upstream of its first hit.
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
    double timeError;
    double variance;
    std::size_t index;
  };

  /// Seed parameters on the perigee at their anchor, and their diagonal prior
  struct Seed {
    Acts::Vector3 anchor;
    Acts::BoundVector parameters;
    Acts::BoundVector errors;
  };

  /// Rotation from global coordinates to the seeding frame: its rows are x' (the beam direction
  /// across the field), y' and z' (the field direction)
  static Acts::RotationMatrix3 frameRotation(const Acts::Vector3& field, const Acts::Vector3& beam);

  /// Whether the later hit can follow the first one in time (see B0TrackSeedingConfig)
  static bool inTimeWindow(const B0TrackSeedingConfig& cfg, const Hit& first, const Hit& later);

  /// Seed parameters of a hit triplet in the field at its middle hit, if it passes minMomentum
  static std::optional<Seed> makeSeed(const B0TrackSeedingConfig& cfg,
                                      const Acts::GeometryContext& gctx, const Hit& a, const Hit& b,
                                      const Hit& c);

  struct Finders;

private:
  const algorithms::ActsSvc& m_actsSvc{algorithms::ActsSvc::instance()};
  const std::shared_ptr<const ActsGeometryProvider> m_geoSvc{m_actsSvc.acts_geometry_provider()};
  Acts::Vector3 m_beam{0., 0., 1.};
  std::shared_ptr<const Finders> m_finders;
};

} // namespace eicrecon
