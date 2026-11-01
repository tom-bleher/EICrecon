// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/ActsVersion.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"
#include "algorithms/tracking/B0Seeding2Config.h"

// The Seeding2 triplet seeder is in Acts 45.3 to 46.x; Acts 47 renamed it to Seeding and its
// containers (SpacePointContainer2 to SpacePointContainer, ...), which this class does not use yet
#define B0SEEDING2_HAS_SEEDING2                                                                    \
  ((Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 3)) &&           \
   (Acts_VERSION_MAJOR < 47))

namespace eicrecon {

using B0Seeding2Algorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Track seeds for the B0 tracker from the Acts Seeding2 triplet seeder (DoubletSeedFinder,
/// TripletSeedFinder, BroadTripletSeedFilter, TripletSeeder).
///
/// The Acts finders assume a uniform field along their z axis, layers ordered in radius and
/// tracks from a beam line. The B0 tracker is four planar stations inside a dipole with a
/// quadrupole gradient, whose field turns by atan(G y / B) off the midplane. So for each middle
/// hit the hits are rotated into a frame whose z' axis is along the field at that hit and whose
/// x' axis is along the ion beam: the stations become "layers" at r' = distance along the beam,
/// and the field is a solenoid-like field along z' of the magnitude there (finders are built for
/// bins of |B|). The beam-line cuts and seed weights are opened up, since B0 tracks may come from
/// displaced decays. In a frame fixed along global y instead, the r-z' compatibility cut loses
/// tracks in proportion to their distance from the midplane.
///
/// The hit-time window, the seed parameters in the effective field along the triplet, the
/// anchor and the prior are those of B0TripletSeeding (shared helpers and configuration), so for
/// the same hit triplet both write identical seeds and the track finding downstream is the same.
class B0Seeding2 : public B0Seeding2Algorithm, public WithPodConfig<B0Seeding2Config> {
public:
  /// Whether this Acts build has the Seeding2 triplet seeder; init() throws otherwise
  static constexpr bool available() { return B0SEEDING2_HAS_SEEDING2; }

  B0Seeding2(std::string_view name)
      : B0Seeding2Algorithm{name,
                            {"inputTrackerHits"},
                            {"outputTrackSeeds", "outputTrackParameters"},
                            "create B0 track seeds with the Acts Seeding2 triplet seeder"} {}
  void init() final;
  void process(const Input& input, const Output& output) const final;

  /// Rotation from global to the primed frame: rows are x', y', z'
  static Acts::RotationMatrix3 frameRotation(const Acts::Vector3& fieldAxis,
                                             const Acts::Vector3& beamAxis);

  struct Impl;

private:
  const algorithms::ActsSvc& m_actsSvc{algorithms::ActsSvc::instance()};
  const std::shared_ptr<const ActsGeometryProvider> m_geoSvc{m_actsSvc.acts_geometry_provider()};
  std::shared_ptr<const Impl> m_impl;
};

} // namespace eicrecon
