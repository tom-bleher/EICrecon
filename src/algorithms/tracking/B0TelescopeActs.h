// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#pragma once

// Acts_VERSION_* is supplied by the Acts CMake target, as in TrackSeeding.h.
#define EICRECON_HAS_B0_TELESCOPE \
  (Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 3))

#if EICRECON_HAS_B0_TELESCOPE
#if Acts_VERSION_MAJOR >= 47
#include <Acts/EventData/SeedContainer.hpp>
#include <Acts/EventData/SpacePointContainer.hpp>
#include <Acts/Seeding/DoubletSeedFinder.hpp>
#include <Acts/Seeding/ITripletSeedFilter.hpp>
#include <Acts/Seeding/TripletSeedFinder.hpp>
#include <Acts/Seeding/TripletSeeder.hpp>
#else
#include <Acts/EventData/SeedContainer2.hpp>
#include <Acts/EventData/SpacePointContainer2.hpp>
#include <Acts/Seeding2/DoubletSeedFinder.hpp>
#include <Acts/Seeding2/ITripletSeedFilter.hpp>
#include <Acts/Seeding2/TripletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeeder.hpp>
#endif

#include <type_traits>

namespace eicrecon::b0telescope {
#if Acts_VERSION_MAJOR >= 47
using SpacePoints = Acts::SpacePointContainer;
using SpacePoint = Acts::ConstSpacePointProxy;
using Seeds = Acts::SeedContainer;
#else
using SpacePoints = Acts::SpacePointContainer2;
using SpacePoint = Acts::ConstSpacePointProxy2;
using Seeds = Acts::SeedContainer2;
#endif
// Released ACTS uses mutable views; development ACTS passes views by value and
// returns triplet cursors. Detect the signatures instead of guessing a version.
template <typename Finder>
inline constexpr bool mutableDoubletViews = requires {
  static_cast<void (Finder::*)(const SpacePoint&, const Acts::MiddleSpInfo&,
                               SpacePoints::ConstSubset&, Acts::DoubletsForMiddleSp&) const>(
      &Finder::createDoublets);
};
template <typename Finder>
inline constexpr bool mutableTripletViews = requires {
  static_cast<void (Finder::*)(const SpacePoints&, const SpacePoint&,
                               const Acts::DoubletsForMiddleSp::Proxy&,
                               Acts::DoubletsForMiddleSp::Range&,
                               Acts::TripletTopCandidates&) const>(
      &Finder::createTripletTopCandidates);
};
template <typename View>
using DoubletView = std::conditional_t<mutableDoubletViews<Acts::DoubletSeedFinder>, View&, View>;
template <typename View>
using TripletView = std::conditional_t<mutableTripletViews<Acts::TripletSeedFinder>, View&, View>;
template <typename View>
using TripletReturn = std::conditional_t<mutableTripletViews<Acts::TripletSeedFinder>, void, View>;
} // namespace eicrecon::b0telescope
#endif
