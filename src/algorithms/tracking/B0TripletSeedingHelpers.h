// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "algorithms/tracking/B0TripletSeedingConfig.h"

/// Steps of B0TripletSeeding, separated from the algorithm so that they can be tested without
/// the geometry service
namespace eicrecon::b0_triplet_seeding {

/// A hit in Acts units, with the magnetic field at its position (zero if the lookup failed)
struct Point {
  Acts::Vector3 position;
  Acts::Vector3 field;
  double time;
  double variance;
  std::size_t index;
};

/// A triplet of hits on three stations and the distance of its middle hit from the chord
struct Candidate {
  double residual;
  std::array<const Point*, 3> points;
};

/// Seed parameters of a triplet and their diagonal prior, on a perigee surface at the anchor
struct Seed {
  Acts::Vector3 anchor;
  Acts::BoundVector parameters;
  Acts::BoundVector errors;
};

/// Sorts the hits with a finite position and time by z and splits them into stations at gaps
/// larger than stationGap
std::vector<std::vector<Point>> groupStations(const edm4eic::TrackerHitCollection& hits,
                                              double stationGap,
                                              const Acts::MagneticFieldProvider& field,
                                              Acts::MagneticFieldProvider::Cache& fieldCache);

/// Distance of the middle hit from the chord of the outer two, along the field at the middle hit
double chordResidual(const Point& a, const Point& b, const Point& c);

/// Seed parameters of a triplet in the field at its middle hit, if the triplet passes the
/// momentum limit and its prior is positive definite
std::optional<Seed> makeSeed(const B0TripletSeedingConfig& cfg, const Acts::GeometryContext& gctx,
                             const Point& a, const Point& b, const Point& c);

/// Writes a seed and its parameters, converted to EDM4eic units
void writeSeed(const Seed& seed, const Candidate& candidate,
               const edm4eic::TrackerHitCollection& hits, edm4eic::TrackSeedCollection& trackSeeds,
               edm4eic::TrackParametersCollection& trackParams);

} // namespace eicrecon::b0_triplet_seeding
