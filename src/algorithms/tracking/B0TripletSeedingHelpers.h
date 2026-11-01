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
#include <functional>
#include <optional>
#include <vector>

#include "algorithms/tracking/B0TripletSeedingConfig.h"

/// Steps of B0TripletSeeding, separated from the algorithm so that they can be tested without
/// the geometry service
namespace eicrecon::b0_triplet_seeding {

/// A hit in Acts units. The magnetic field at its position and the field direction are only
/// looked up for hits that can be the middle hit of a triplet, and are zero otherwise or if the
/// lookup failed.
struct Point {
  Acts::Vector3 position;
  Acts::Vector3 field;
  Acts::Vector3 fieldDirection;
  double time;
  double variance;
  std::size_t index;
  double timeError = 0.;
};

/// A triplet of hits on three stations, the distance of its middle hit from the chord, the
/// effective magnetic field and the q/p estimated in it
struct Candidate {
  double residual;
  std::array<const Point*, 3> points;
  Acts::Vector3 field = Acts::Vector3::Zero();
  double qOverP       = 0.;
};

/// The magnetic field at a position, zero if the lookup fails
using FieldLookup = std::function<Acts::Vector3(const Acts::Vector3&)>;

/// Effective field of a triplet, and the residual that the field change along it predicts per
/// unit q/p
struct FieldIntegral {
  Acts::Vector3 field;
  double residualPerQOverP;
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
                                              double stationGap);

/// Looks up the magnetic field at the hits of the stations that can hold the middle hit of a
/// triplet (all but the first and the last), and returns the number of failed lookups
std::size_t lookUpField(std::vector<std::vector<Point>>& stations,
                        const Acts::MagneticFieldProvider& field,
                        Acts::MagneticFieldProvider::Cache& fieldCache);

/// Signed distance of the middle hit from the chord of the outer two, along the field at the
/// middle hit, and its magnitude
double signedChordResidual(const Point& a, const Point& b, const Point& c);
double chordResidual(const Point& a, const Point& b, const Point& c);

/// Samples the field on the parabola through the three hits, at two Gauss-Legendre points in each
/// of (samples + 3) / 4 intervals on either side of the middle hit, and weights each by its effect
/// on the position of the middle hit relative to the chord (a triangle in z, zero at the outer
/// hits). The weighted mean is the effective field. The curvature d2r/dz2 it causes along the
/// field at the middle hit, integrated with the same weights, predicts the chord residual. Falls
/// back to the field at the middle hit if all samples are zero.
FieldIntegral fieldIntegral(const Point& a, const Point& b, const Point& c, std::size_t samples,
                            const FieldLookup& fieldAt);

/// Whether the later hit can follow the first one in time: it may lag by up to the flight time
/// of the slowest seeded particle over their distance and may not be early, both within the time
/// window tolerance. Always true when the window is off.
bool inTimeWindow(const B0TripletSeedingConfig& cfg, const Point& first, const Point& later);

/// All triplets of hits on three different stations, in time with each other, with an estimated
/// momentum of at least minMomentum, whose middle hit lies within maxResidual of the chord of the
/// outer two after subtracting the predicted field term, in no particular order. Without a field
/// lookup, or with fieldSamples = 0, the field at the middle hit is used and nothing is subtracted.
std::vector<Candidate> findCandidates(const B0TripletSeedingConfig& cfg,
                                      const std::vector<std::vector<Point>>& stations,
                                      const FieldLookup& fieldAt = {});

/// Passes the candidates to makeSeed in order of increasing residual, ties broken by the hit
/// indices, until makeSeed has returned true maxSeeds times. Only the candidates that are reached
/// are ordered. Returns the number of candidates that were not reached.
std::size_t selectSeeds(const B0TripletSeedingConfig& cfg, std::vector<Candidate>& candidates,
                        const std::function<bool(const Candidate&)>& makeSeed);

/// Seed parameters of a triplet in a field, by default the one at its middle hit, if the triplet
/// passes the momentum limit and its prior is positive definite
std::optional<Seed> makeSeed(const B0TripletSeedingConfig& cfg, const Acts::GeometryContext& gctx,
                             const Point& a, const Point& b, const Point& c,
                             const Acts::Vector3& field);
std::optional<Seed> makeSeed(const B0TripletSeedingConfig& cfg, const Acts::GeometryContext& gctx,
                             const Point& a, const Point& b, const Point& c);

/// Writes a seed and its parameters, converted to EDM4eic units
void writeSeed(const Seed& seed, const Candidate& candidate,
               const edm4eic::TrackerHitCollection& hits, edm4eic::TrackSeedCollection& trackSeeds,
               edm4eic::TrackParametersCollection& trackParams);

} // namespace eicrecon::b0_triplet_seeding
