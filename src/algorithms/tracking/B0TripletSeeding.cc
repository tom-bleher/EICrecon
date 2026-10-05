// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0TripletSeeding.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/TransformationHelpers.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "ActsGeometryProvider.h"
#include "B0TripletSeedingConfig.h"
#include "extensions/edm4eic/EDM4eicToActs.h"

namespace eicrecon {

namespace {

  /// A hit in Acts units, with the magnetic field at its position (zero if the lookup failed)
  struct Point {
    Acts::Vector3 position;
    Acts::Vector3 field;
    double time;
    double variance;
    std::size_t index;
    std::size_t rank;
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

  struct SeedCandidate {
    Candidate candidate;
    Seed seed;
  };

  bool candidateLess(const Candidate& lhs, const Candidate& rhs) {
    return std::tie(lhs.residual, lhs.points[0]->rank, lhs.points[1]->rank, lhs.points[2]->rank) <
           std::tie(rhs.residual, rhs.points[0]->rank, rhs.points[1]->rank, rhs.points[2]->rank);
  }

  /// Sorts the hits with a finite position and time by z and splits them into stations at gaps
  /// larger than stationGap. Each covariance component is validated
  /// individually: a maximum can hide a negative component or ignore a NaN.
  std::vector<std::vector<Point>> groupStations(const edm4eic::TrackerHitCollection& hits,
                                                double stationGap,
                                                const Acts::MagneticFieldProvider& field,
                                                Acts::MagneticFieldProvider::Cache& fieldCache) {
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < hits.size(); ++i) {
      const auto pos = hits[i].getPosition();
      const auto cov = hits[i].getPositionError();
      if (std::isfinite(pos.x) && std::isfinite(pos.y) && std::isfinite(pos.z) &&
          std::isfinite(hits[i].getTime()) && std::isfinite(cov.xx) && cov.xx >= 0.f &&
          std::isfinite(cov.yy) && cov.yy >= 0.f && std::isfinite(cov.zz) && cov.zz >= 0.f) {
        order.push_back(i);
      }
    }
    auto sortKey = [&](std::size_t i) {
      const auto pos = hits[i].getPosition();
      const auto cov = hits[i].getPositionError();
      return std::tuple{pos.z,
                        pos.x,
                        pos.y,
                        hits[i].getCellID(),
                        hits[i].getTime(),
                        std::max({cov.xx, cov.yy, cov.zz})};
    };
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return sortKey(a) < sortKey(b); });

    std::vector<std::vector<Point>> stations;
    double lastZ = -std::numeric_limits<double>::infinity();
    for (std::size_t rank = 0; rank < order.size(); ++rank) {
      const auto i   = order[rank];
      const auto pos = hits[i].getPosition();
      const Acts::Vector3 position{pos.x / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                   pos.y / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                   pos.z / edm4eic::unit::mm * Acts::UnitConstants::mm};
      const double time = hits[i].getTime() / edm4eic::unit::ns * Acts::UnitConstants::ns;
      // The hit covariance is diagonal in local sensor coordinates; use its
      // largest variance as an isotropic position variance
      const auto cov        = hits[i].getPositionError();
      const double variance = std::max({cov.xx, cov.yy, cov.zz}) /
                              (edm4eic::unit::mm * edm4eic::unit::mm) *
                              (Acts::UnitConstants::mm * Acts::UnitConstants::mm);
      const auto hitField = field.getField(position, fieldCache);
      if (stations.empty() || position.z() - lastZ > stationGap) {
        stations.emplace_back();
      }
      lastZ = position.z();
      stations.back().push_back(
          {position, hitField.ok() ? *hitField : Acts::Vector3::Zero(), time, variance, i, rank});
    }
    return stations;
  }

  /// Distance of the middle hit from the chord of the outer two, along the field at the middle hit
  double chordResidual(const Point& a, const Point& b, const Point& c) {
    const Acts::Vector3 chord = c.position - a.position;
    const double fraction     = (b.position.z() - a.position.z()) / chord.z();
    return std::abs((b.position - a.position - fraction * chord).dot(b.field.normalized()));
  }

  /// Seed parameters of a triplet in the field at its middle hit, if the triplet passes the
  /// momentum limit and its prior is positive definite
  std::optional<Seed> makeSeed(const B0TripletSeedingConfig& cfg, const Acts::GeometryContext& gctx,
                               const Point& a, const Point& b, const Point& c) {
    const Acts::Vector3& bField = b.field;
#if Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 2)
    Acts::FreeVector freeParams =
        Acts::estimateTrackParamsFromSeed(a.position, a.time, b.position, c.position, bField);
#else
    Acts::FreeVector freeParams =
        Acts::estimateTrackParamsFromSeed(a.position, b.position, c.position, bField);
    freeParams[Acts::eFreeTime] = a.time;
#endif
    if (!freeParams.allFinite() || std::abs(freeParams[Acts::eFreeQOverP]) * cfg.minMomentum > 1.) {
      return std::nullopt;
    }

    // A relative q/p error vanishes for a straight triplet. The three hits
    // resolve the curvature only up to their sagitta error over the lever arms
    // across the field, and |q/p| = curvature * sin(angle to field) / |B|.
    const Acts::Vector3 fieldDirection = bField.normalized();
    auto across                        = [&](const Acts::Vector3& v) {
      return (v - v.dot(fieldDirection) * fieldDirection).norm();
    };
    const double d1 = across(b.position - a.position);
    const double d2 = across(c.position - b.position);
    const double sagittaVariance =
        b.variance + (d2 * d2 * a.variance + d1 * d1 * c.variance) / ((d1 + d2) * (d1 + d2));
    const Acts::Vector3 direction = freeParams.segment<3>(Acts::eFreeDir0);
    const double qOverPResolution = 2. * std::sqrt(sagittaVariance) / (d1 * d2) *
                                    direction.cross(fieldDirection).norm() / bField.norm();

    // Move the parameters estimated at the first hit back along the tangent to a
    // perigee anchorDistance upstream, where CKFTracking starts the track finding
    const Acts::Vector3 anchor             = a.position - cfg.anchorDistance * direction;
    freeParams.segment<3>(Acts::eFreePos0) = anchor;
    // Acts measures time in length units (c = 1); the shift assumes beta = 1
    freeParams[Acts::eFreeTime] -= cfg.anchorDistance;
    const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(anchor);
    const auto local   = Acts::transformFreeToBoundParameters(freeParams, *perigee, gctx);
    if (!local.ok() || !local->allFinite()) {
      return std::nullopt;
    }

    // The azimuth and the position along the perigee line of a forward track
    // are poorly defined, so scale their errors
    const double theta = (*local)[Acts::eBoundTheta];
    Acts::BoundVector errors;
    errors << cfg.positionError, cfg.positionError / std::tan(theta),
        cfg.angleError / std::sin(theta), cfg.angleError,
        std::max(cfg.qOverPRelativeError * std::abs((*local)[Acts::eBoundQOverP]),
                 qOverPResolution),
        cfg.timeError;
    if (!errors.allFinite() || (errors.array() <= 0.).any()) {
      return std::nullopt;
    }
    return Seed{.anchor = anchor, .parameters = *local, .errors = errors};
  }

  /// Writes a seed and its parameters, converted to EDM4eic units
  void writeSeed(const Seed& seed, const Candidate& candidate,
                 const edm4eic::TrackerHitCollection& hits,
                 edm4eic::TrackSeedCollection& trackSeeds,
                 edm4eic::TrackParametersCollection& trackParams) {
    const auto& parameter = seed.parameters;
    auto pars             = trackParams.create();
    pars.setType(-1); // type --> seed(-1)
    pars.setLoc({static_cast<float>(parameter[Acts::eBoundLoc0] / Acts::UnitConstants::mm *
                                    edm4eic::unit::mm),
                 static_cast<float>(parameter[Acts::eBoundLoc1] / Acts::UnitConstants::mm *
                                    edm4eic::unit::mm)});
    pars.setPhi(static_cast<float>(parameter[Acts::eBoundPhi] / Acts::UnitConstants::rad));
    pars.setTheta(static_cast<float>(parameter[Acts::eBoundTheta] / Acts::UnitConstants::rad));
    pars.setQOverP(static_cast<float>(parameter[Acts::eBoundQOverP] * Acts::UnitConstants::GeV /
                                      edm4eic::unit::GeV));
    pars.setTime(static_cast<float>(parameter[Acts::eBoundTime] / Acts::UnitConstants::ns *
                                    edm4eic::unit::ns));
    // The prior is diagonal
    edm4eic::Cov6f cov;
    for (std::size_t i = 0; const auto& [ia, x] : edm4eic_indexed_units) {
      cov(i, i) = seed.errors[ia] * seed.errors[ia] / (x * x);
      ++i;
    }
    pars.setCovariance(cov);

    auto trackSeed = trackSeeds.create();
    trackSeed.setPerigee(
        {static_cast<float>(seed.anchor.x() / Acts::UnitConstants::mm * edm4eic::unit::mm),
         static_cast<float>(seed.anchor.y() / Acts::UnitConstants::mm * edm4eic::unit::mm),
         static_cast<float>(seed.anchor.z() / Acts::UnitConstants::mm * edm4eic::unit::mm)});
    // The quality is the negative residual in mm, so that a higher quality is a better seed
    trackSeed.setQuality(static_cast<float>(-candidate.residual / Acts::UnitConstants::mm));
    trackSeed.setParams(pars);
    for (const auto* point : candidate.points) {
      trackSeed.addToHits(hits[point->index]);
    }
  }

} // namespace

void B0TripletSeeding::init() { validateConfig(m_cfg); }

void B0TripletSeeding::validateConfig(const B0TripletSeedingConfig& cfg) {
  // Positive comparisons reject NaN; isfinite additionally rejects infinity.
  // No parameter permits an infinite value.
  // NOLINTNEXTLINE(readability-simplify-boolean-expr)
  if (!(std::isfinite(cfg.stationGap) && cfg.stationGap > 0. && std::isfinite(cfg.maxResidual) &&
        cfg.maxResidual > 0. && std::isfinite(cfg.minMomentum) && cfg.minMomentum > 0. &&
        std::isfinite(cfg.anchorDistance) && cfg.anchorDistance >= 0. && cfg.maxSeeds > 0 &&
        std::isfinite(cfg.positionError) && cfg.positionError > 0. &&
        std::isfinite(cfg.angleError) && cfg.angleError > 0. &&
        std::isfinite(cfg.qOverPRelativeError) && cfg.qOverPRelativeError > 0. &&
        std::isfinite(cfg.timeError) && cfg.timeError > 0.)) {
    throw std::runtime_error("B0TripletSeeding: configuration values must be positive and finite "
                             "(anchorDistance may be zero)");
  }
}

void B0TripletSeeding::process(const Input& input, const Output& output) const {
  const auto [hits]                = input;
  auto [track_seeds, track_params] = output;

  const auto& gctx = m_geoSvc->getActsGeometryContext();
  const auto& mctx = m_geoSvc->getActsMagneticFieldContext();
  const auto field = m_geoSvc->getFieldProvider();
  const auto stats = seedHits(m_cfg, gctx, mctx, *field, *hits, *track_seeds, *track_params);
  debug("{} hits, {} stations, {} seeds", hits->size(), stats.stations, track_seeds->size());
  debug("triplets enumerated={}, residual-passing={}, rank-pruned={}, estimated={}, "
        "estimation-failures={}, heap-replacements={}, retained={}",
        stats.enumerated, stats.residualPassing, stats.rankPruned, stats.estimated,
        stats.estimationFailures, stats.heapReplacements, stats.retained);
}

B0TripletSeeding::Stats B0TripletSeeding::seedHits(
    const B0TripletSeedingConfig& cfg, const Acts::GeometryContext& gctx,
    const Acts::MagneticFieldContext& mctx, const Acts::MagneticFieldProvider& field,
    const edm4eic::TrackerHitCollection& hits, edm4eic::TrackSeedCollection& trackSeeds,
    edm4eic::TrackParametersCollection& trackParams) {
  validateConfig(cfg);
  Stats stats;
  if (hits.empty() || trackSeeds.size() >= cfg.maxSeeds) {
    return stats;
  }

  auto fieldCache = field.makeCache(mctx);

  const auto stations = groupStations(hits, cfg.stationGap, field, fieldCache);
  stats.stations      = stations.size();

  // Triplets of hits on three stations. In a uniform field the middle hit of a
  // forward track lies close to the chord of the outer two along the field. The
  // quadrupole component of the B0 magnet breaks this off its midplane, so a
  // soft track whose triplet skips a station can exceed maxResidual.
  // Keep only the best valid seeds. The heap's front is the worst retained
  // candidate, allowing worse triplets to skip parameter estimation. Triplet
  // enumeration is still cubic in station occupancy, but storage is bounded.
  const auto limit = cfg.maxSeeds - trackSeeds.size();
  std::vector<SeedCandidate> candidates;
  candidates.reserve(limit);
  const auto less = [](const SeedCandidate& lhs, const SeedCandidate& rhs) {
    return candidateLess(lhs.candidate, rhs.candidate);
  };
  for (std::size_t i = 0; i < stations.size(); ++i) {
    for (std::size_t j = i + 1; j < stations.size(); ++j) {
      for (std::size_t k = j + 1; k < stations.size(); ++k) {
        for (const auto& a : stations[i]) {
          for (const auto& b : stations[j]) {
            if (b.field.norm() == 0.) {
              continue;
            }
            for (const auto& c : stations[k]) {
              ++stats.enumerated;
              const double residual = chordResidual(a, b, c);
              if (!(residual < cfg.maxResidual)) {
                continue;
              }
              ++stats.residualPassing;
              const Candidate candidate{residual, {&a, &b, &c}};
              if (candidates.size() == limit &&
                  !candidateLess(candidate, candidates.front().candidate)) {
                ++stats.rankPruned;
                continue;
              }
              ++stats.estimated;
              if (const auto seed = makeSeed(cfg, gctx, a, b, c)) {
                if (candidates.size() == limit) {
                  std::pop_heap(candidates.begin(), candidates.end(), less);
                  candidates.pop_back();
                  ++stats.heapReplacements;
                }
                candidates.push_back({candidate, *seed});
                std::push_heap(candidates.begin(), candidates.end(), less);
              } else {
                ++stats.estimationFailures;
              }
            }
          }
        }
      }
    }
  }
  std::sort_heap(candidates.begin(), candidates.end(), less);
  for (const auto& candidate : candidates) {
    writeSeed(candidate.seed, candidate.candidate, hits, trackSeeds, trackParams);
  }
  stats.retained = candidates.size();
  return stats;
}

} // namespace eicrecon
