// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0TripletSeeding.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/TransformationHelpers.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Result.hpp>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <fmt/format.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <compare>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "ActsGeometryProvider.h"
#include "B0TripletSeedingConfig.h"
#include "B0TripletSeedingHelpers.h"
#include "extensions/edm4eic/EDM4eicToActs.h"

namespace eicrecon {

namespace b0_triplet_seeding {

  std::vector<std::vector<Point>> groupStations(const edm4eic::TrackerHitCollection& hits,
                                                double stationGap) {
    // The stations are found from z gaps, not from the cellID. The layer field of the B0 readout
    // means different things in different geometries: one layer per station (1-4) in the current
    // geometry, a front and a back layer per station (1-8) in the realistic one (epic PR #1194),
    // so decoding it would tie the seeder to one of them. The stations are at least 270 mm apart
    // and at most 12 mm deep in both, so a gap in z separates them in both.
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < hits.size(); ++i) {
      const auto pos = hits[i].getPosition();
      if (std::isfinite(pos.x) && std::isfinite(pos.y) && std::isfinite(pos.z) &&
          std::isfinite(hits[i].getTime())) {
        order.push_back(i);
      }
    }
    auto sortKey = [&](std::size_t i) {
      const auto pos = hits[i].getPosition();
      return std::tuple{pos.z, pos.x, pos.y, hits[i].getCellID()};
    };
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return sortKey(a) < sortKey(b); });

    std::vector<std::vector<Point>> stations;
    double lastZ = -std::numeric_limits<double>::infinity();
    for (const auto i : order) {
      const auto pos = hits[i].getPosition();
      const Acts::Vector3 position{pos.x / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                   pos.y / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                   pos.z / edm4eic::unit::mm * Acts::UnitConstants::mm};
      const double time = hits[i].getTime() / edm4eic::unit::ns * Acts::UnitConstants::ns;
      // The variance enters the sagitta error, across the field in the sensor plane. The hit
      // covariance is diagonal in local sensor coordinates (xx, yy in the plane; zz along the
      // normal, zero for 2-D pixels). The local axes are not aligned with the bending direction in
      // general, so take the larger in-plane variance; for square pixels, as in the B0 tracker,
      // both are equal and this is exact. The normal is nearly along the track, where a shift does
      // not change the sagitta, so zz is left out.
      const auto cov        = hits[i].getPositionError();
      const double variance = std::max(cov.xx, cov.yy) / (edm4eic::unit::mm * edm4eic::unit::mm) *
                              (Acts::UnitConstants::mm * Acts::UnitConstants::mm);
      if (stations.empty() || position.z() - lastZ > stationGap) {
        stations.emplace_back();
      }
      lastZ = position.z();
      stations.back().push_back({.position       = position,
                                 .field          = Acts::Vector3::Zero(),
                                 .fieldDirection = Acts::Vector3::Zero(),
                                 .time           = time,
                                 .variance       = variance,
                                 .index          = i});
    }
    return stations;
  }

  std::size_t lookUpField(std::vector<std::vector<Point>>& stations,
                          const Acts::MagneticFieldProvider& field,
                          Acts::MagneticFieldProvider::Cache& fieldCache) {
    std::size_t failed = 0;
    for (std::size_t j = 1; j + 1 < stations.size(); ++j) {
      for (auto& point : stations[j]) {
        const auto hitField = field.getField(point.position, fieldCache);
        if (hitField.ok() && hitField->norm() > 0.) {
          point.field          = *hitField;
          point.fieldDirection = hitField->normalized();
        } else if (!hitField.ok()) {
          ++failed;
        }
      }
    }
    return failed;
  }

  double chordResidual(const Point& a, const Point& b, const Point& c) {
    const Acts::Vector3 chord = c.position - a.position;
    const double fraction     = (b.position.z() - a.position.z()) / chord.z();
    return std::abs((b.position - a.position - fraction * chord).dot(b.fieldDirection));
  }

  std::vector<Candidate> findCandidates(const B0TripletSeedingConfig& cfg,
                                        const std::vector<std::vector<Point>>& stations) {
    // Triplets of hits on three stations. In a uniform field the middle hit of a
    // forward track lies close to the chord of the outer two along the field. The
    // quadrupole component of the B0 magnet breaks this off its midplane, so a
    // soft track whose triplet skips a station can exceed maxResidual.
    // Hits of one track are separated in time by their flight time, in Acts units (c = 1) the
    // distance between them. Checking the first hit against the middle one before the loop over
    // the last hit rejects out-of-time hits cheaply.
    const double maxTimeResidual =
        cfg.maxTimeResidual / edm4eic::unit::ns * Acts::UnitConstants::ns;
    auto inTime = [&](const Point& first, const Point& later) {
      return maxTimeResidual <= 0. ||
             std::abs(later.time - first.time - (later.position - first.position).norm()) <
                 maxTimeResidual;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < stations.size(); ++i) {
      for (std::size_t j = i + 1; j < stations.size(); ++j) {
        for (std::size_t k = j + 1; k < stations.size(); ++k) {
          for (const auto& a : stations[i]) {
            for (const auto& b : stations[j]) {
              if (b.fieldDirection.isZero() || !inTime(a, b)) {
                continue;
              }
              for (const auto& c : stations[k]) {
                if (!inTime(a, c)) {
                  continue;
                }
                const double residual = chordResidual(a, b, c);
                if (residual < cfg.maxResidual) {
                  candidates.push_back({residual, {&a, &b, &c}});
                }
              }
            }
          }
        }
      }
    }
    return candidates;
  }

  std::size_t selectSeeds(const B0TripletSeedingConfig& cfg, std::vector<Candidate>& candidates,
                          const std::function<bool(const Candidate&)>& makeSeed) {
    // The hit indices make the order total, so the seeds do not depend on the order of the hits
    auto ranksBefore = [](const Candidate& lhs, const Candidate& rhs) {
      return std::tie(lhs.residual, lhs.points[0]->index, lhs.points[1]->index,
                      lhs.points[2]->index) < std::tie(rhs.residual, rhs.points[0]->index,
                                                       rhs.points[1]->index, rhs.points[2]->index);
    };
    // A heap with the best candidate on top: building it is linear, and each candidate taken costs
    // a logarithm, so a busy event is not fully sorted only to be cut at maxSeeds
    auto ranksAfter = [&](const Candidate& lhs, const Candidate& rhs) {
      return ranksBefore(rhs, lhs);
    };
    std::ranges::make_heap(candidates, ranksAfter);
    // Hits of the seeds made so far, by hit index
    std::vector<bool> used;
    auto isUsed = [&](const Point* point) {
      return point->index < used.size() && used[point->index];
    };
    auto end          = candidates.end();
    std::size_t seeds = 0;
    while (end != candidates.begin() && seeds < cfg.maxSeeds) {
      std::ranges::pop_heap(candidates.begin(), end, ranksAfter);
      --end;
      if (cfg.dropSeedsWithUsedHits && std::ranges::all_of(end->points, isUsed)) {
        continue;
      }
      if (makeSeed(*end)) {
        ++seeds;
        if (cfg.dropSeedsWithUsedHits) {
          for (const auto* point : end->points) {
            used.resize(std::max(used.size(), point->index + 1));
            used[point->index] = true;
          }
        }
      }
    }
    return static_cast<std::size_t>(std::distance(candidates.begin(), end));
  }

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
    const Acts::Vector3& fieldDirection = b.fieldDirection;
    auto across                         = [&](const Acts::Vector3& v) {
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
    // Acts measures time in length units (c = 1); the shift assumes beta = 1. For a 1 GeV proton
    // (beta = 0.73) it is off by 0.37 anchorDistance, 0.01 ns at the default, far below timeError.
    freeParams[Acts::eFreeTime] -= cfg.anchorDistance;
    // The track bends away from the straight tangent by about the sagitta anchorDistance^2 / (2 R),
    // with 1 / R = |q/p| |direction x B| in Acts units. Add it to the position error, so that the
    // prior covers it for any anchorDistance; at the default of 10 mm it is at most 0.02 mm.
    const double anchorSagitta = 0.5 * cfg.anchorDistance * cfg.anchorDistance *
                                 std::abs(freeParams[Acts::eFreeQOverP]) *
                                 direction.cross(bField).norm();
    const double positionError = std::hypot(cfg.positionError, anchorSagitta);

    // A perigee surface, whose line is along global z, because a TrackSeed carries only a perigee
    // point and CKFTracking builds a perigee surface from it, and reports the track there. A
    // forward track is nearly parallel to that line, so its position along the line and its
    // azimuth are poorly defined; the prior below scales their errors to match. A plane across the
    // track would be better conditioned, but needs changes to the TrackSeed data model and to
    // CKFTracking, which the central tracking shares.
    const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(anchor);
    const auto local   = Acts::transformFreeToBoundParameters(freeParams, *perigee, gctx);
    if (!local.ok() || !local->allFinite()) {
      return std::nullopt;
    }

    // The azimuth and the position along the perigee line of a forward track
    // are poorly defined, so scale their errors
    const double theta = (*local)[Acts::eBoundTheta];
    Acts::BoundVector errors;
    errors << positionError, positionError / std::abs(std::tan(theta)),
        cfg.angleError / std::sin(theta), cfg.angleError,
        std::max(cfg.qOverPRelativeError * std::abs((*local)[Acts::eBoundQOverP]),
                 qOverPResolution),
        cfg.timeError / edm4eic::unit::ns * Acts::UnitConstants::ns;
    if (!errors.allFinite() || (errors.array() <= 0.).any()) {
      return std::nullopt;
    }
    return Seed{.anchor = anchor, .parameters = *local, .errors = errors};
  }

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

} // namespace b0_triplet_seeding

void B0TripletSeeding::init() {
  struct Named {
    std::string_view name;
    double value;
  };
  // The comparisons are written so that NaN values fail them too
  const std::array<Named, 7> positive{{
      {"stationGap", m_cfg.stationGap},
      {"maxResidual", m_cfg.maxResidual},
      {"minMomentum", m_cfg.minMomentum},
      {"positionError", m_cfg.positionError},
      {"angleError", m_cfg.angleError},
      {"qOverPRelativeError", m_cfg.qOverPRelativeError},
      {"timeError", m_cfg.timeError},
  }};
  for (const auto& [name, value] : positive) {
    if (!(value > 0.)) {
      throw std::runtime_error(
          fmt::format("B0TripletSeeding: {} must be positive, got {}", name, value));
    }
  }
  if (!(m_cfg.maxTimeResidual >= 0.)) {
    throw std::runtime_error(fmt::format(
        "B0TripletSeeding: maxTimeResidual must not be negative, got {}", m_cfg.maxTimeResidual));
  }
  if (!(m_cfg.anchorDistance >= 0.)) {
    throw std::runtime_error(fmt::format(
        "B0TripletSeeding: anchorDistance must not be negative, got {}", m_cfg.anchorDistance));
  }
  if (m_cfg.maxSeeds == 0) {
    throw std::runtime_error("B0TripletSeeding: maxSeeds must be positive");
  }
  if (!m_geoSvc || !m_geoSvc->getFieldProvider()) {
    throw std::runtime_error("B0TripletSeeding: no magnetic field provider");
  }
}

void B0TripletSeeding::process(const Input& input, const Output& output) const {
  using namespace b0_triplet_seeding;

  const auto [hits]                = input;
  auto [track_seeds, track_params] = output;

  if (hits->empty()) {
    return;
  }

  const auto& gctx = m_geoSvc->getActsGeometryContext();
  const auto& mctx = m_geoSvc->getActsMagneticFieldContext();
  const auto field = m_geoSvc->getFieldProvider();
  auto fieldCache  = field->makeCache(mctx);

  auto stations = groupStations(*hits, m_cfg.stationGap);
  // A hit without a field cannot be a middle hit. If the field map does not cover the B0 tracker,
  // for example with a mismatched geometry, no seeds are made, so say so instead of failing quietly
  if (const std::size_t failed = lookUpField(stations, *field, fieldCache); failed > 0) {
    const std::size_t before = m_fieldLookupFailures.fetch_add(failed);
    if (std::bit_width(before) < std::bit_width(before + failed)) {
      warning(
          "The magnetic field lookup failed at {} B0 hits in this event ({} so far); those hits "
          "cannot be the middle hit of a seed",
          failed, before + failed);
    }
  }

  auto candidates = findCandidates(m_cfg, stations);
  // A track crossing all stations gives one candidate per station triple, and all of them are kept.
  // The track finding does not succeed from every seed of a track, so keeping only the best one
  // per track loses tracks; the ambiguity solver removes the duplicate tracks instead.
  const std::size_t dropped = selectSeeds(m_cfg, candidates, [&](const Candidate& candidate) {
    const auto& [a, b, c] = candidate.points;
    const auto seed       = makeSeed(m_cfg, gctx, *a, *b, *c);
    if (seed) {
      writeSeed(*seed, candidate, *hits, *track_seeds, *track_params);
    }
    return seed.has_value();
  });
  if (dropped > 0) {
    // The candidates are ranked by residual alone, so genuine seeds can be among those dropped
    const std::size_t truncated = m_truncatedEvents.fetch_add(1) + 1;
    if (std::has_single_bit(truncated)) {
      warning(
          "Reached maxSeeds={}: dropped the {} lowest-ranked of {} seed candidates in this event "
          "({} events so far); genuine seeds can be among them",
          m_cfg.maxSeeds, dropped, candidates.size(), truncated);
    } else {
      debug("Reached maxSeeds={}: dropped the {} lowest-ranked of {} seed candidates",
            m_cfg.maxSeeds, dropped, candidates.size());
    }
  }

  debug("{} stations, {} triplet candidates, {} seeds", stations.size(), candidates.size(),
        track_seeds->size());
}

} // namespace eicrecon
