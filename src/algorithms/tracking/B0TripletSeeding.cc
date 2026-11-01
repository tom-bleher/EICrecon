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
      const double time      = hits[i].getTime() / edm4eic::unit::ns * Acts::UnitConstants::ns;
      const double timeError = hits[i].getTimeError() / edm4eic::unit::ns * Acts::UnitConstants::ns;
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
                                 .index          = i,
                                 .timeError      = timeError});
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

  double signedChordResidual(const Point& a, const Point& b, const Point& c) {
    const Acts::Vector3 chord = c.position - a.position;
    const double fraction     = (b.position.z() - a.position.z()) / chord.z();
    return (b.position - a.position - fraction * chord).dot(b.fieldDirection);
  }

  double chordResidual(const Point& a, const Point& b, const Point& c) {
    return std::abs(signedChordResidual(a, b, c));
  }

  FieldIntegral fieldIntegral(const Point& a, const Point& b, const Point& c, std::size_t samples,
                              const FieldLookup& fieldAt) {
    const double za     = a.position.z();
    const double zb     = b.position.z();
    const double zc     = c.position.z();
    const double length = zc - za;
    // The Green's function below has a kink at the middle hit, so the samples are placed on each
    // side of it separately: each segment is split into intervals with two Gauss-Legendre points
    // each, exact for an integrand cubic in z. A midpoint rule across the kink is off by a few %.
    const std::size_t intervals = std::max<std::size_t>(1, (samples + 3) / 4);
    const double node           = 0.5 / std::sqrt(3.);
    Acts::Vector3 weightedField = Acts::Vector3::Zero();
    double weights              = 0.;
    double residualPerQOverP    = 0.;
    for (const auto& [lower, upper] : {std::pair{za, zb}, std::pair{zb, zc}}) {
      const double width = (upper - lower) / static_cast<double>(intervals);
      for (std::size_t i = 0; i < intervals; ++i) {
        const double centre = lower + (static_cast<double>(i) + 0.5) * width;
        for (const double z : {centre - node * width, centre + node * width}) {
          // The parabola through the hits, x(z) = sum of hit * Lagrange basis, and its slope
          const std::array<double, 3> basis{(z - zb) * (z - zc) / ((za - zb) * (za - zc)),
                                            (z - za) * (z - zc) / ((zb - za) * (zb - zc)),
                                            (z - za) * (z - zb) / ((zc - za) * (zc - zb))};
          const std::array<double, 3> slope{(2. * z - zb - zc) / ((za - zb) * (za - zc)),
                                            (2. * z - za - zc) / ((zb - za) * (zb - zc)),
                                            (2. * z - za - zb) / ((zc - za) * (zc - zb))};
          const Acts::Vector3 position =
              basis[0] * a.position + basis[1] * b.position + basis[2] * c.position;
          const Acts::Vector3 direction =
              (slope[0] * a.position + slope[1] * b.position + slope[2] * c.position).normalized();
          // Green's function of the middle hit in z: a curvature d2x/dz2 = g at z moves the
          // middle hit off the chord by -g times this weight per unit length in z
          const double weight =
              (z <= zb ? (z - za) * (zc - zb) / length : (zb - za) * (zc - z) / length) * 0.5 *
              width;
          const Acts::Vector3 field = fieldAt(position);
          weightedField += weight * field;
          weights += weight;
          // The curvature per unit q/p is d2r/ds2 = direction x field in Acts units. The weight
          // is a Green's function in z, so convert it to d2r/dz2: divide by t_z^2 and drop the
          // part along the track, which only changes dz/ds. For a 2 GeV track with slopes of up
          // to 0.23 this changes the predicted residual by up to 0.11 mm (8 % of it).
          const Acts::Vector3 bend = direction.cross(field);
          const Acts::Vector3 curvature =
              (bend - direction * (bend.z() / direction.z())) / (direction.z() * direction.z());
          residualPerQOverP -= weight * curvature.dot(b.fieldDirection);
        }
      }
    }
    if (!(weights > 0.) || weightedField.isZero()) {
      return {.field = b.field, .residualPerQOverP = 0.};
    }
    return {.field = weightedField / weights, .residualPerQOverP = residualPerQOverP};
  }

  namespace {
    Acts::FreeVector estimateParameters(const Point& a, const Point& b, const Point& c,
                                        const Acts::Vector3& field) {
#if Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 2)
      return Acts::estimateTrackParamsFromSeed(a.position, a.time, b.position, c.position, field);
#else
      Acts::FreeVector freeParams =
          Acts::estimateTrackParamsFromSeed(a.position, b.position, c.position, field);
      freeParams[Acts::eFreeTime] = a.time;
      return freeParams;
#endif
    }
  } // namespace

  bool inTimeWindow(const B0TripletSeedingConfig& cfg, const Point& first, const Point& later) {
    if (!(cfg.timeWindowSigmas > 0.)) {
      return true;
    }
    // Acts measures time as a length (c = 1), so the lag behind light is in mm
    const double distance = (later.position - first.position).norm();
    const double lag      = later.time - first.time - distance;
    // The flight time of the slowest seeded particle exceeds that of light by distance * (1 / beta
    // - 1), with 1 / beta = E / p at the minimum momentum
    const double slowLag =
        distance * (std::hypot(cfg.minMomentum, cfg.slowParticleMass) / cfg.minMomentum - 1.);
    // A hit without a valid time error cannot be rejected
    const double floor = cfg.timeErrorFloor / edm4eic::unit::ns * Acts::UnitConstants::ns;
    auto error         = [&](const Point& point) {
      return std::isfinite(point.timeError) ? std::max(point.timeError, floor)
                                            : std::numeric_limits<double>::infinity();
    };
    const double tolerance = cfg.timeWindowSigmas * std::hypot(error(first), error(later)) +
                             cfg.timeWindowTolerance / edm4eic::unit::ns * Acts::UnitConstants::ns;
    return lag >= -tolerance && lag <= slowLag + tolerance;
  }

  std::vector<Candidate> findCandidates(const B0TripletSeedingConfig& cfg,
                                        const std::vector<std::vector<Point>>& stations,
                                        const FieldLookup& fieldAt) {
    // Triplets of hits on three stations. In a uniform field the middle hit of a forward track
    // lies close to the chord of the outer two along the field. The quadrupole component of the
    // B0 magnet breaks this off its midplane: the field turns along the track, so a soft track
    // whose triplet skips a station has a residual of up to a few mm. The field integral predicts
    // that term, and it is subtracted before the cut.
    const bool integrate   = cfg.fieldSamples > 0 && fieldAt;
    const double prefilter = cfg.maxResidual + (integrate ? cfg.fieldTermAllowance : 0.);
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < stations.size(); ++i) {
      for (std::size_t j = i + 1; j < stations.size(); ++j) {
        for (std::size_t k = j + 1; k < stations.size(); ++k) {
          for (const auto& a : stations[i]) {
            for (const auto& b : stations[j]) {
              if (b.fieldDirection.isZero() || !inTimeWindow(cfg, a, b)) {
                continue;
              }
              for (const auto& c : stations[k]) {
                if (!inTimeWindow(cfg, b, c) || !inTimeWindow(cfg, a, c)) {
                  continue;
                }
                const double rawResidual = signedChordResidual(a, b, c);
                if (!(std::abs(rawResidual) < prefilter)) {
                  continue;
                }
                const auto integral =
                    integrate ? fieldIntegral(a, b, c, cfg.fieldSamples, fieldAt)
                              : FieldIntegral{.field = b.field, .residualPerQOverP = 0.};
                const double qOverP =
                    estimateParameters(a, b, c, integral.field)[Acts::eFreeQOverP];
                if (!std::isfinite(qOverP) || std::abs(qOverP) * cfg.minMomentum > 1.) {
                  continue;
                }
                const double residual = std::abs(rawResidual - qOverP * integral.residualPerQOverP);
                if (residual < cfg.maxResidual) {
                  candidates.push_back({.residual = residual,
                                        .points   = {&a, &b, &c},
                                        .field    = integral.field,
                                        .qOverP   = qOverP});
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
    auto ranksAfter = [&](const Candidate& first, const Candidate& second) {
      return ranksBefore(second, first);
    };
    std::ranges::make_heap(candidates, ranksAfter);
    auto end          = candidates.end();
    std::size_t seeds = 0;
    while (end != candidates.begin() && seeds < cfg.maxSeeds) {
      std::ranges::pop_heap(candidates.begin(), end, ranksAfter);
      --end;
      if (makeSeed(*end)) {
        ++seeds;
      }
    }
    return static_cast<std::size_t>(std::distance(candidates.begin(), end));
  }

  std::optional<Seed> makeSeed(const B0TripletSeedingConfig& cfg, const Acts::GeometryContext& gctx,
                               const Point& a, const Point& b, const Point& c) {
    return makeSeed(cfg, gctx, a, b, c, b.field);
  }

  std::optional<Seed> makeSeed(const B0TripletSeedingConfig& cfg, const Acts::GeometryContext& gctx,
                               const Point& a, const Point& b, const Point& c,
                               const Acts::Vector3& bField) {
    Acts::FreeVector freeParams = estimateParameters(a, b, c, bField);
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
  const std::array<Named, 8> positive{{
      {.name = "stationGap", .value = m_cfg.stationGap},
      {.name = "maxResidual", .value = m_cfg.maxResidual},
      {.name = "minMomentum", .value = m_cfg.minMomentum},
      {.name = "positionError", .value = m_cfg.positionError},
      {.name = "angleError", .value = m_cfg.angleError},
      {.name = "qOverPRelativeError", .value = m_cfg.qOverPRelativeError},
      {.name = "timeError", .value = m_cfg.timeError},
      {.name = "timeErrorFloor", .value = m_cfg.timeErrorFloor},
  }};
  for (const auto& [name, value] : positive) {
    if (!(value > 0.)) {
      throw std::runtime_error(
          fmt::format("B0TripletSeeding: {} must be positive, got {}", name, value));
    }
  }
  const std::array<Named, 5> nonNegative{{
      {.name = "anchorDistance", .value = m_cfg.anchorDistance},
      {.name = "fieldTermAllowance", .value = m_cfg.fieldTermAllowance},
      {.name = "timeWindowSigmas", .value = m_cfg.timeWindowSigmas},
      {.name = "timeWindowTolerance", .value = m_cfg.timeWindowTolerance},
      {.name = "slowParticleMass", .value = m_cfg.slowParticleMass},
  }};
  for (const auto& [name, value] : nonNegative) {
    if (!(value >= 0.)) {
      throw std::runtime_error(
          fmt::format("B0TripletSeeding: {} must not be negative, got {}", name, value));
    }
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

  // The field along the triplets, from the same cache
  const FieldLookup fieldAt = [&](const Acts::Vector3& position) -> Acts::Vector3 {
    const auto value = field->getField(position, fieldCache);
    return value.ok() ? *value : Acts::Vector3::Zero();
  };
  auto candidates = findCandidates(m_cfg, stations, fieldAt);
  // A track crossing all stations gives one candidate per station triple, and all of them are kept.
  // The track finding does not succeed from every seed of a track, so keeping only the best one
  // per track loses tracks; the ambiguity solver removes the duplicate tracks instead.
  auto& seedsOut            = *track_seeds;
  auto& paramsOut           = *track_params;
  const auto& hitsIn        = *hits;
  const std::size_t dropped = selectSeeds(m_cfg, candidates, [&](const Candidate& candidate) {
    const auto& [a, b, c] = candidate.points;
    const auto seed       = makeSeed(m_cfg, gctx, *a, *b, *c, candidate.field);
    if (seed) {
      writeSeed(*seed, candidate, hitsIn, seedsOut, paramsOut);
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
