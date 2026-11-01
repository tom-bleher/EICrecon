// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0TrackSeeding.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Direction.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/SeedContainer2.hpp>
#include <Acts/EventData/SpacePointColumns.hpp>
#include <Acts/EventData/SpacePointContainer2.hpp>
#include <Acts/EventData/SpacePointProxy2.hpp>
#include <Acts/EventData/TransformationHelpers.hpp>
#include <Acts/EventData/Types.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Seeding2/BroadTripletSeedFilter.hpp>
#include <Acts/Seeding2/DoubletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeeder.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Result.hpp>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <fmt/format.h>
#include <spdlog/logger.h>
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "extensions/edm4eic/EDM4eicToActs.h"
#include "extensions/spdlog/SpdlogToActs.h"

namespace eicrecon {

namespace {

  /// The Acts finders for one field magnitude
  struct FinderSet {
    std::unique_ptr<Acts::DoubletSeedFinder> bottom;
    std::unique_ptr<Acts::DoubletSeedFinder> top;
    std::unique_ptr<Acts::TripletSeedFinder> triplet;
  };

  /// The time window between the bottom and the top hit of a triplet; the pairs with the middle
  /// hit are checked when the doublet candidates are selected
  class TimeWindowCut final : public Acts::ITripletSeedCuts {
  public:
    TimeWindowCut(const B0TrackSeedingConfig& cfg, const std::vector<B0TrackSeeding::Hit>& hits)
        : m_cfg(cfg), m_hits(hits) {}

    float seedWeight(const Acts::ConstSpacePointProxy2& /*bottom*/,
                     const Acts::ConstSpacePointProxy2& /*middle*/,
                     const Acts::ConstSpacePointProxy2& /*top*/) const override {
      return 0.F;
    }
    bool singleSeedCut(float /*weight*/, const Acts::ConstSpacePointProxy2& bottom,
                       const Acts::ConstSpacePointProxy2& /*middle*/,
                       const Acts::ConstSpacePointProxy2& top) const override {
      return B0TrackSeeding::inTimeWindow(m_cfg, m_hits[bottom.copyFromIndex()],
                                          m_hits[top.copyFromIndex()]);
    }
    void cutPerMiddleSp(std::span<Acts::TripletCandidate2>& /*candidates*/) const override {}

  private:
    const B0TrackSeedingConfig& m_cfg;
    const std::vector<B0TrackSeeding::Hit>& m_hits;
  };

} // namespace

struct B0TrackSeeding::Finders {
  /// Finders for |B| in bins of fieldBinWidth
  std::vector<FinderSet> bins;
  Acts::BroadTripletSeedFilter::Config filter;
  std::unique_ptr<const Acts::Logger> logger;
  Acts::TripletSeeder seeder;
};

Acts::RotationMatrix3 B0TrackSeeding::frameRotation(const Acts::Vector3& field,
                                                    const Acts::Vector3& beam) {
  const Acts::Vector3 z = field.normalized();
  const Acts::Vector3 x = (beam - beam.dot(z) * z).normalized();
  Acts::RotationMatrix3 rotation;
  rotation.row(0) = x.transpose();
  rotation.row(1) = z.cross(x).transpose();
  rotation.row(2) = z.transpose();
  return rotation;
}

bool B0TrackSeeding::inTimeWindow(const B0TrackSeedingConfig& cfg, const Hit& first,
                                  const Hit& later) {
  if (!(cfg.timeWindowSigmas > 0.)) {
    return true;
  }
  // Acts measures time as a length (c = 1), so the lag behind light is in mm. The slowest
  // seeded particle lags by distance * (1 / beta - 1), with 1 / beta = E / p at minMomentum.
  const double distance = (later.position - first.position).norm();
  const double lag      = later.time - first.time - distance;
  const double slowLag =
      distance * (std::hypot(cfg.minMomentum, cfg.slowParticleMass) / cfg.minMomentum - 1.);
  const double floor = cfg.timeErrorFloor / edm4eic::unit::ns * Acts::UnitConstants::ns;
  // A hit without a valid time error cannot be rejected
  auto error = [&](const Hit& hit) {
    return std::isfinite(hit.timeError) ? std::max(hit.timeError, floor)
                                        : std::numeric_limits<double>::infinity();
  };
  const double tolerance = cfg.timeWindowSigmas * std::hypot(error(first), error(later));
  return lag >= -tolerance && lag <= slowLag + tolerance;
}

std::optional<B0TrackSeeding::Seed> B0TrackSeeding::makeSeed(const B0TrackSeedingConfig& cfg,
                                                             const Acts::GeometryContext& gctx,
                                                             const Hit& a, const Hit& b,
                                                             const Hit& c) {
  Acts::FreeVector freeParams =
      Acts::estimateTrackParamsFromSeed(a.position, a.time, b.position, c.position, b.field);
  if (!freeParams.allFinite() || std::abs(freeParams[Acts::eFreeQOverP]) * cfg.minMomentum > 1.) {
    return std::nullopt;
  }

  // Move the parameters estimated at the first hit back along the tangent to the anchor, where
  // CKFTracking starts the track finding. Acts measures time as a length; the shift assumes
  // beta = 1, which is off by 0.01 ns for a 1 GeV proton at the default anchorDistance.
  const Acts::Vector3 direction          = freeParams.segment<3>(Acts::eFreeDir0);
  const Acts::Vector3 anchor             = a.position - cfg.anchorDistance * direction;
  freeParams.segment<3>(Acts::eFreePos0) = anchor;
  freeParams[Acts::eFreeTime] -= cfg.anchorDistance;

  // A TrackSeed carries only a perigee point, and CKFTracking builds a perigee surface (along
  // global z) from it. A forward track is nearly parallel to that line, so the position along it
  // and the azimuth are poorly defined, and their errors are scaled to match.
  const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(anchor);
  const auto bound   = Acts::transformFreeToBoundParameters(freeParams, *perigee, gctx);
  if (!bound.ok() || !bound->allFinite()) {
    return std::nullopt;
  }
  // The three hits resolve the curvature only up to their sagitta error across the field, which
  // limits the q/p of a stiff track: |q/p| = curvature * sin(angle to the field) / |B|
  const Acts::Vector3 fieldDirection = b.field.normalized();
  auto across                        = [&](const Acts::Vector3& v) {
    return (v - v.dot(fieldDirection) * fieldDirection).norm();
  };
  const double d1 = across(b.position - a.position);
  const double d2 = across(c.position - b.position);
  const double sagittaVariance =
      b.variance + (d2 * d2 * a.variance + d1 * d1 * c.variance) / ((d1 + d2) * (d1 + d2));
  const double qOverPResolution = 2. * std::sqrt(sagittaVariance) / (d1 * d2) *
                                  direction.cross(fieldDirection).norm() / b.field.norm();

  const double theta = (*bound)[Acts::eBoundTheta];
  Acts::BoundVector errors;
  errors << cfg.positionError, cfg.positionError / std::abs(std::tan(theta)),
      cfg.angleError / std::sin(theta), cfg.angleError,
      std::max(cfg.qOverPRelativeError * std::abs((*bound)[Acts::eBoundQOverP]), qOverPResolution),
      cfg.timeError / edm4eic::unit::ns * Acts::UnitConstants::ns;
  if (!errors.allFinite() || (errors.array() <= 0.).any()) {
    return std::nullopt;
  }
  return Seed{.anchor = anchor, .parameters = *bound, .errors = errors};
}

void B0TrackSeeding::init() {
  const std::array<std::tuple<std::string_view, double>, 18> positive{{
      {"rMin", m_cfg.rMin},
      {"deltaRMin", m_cfg.deltaRMin},
      {"hitErrorScale", m_cfg.hitErrorScale},
      {"rMax", m_cfg.rMax},
      {"zMax", m_cfg.zMax},
      {"rMaxMiddle", m_cfg.rMaxMiddle},
      {"deltaRMax", m_cfg.deltaRMax},
      {"cotThetaMax", m_cfg.cotThetaMax},
      {"minPt", m_cfg.minPt},
      {"sigmaScattering", m_cfg.sigmaScattering},
      {"radLengthPerSeed", m_cfg.radLengthPerSeed},
      {"fieldBinWidth", m_cfg.fieldBinWidth},
      {"timeErrorFloor", m_cfg.timeErrorFloor},
      {"minMomentum", m_cfg.minMomentum},
      {"positionError", m_cfg.positionError},
      {"angleError", m_cfg.angleError},
      {"qOverPRelativeError", m_cfg.qOverPRelativeError},
      {"timeError", m_cfg.timeError},
  }};
  // Written so that NaN fails too
  for (const auto& [name, value] : positive) {
    if (!(value > 0.)) {
      throw std::runtime_error(
          fmt::format("B0TrackSeeding: {} must be positive, got {}", name, value));
    }
  }
  if (!(m_cfg.timeWindowSigmas >= 0.) || !(m_cfg.anchorDistance >= 0.) ||
      !(m_cfg.slowParticleMass >= 0.) || m_cfg.maxSeeds == 0 || m_cfg.maxSeedsPerMiddle == 0) {
    throw std::runtime_error("B0TrackSeeding: timeWindowSigmas, anchorDistance and "
                             "slowParticleMass must not be negative, maxSeeds and "
                             "maxSeedsPerMiddle must be positive");
  }
  if (!(m_cfg.rMin < m_cfg.rMax) || !(m_cfg.rMinMiddle < m_cfg.rMaxMiddle) ||
      !(m_cfg.deltaRMin < m_cfg.deltaRMax)) {
    throw std::runtime_error("B0TrackSeeding: rMin, rMinMiddle and deltaRMin must be below rMax, "
                             "rMaxMiddle and deltaRMax");
  }
  if (!m_geoSvc || !m_geoSvc->getFieldProvider()) {
    throw std::runtime_error("B0TrackSeeding: no magnetic field provider");
  }

  // x' of the seeding frame along the ion beam, at the crossing angle about global y
  const double angle = m_geoSvc->dd4hepDetector()->constant<double>("CrossingAngle") / dd4hep::rad;
  m_beam             = Acts::Vector3{std::sin(angle), 0., std::cos(angle)};

  // Limits that do not cut on B0 tracks, finite as the Acts finders compute with them
  constexpr float noLimitLength    = 1000. * Acts::UnitConstants::mm;
  constexpr float noLimitCollision = 10000. * Acts::UnitConstants::mm;
  constexpr float noLimitImpact    = 1.e5 * Acts::UnitConstants::mm;

  auto finders = std::make_shared<Finders>();
  // Finders for |B| up to 4 T, each at the centre of its bin
  const auto bins = static_cast<std::size_t>(std::ceil(4. / m_cfg.fieldBinWidth));
  for (std::size_t i = 0; i < bins; ++i) {
    const auto field = static_cast<float>((static_cast<double>(i) + 0.5) * m_cfg.fieldBinWidth *
                                          Acts::UnitConstants::T);
    // No interaction point cut, and collision region, impact and slope limits far beyond the B0
    // acceptance: tracks may come from decays anywhere upstream
    Acts::DoubletSeedFinder::Config bottom;
    bottom.spacePointsSortedByRadius    = false;
    bottom.candidateDirection           = Acts::Direction::Backward();
    bottom.deltaRMin                    = m_cfg.deltaRMin;
    bottom.deltaRMax                    = m_cfg.deltaRMax;
    bottom.deltaZMin                    = -noLimitLength;
    bottom.deltaZMax                    = noLimitLength;
    bottom.interactionPointCut          = false;
    bottom.impactMax                    = noLimitImpact;
    bottom.collisionRegionMin           = -noLimitCollision;
    bottom.collisionRegionMax           = noLimitCollision;
    bottom.cotThetaMax                  = m_cfg.cotThetaMax;
    bottom.minPt                        = m_cfg.minPt;
    Acts::DoubletSeedFinder::Config top = bottom;
    top.candidateDirection              = Acts::Direction::Forward();
    Acts::TripletSeedFinder::Config triplet;
    triplet.useStripInfo     = false;
    triplet.sortedByCotTheta = true;
    triplet.minPt            = m_cfg.minPt;
    triplet.sigmaScattering  = m_cfg.sigmaScattering;
    triplet.radLengthPerSeed = m_cfg.radLengthPerSeed;
    triplet.impactMax        = noLimitImpact;
    finders->bins.push_back(
        {.bottom =
             Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(bottom, field)),
         .top = Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(top, field)),
         .triplet = Acts::TripletSeedFinder::create(
             Acts::TripletSeedFinder::DerivedConfig(triplet, field))});
  }
  // Seeds are not weighted by their impact or origin along the beam
  auto& filter                     = finders->filter;
  filter.deltaInvHelixDiameter     = 0.00003 / Acts::UnitConstants::mm;
  filter.deltaRMin                 = m_cfg.deltaRMin;
  filter.compatSeedWeight          = 200.;
  filter.compatSeedLimit           = 2;
  filter.impactWeightFactor        = 0.;
  filter.zOriginWeightFactor       = 0.;
  filter.maxSeedsPerSpM            = m_cfg.maxSeedsPerMiddle;
  filter.maxSeedsPerSpMConf        = 10;
  filter.maxQualitySeedsPerSpMConf = 10;
  filter.seedConfirmation          = false;

  auto logger = spdlog::default_logger()->clone(std::string(name()));
  logger->set_level(static_cast<spdlog::level::level_enum>(level()));
  finders->logger = getSpdlogLogger(std::string(name()), logger);
  m_finders       = finders;
}

void B0TrackSeeding::process(const Input& input, const Output& output) const {
  const auto [hitCollection]       = input;
  auto [track_seeds, track_params] = output;

  const auto& gctx = m_geoSvc->getActsGeometryContext();
  const auto field = m_geoSvc->getFieldProvider();
  auto fieldCache  = field->makeCache(m_geoSvc->getActsMagneticFieldContext());

  // The hits in Acts units with the field at each, sorted so that the seeds do not depend on the
  // order of the input. The variance is the larger in-plane one of the sensor.
  std::vector<Hit> hits;
  for (std::size_t i = 0; i < hitCollection->size(); ++i) {
    const auto& hit = (*hitCollection)[i];
    const auto pos  = hit.getPosition();
    if (!std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z) ||
        !std::isfinite(hit.getTime())) {
      continue;
    }
    const Acts::Vector3 position{pos.x / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                 pos.y / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                 pos.z / edm4eic::unit::mm * Acts::UnitConstants::mm};
    const auto value = field->getField(position, fieldCache);
    const auto cov   = hit.getPositionError();
    hits.push_back({.position  = position,
                    .field     = value.ok() ? *value : Acts::Vector3::Zero(),
                    .time      = hit.getTime() / edm4eic::unit::ns * Acts::UnitConstants::ns,
                    .timeError = hit.getTimeError() / edm4eic::unit::ns * Acts::UnitConstants::ns,
                    .variance = std::max(cov.xx, cov.yy) / (edm4eic::unit::mm * edm4eic::unit::mm) *
                                (Acts::UnitConstants::mm * Acts::UnitConstants::mm),
                    .index    = i});
  }
  std::ranges::sort(hits, [](const Hit& lhs, const Hit& rhs) {
    return std::tuple{lhs.position.z(), lhs.position.x(), lhs.position.y(), lhs.index} <
           std::tuple{rhs.position.z(), rhs.position.x(), rhs.position.y(), rhs.index};
  });

  const auto& finders         = *m_finders;
  auto filterConfig           = finders.filter;
  filterConfig.experimentCuts = std::make_shared<TimeWindowCut>(m_cfg, hits);
  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter filter(filterConfig, filterState, filterCache, *finders.logger);
  Acts::TripletSeeder::Cache seederCache;

  const auto columns = Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
                       Acts::SpacePointColumns::Phi | Acts::SpacePointColumns::VarianceZ |
                       Acts::SpacePointColumns::VarianceR | Acts::SpacePointColumns::CopyFromIndex;
  Acts::SpacePointContainer2 spacePoints(columns);
  std::vector<Acts::SpacePointIndex2> bottoms;
  std::vector<Acts::SpacePointIndex2> tops;

  // Hit errors for the Acts finders, scaled to floor their slope cut (see B0TrackSeedingConfig)
  const double hitVariance = m_cfg.hitErrorScale * m_cfg.hitErrorScale;

  struct Found {
    float weight;
    std::array<std::uint32_t, 3> hits;
  };
  std::vector<Found> found;

  for (std::uint32_t m = 0; m < hits.size(); ++m) {
    const Hit& middle   = hits[m];
    const double along  = middle.position.dot(m_beam);
    const double bField = middle.field.norm() / Acts::UnitConstants::T;
    if (along < m_cfg.rMinMiddle || along > m_cfg.rMaxMiddle || !(bField > 0.)) {
      continue;
    }

    // The hits in the frame along the field at the middle hit
    const auto rotation = frameRotation(middle.field, m_beam);
    spacePoints.clear();
    spacePoints.reserve(hits.size());
    std::optional<Acts::SpacePointIndex2> middleIndex;
    for (std::uint32_t i = 0; i < hits.size(); ++i) {
      const Acts::Vector3 primed = rotation * hits[i].position;
      const double r             = std::hypot(primed.x(), primed.y());
      if (r < m_cfg.rMin || r > m_cfg.rMax || std::abs(primed.z()) > m_cfg.zMax) {
        continue;
      }
      if (i == m) {
        middleIndex = spacePoints.size();
      }
      auto sp            = spacePoints.createSpacePoint();
      sp.xy()            = {static_cast<float>(primed.x()), static_cast<float>(primed.y())};
      sp.zr()            = {static_cast<float>(primed.z()), static_cast<float>(r)};
      sp.phi()           = static_cast<float>(std::atan2(primed.y(), primed.x()));
      sp.varianceZ()     = static_cast<float>(hitVariance * hits[i].variance);
      sp.varianceR()     = static_cast<float>(hitVariance * hits[i].variance);
      sp.copyFromIndex() = i;
    }
    if (!middleIndex) {
      continue;
    }
    const Acts::SpacePointContainer2& points = spacePoints;
    const auto spM                           = points.at(*middleIndex);

    // Doublet candidates before and after the middle hit along the beam, in time with it
    bottoms.clear();
    tops.clear();
    for (Acts::SpacePointIndex2 o = 0; o < spacePoints.size(); ++o) {
      if (o == *middleIndex) {
        continue;
      }
      const Hit& other = hits[points.at(o).copyFromIndex()];
      if (points.at(o).zr()[1] < spM.zr()[1]) {
        if (inTimeWindow(m_cfg, other, middle)) {
          bottoms.push_back(o);
        }
      } else if (inTimeWindow(m_cfg, middle, other)) {
        tops.push_back(o);
      }
    }
    if (bottoms.empty() || tops.empty()) {
      continue;
    }

    const auto& set = finders.bins[std::min(static_cast<std::size_t>(bField / m_cfg.fieldBinWidth),
                                            finders.bins.size() - 1)];
    auto bottomSps  = points.subset(std::span<const Acts::SpacePointIndex2>(bottoms));
    auto topSps     = points.subset(std::span<const Acts::SpacePointIndex2>(tops));
    Acts::SeedContainer2 seeds;
    seeds.assignSpacePointContainer(points);
    finders.seeder.createSeedsFromGroup(seederCache, *set.bottom, *set.top, *set.triplet, filter,
                                        points, bottomSps, spM, topSps, seeds);
    for (const auto& seed : seeds) {
      const auto& indices = seed.spacePointIndices();
      found.push_back(
          {.weight = seed.quality(),
           .hits   = {points.at(indices[0]).copyFromIndex(), points.at(indices[1]).copyFromIndex(),
                      points.at(indices[2]).copyFromIndex()}});
    }
  }

  // Highest weight first, ties broken by the hits, so that the result is deterministic
  std::ranges::sort(found, [&](const Found& lhs, const Found& rhs) {
    return std::tuple{-lhs.weight, hits[lhs.hits[0]].index, hits[lhs.hits[1]].index,
                      hits[lhs.hits[2]].index} < std::tuple{-rhs.weight, hits[rhs.hits[0]].index,
                                                            hits[rhs.hits[1]].index,
                                                            hits[rhs.hits[2]].index};
  });

  std::size_t written = 0;
  for (const auto& f : found) {
    if (written == m_cfg.maxSeeds) {
      // Logged the 1st, 2nd, 4th, ... time, as busy samples can reach it in many events
      const std::size_t truncated = m_truncatedEvents.fetch_add(1) + 1;
      if (std::has_single_bit(truncated)) {
        warning("Reached maxSeeds={} with {} seed candidates in this event ({} events so far)",
                m_cfg.maxSeeds, found.size(), truncated);
      }
      break;
    }
    const Hit& a    = hits[f.hits[0]];
    const Hit& b    = hits[f.hits[1]];
    const Hit& c    = hits[f.hits[2]];
    const auto seed = makeSeed(m_cfg, gctx, a, b, c);
    if (!seed) {
      continue;
    }

    const auto& par = seed->parameters;
    auto params     = track_params->create();
    params.setType(-1); // seed
    params.setLoc(
        {static_cast<float>(par[Acts::eBoundLoc0] / Acts::UnitConstants::mm * edm4eic::unit::mm),
         static_cast<float>(par[Acts::eBoundLoc1] / Acts::UnitConstants::mm * edm4eic::unit::mm)});
    params.setPhi(static_cast<float>(par[Acts::eBoundPhi] / Acts::UnitConstants::rad));
    params.setTheta(static_cast<float>(par[Acts::eBoundTheta] / Acts::UnitConstants::rad));
    params.setQOverP(static_cast<float>(par[Acts::eBoundQOverP] * Acts::UnitConstants::GeV /
                                        edm4eic::unit::GeV));
    params.setTime(
        static_cast<float>(par[Acts::eBoundTime] / Acts::UnitConstants::ns * edm4eic::unit::ns));
    edm4eic::Cov6f cov;
    for (std::size_t i = 0; const auto& [ia, unit] : edm4eic_indexed_units) {
      cov(i, i) = static_cast<float>(seed->errors[ia] * seed->errors[ia] / (unit * unit));
      ++i;
    }
    params.setCovariance(cov);

    auto trackSeed = track_seeds->create();
    trackSeed.setPerigee(
        {static_cast<float>(seed->anchor.x() / Acts::UnitConstants::mm * edm4eic::unit::mm),
         static_cast<float>(seed->anchor.y() / Acts::UnitConstants::mm * edm4eic::unit::mm),
         static_cast<float>(seed->anchor.z() / Acts::UnitConstants::mm * edm4eic::unit::mm)});
    trackSeed.setQuality(f.weight);
    trackSeed.setParams(params);
    for (const Hit* hit : {&a, &b, &c}) {
      trackSeed.addToHits((*hitCollection)[hit->index]);
    }
    ++written;
  }

  debug("{} hits, {} seed candidates, {} seeds", hits.size(), found.size(), written);
}

} // namespace eicrecon
