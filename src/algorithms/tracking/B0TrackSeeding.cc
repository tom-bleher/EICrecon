// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0TrackSeeding.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Direction.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/SeedContainer2.hpp>
#include <Acts/EventData/SeedProxy2.hpp>
#include <Acts/EventData/SpacePointColumns.hpp>
#include <Acts/EventData/SpacePointContainer2.hpp>
#include <Acts/EventData/SpacePointProxy2.hpp>
#include <Acts/EventData/Types.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>
#include <DD4hep/Detector.h>
#include <Evaluator/DD4hepUnits.h>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <spdlog/common.h>
#include <spdlog/logger.h>
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "ActsGeometryProvider.h"
#include "B0TrackSeedingConfig.h"
#include "extensions/edm4eic/EDM4eicToActs.h"
#include "extensions/spdlog/SpdlogToActs.h"

namespace eicrecon {

namespace {

  /// The three hits of a triplet are in time with each other
  class TimeWindowCut final : public Acts::ITripletSeedCuts {
  public:
    TimeWindowCut(const B0TrackSeedingConfig& cfg, std::span<const B0TrackSeeding::Hit> hits)
        : m_cfg(cfg), m_hits(hits) {}

    float seedWeight(const Acts::ConstSpacePointProxy2& /*bottom*/,
                     const Acts::ConstSpacePointProxy2& /*middle*/,
                     const Acts::ConstSpacePointProxy2& /*top*/) const override {
      return 0.F;
    }
    bool singleSeedCut(float /*weight*/, const Acts::ConstSpacePointProxy2& bottom,
                       const Acts::ConstSpacePointProxy2& middle,
                       const Acts::ConstSpacePointProxy2& top) const override {
      const auto& b = m_hits[bottom.index()];
      const auto& m = m_hits[middle.index()];
      const auto& t = m_hits[top.index()];
      return B0TrackSeeding::inTimeWindow(m_cfg, b, m) &&
             B0TrackSeeding::inTimeWindow(m_cfg, m, t) && B0TrackSeeding::inTimeWindow(m_cfg, b, t);
    }
    void cutPerMiddleSp(std::span<Acts::TripletCandidate2>& /*candidates*/) const override {}

  private:
    B0TrackSeedingConfig m_cfg;
    std::span<const B0TrackSeeding::Hit> m_hits;
  };

} // namespace

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
  // Acts time is a length (c = 1). The slowest seeded particle lags light by
  // distance * (1 / beta - 1), with 1 / beta = E / p at minMomentum.
  const double distance = (later.position - first.position).norm();
  const double lag      = later.time - first.time - distance;
  const double slowLag =
      distance * (std::hypot(cfg.minMomentum, cfg.slowParticleMass) / cfg.minMomentum - 1.);
  const double floor = cfg.timeErrorFloor / edm4eic::unit::ns * Acts::UnitConstants::ns;
  // A hit without a valid time error is never rejected
  auto error = [&](const Hit& hit) {
    return std::isfinite(hit.timeError) ? std::max(hit.timeError, floor)
                                        : std::numeric_limits<double>::infinity();
  };
  const double tolerance = cfg.timeWindowSigmas * std::hypot(error(first), error(later));
  return lag >= -tolerance && lag <= slowLag + tolerance;
}

std::optional<B0TrackSeeding::Seed> B0TrackSeeding::makeSeed(const B0TrackSeedingConfig& cfg,
                                                             const Hit& a, const Hit& b,
                                                             const Hit& c) {
  Acts::FreeVector freeParams =
      Acts::estimateTrackParamsFromSeed(a.position, a.time, b.position, c.position, b.field);
  if (!freeParams.allFinite() || std::abs(freeParams[Acts::eFreeQOverP]) * cfg.minMomentum > 1.) {
    return std::nullopt;
  }

  // Move back along the tangent to the anchor, where CKFTracking starts. The time shift assumes
  // beta = 1: 0.01 ns off for a 1 GeV proton.
  const Acts::Vector3 direction          = freeParams.segment<3>(Acts::eFreeDir0);
  const Acts::Vector3 anchor             = a.position - cfg.anchorDistance * direction;
  freeParams.segment<3>(Acts::eFreePos0) = anchor;
  freeParams[Acts::eFreeTime] -= cfg.anchorDistance;

  // CKFTracking rebuilds the perigee along global z from the seed's point, so the seed is at
  // loc0 = loc1 = 0. Forward tracks run nearly parallel to it, so the loc1 and phi errors are
  // scaled up.
  const double theta             = Acts::VectorHelpers::theta(direction);
  Acts::BoundVector parameters   = Acts::BoundVector::Zero();
  parameters[Acts::eBoundPhi]    = Acts::VectorHelpers::phi(direction);
  parameters[Acts::eBoundTheta]  = theta;
  parameters[Acts::eBoundQOverP] = freeParams[Acts::eFreeQOverP];
  parameters[Acts::eBoundTime]   = freeParams[Acts::eFreeTime];

  Acts::BoundVector errors;
  errors << cfg.positionError, cfg.positionError / std::abs(std::tan(theta)),
      cfg.angleError / std::sin(theta), cfg.angleError,
      cfg.qOverPRelativeError * std::abs(parameters[Acts::eBoundQOverP]),
      cfg.timeError / edm4eic::unit::ns * Acts::UnitConstants::ns;
  if (!errors.allFinite() || (errors.array() <= 0.).any()) {
    return std::nullopt;
  }
  return Seed{.anchor = anchor, .parameters = parameters, .errors = errors};
}

void B0TrackSeeding::init() {
  const double angle = m_geoSvc->dd4hepDetector()->constant<double>("CrossingAngle") / dd4hep::rad;
  m_beam             = Acts::Vector3{std::sin(angle), 0., std::cos(angle)};

  // Limits far beyond the B0 acceptance
  constexpr float noLimitCollision = 10000. * Acts::UnitConstants::mm;
  constexpr float noLimitImpact    = 1.e5 * Acts::UnitConstants::mm;

  // Tracks may come from decays upstream: the collision-region and impact cuts are opened
  auto& bottom              = m_bottomConfig;
  bottom.candidateDirection = Acts::Direction::Backward();
  bottom.deltaRMin          = m_cfg.deltaRMin;
  bottom.deltaRMax          = m_cfg.deltaRMax;
  bottom.collisionRegionMin = -noLimitCollision;
  bottom.collisionRegionMax = noLimitCollision;
  bottom.cotThetaMax        = m_cfg.cotThetaMax;

  m_topConfig                    = bottom;
  m_topConfig.candidateDirection = Acts::Direction::Forward();

  auto& triplet            = m_tripletConfig;
  triplet.minPt            = m_cfg.minPt;
  triplet.sigmaScattering  = m_cfg.sigmaScattering;
  triplet.radLengthPerSeed = m_cfg.radLengthPerSeed;
  triplet.impactMax        = noLimitImpact;

  auto& filter              = m_filterConfig;
  filter.deltaRMin          = m_cfg.deltaRMin;
  filter.impactWeightFactor = 0.;
  // Candidates kept per middle hit, used even without seed confirmation. Of these, at most
  // maxSeedsPerSpM + 1 (6) become seeds.
  filter.maxSeedsPerSpMConf        = 10;
  filter.maxQualitySeedsPerSpMConf = 10;

  auto logger = spdlog::default_logger()->clone(std::string(name()));
  logger->set_level(static_cast<spdlog::level::level_enum>(level()));
  m_filterLogger = getSpdlogLogger(std::string(name()), logger);
  m_seeder =
      std::make_unique<const Acts::TripletSeeder>(getSpdlogLogger(std::string(name()), logger));
}

void B0TrackSeeding::process(const Input& input, const Output& output) const {
  const auto [trk_hits]        = input;
  auto [trk_seeds, trk_params] = output;
  if (trk_hits->empty()) {
    return;
  }

  const auto field = m_geoSvc->getFieldProvider();
  auto fieldCache  = field->makeCache(m_geoSvc->getActsMagneticFieldContext());

  // Hits sorted by position, so that the seeds do not depend on the input order
  std::vector<Hit> hits;
  hits.reserve(trk_hits->size());
  for (std::size_t i = 0; i < trk_hits->size(); ++i) {
    const auto& hit = (*trk_hits)[i];
    const auto pos  = hit.getPosition();
    const Acts::Vector3 position =
        Acts::Vector3(pos.x, pos.y, pos.z) / edm4eic::unit::mm * Acts::UnitConstants::mm;
    const auto cov = hit.getPositionError();
    hits.push_back({.position  = position,
                    .field     = field->getField(position, fieldCache).value(),
                    .time      = hit.getTime() / edm4eic::unit::ns * Acts::UnitConstants::ns,
                    .timeError = hit.getTimeError() / edm4eic::unit::ns * Acts::UnitConstants::ns,
                    .variance = std::max(cov.xx, cov.yy) / (edm4eic::unit::mm * edm4eic::unit::mm) *
                                (Acts::UnitConstants::mm * Acts::UnitConstants::mm),
                    .index    = i});
  }
  std::ranges::sort(hits, {}, [](const Hit& hit) {
    return std::tuple{hit.position.z(), hit.position.x(), hit.position.y(), hit.index};
  });

  auto filterConfig           = m_filterConfig;
  filterConfig.experimentCuts = std::make_shared<TimeWindowCut>(m_cfg, hits);
  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter filter(filterConfig, filterState, filterCache, *m_filterLogger);
  Acts::TripletSeeder::Cache seederCache;

  // Every hit is a space point, in the same order: space point i is hit i, with its position in the
  // frame of the current middle hit. All are bottom and top candidates; the doublet finders keep
  // those deltaRMin to deltaRMax before or after the middle hit in r'.
  Acts::SpacePointContainer2 spacePoints(
      Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  spacePoints.reserve(hits.size());
  const double varianceScale = m_cfg.hitErrorScale * m_cfg.hitErrorScale;
  for (const Hit& hit : hits) {
    auto sp        = spacePoints.createSpacePoint();
    sp.varianceZ() = static_cast<float>(varianceScale * hit.variance);
    sp.varianceR() = sp.varianceZ();
  }
  std::vector<Acts::SpacePointIndex2> indices(hits.size());
  std::iota(indices.begin(), indices.end(), Acts::SpacePointIndex2{0});
  auto candidates = spacePoints.subset(indices).asConst();

  Acts::SeedContainer2 actsSeeds;

  for (std::uint32_t m = 0; m < hits.size(); ++m) {
    const Hit& middle = hits[m];
    const auto bField = static_cast<float>(middle.field.norm());
    // No field outside the field map
    if (!(bField > 0.F)) {
      continue;
    }

    const auto rotation = frameRotation(middle.field, m_beam);
    for (std::uint32_t i = 0; i < hits.size(); ++i) {
      const Acts::Vector3 primed = rotation * hits[i].position;
      auto sp                    = spacePoints[i];
      sp.xy()                    = {static_cast<float>(primed.x()), static_cast<float>(primed.y())};
      sp.zr()                    = {static_cast<float>(primed.z()),
                                    static_cast<float>(std::hypot(primed.x(), primed.y()))};
    }

    // The Acts finders for the field magnitude at the middle hit, along z' in its frame
    const auto bottomFinder  = Acts::DoubletSeedFinder::create({m_bottomConfig, bField});
    const auto topFinder     = Acts::DoubletSeedFinder::create({m_topConfig, bField});
    const auto tripletFinder = Acts::TripletSeedFinder::create({m_tripletConfig, bField});
    m_seeder->createSeedsFromGroup(seederCache, *bottomFinder, *topFinder, *tripletFinder, filter,
                                   spacePoints, candidates, spacePoints.at(m).asConst(), candidates,
                                   actsSeeds);
  }

  const auto edmLength = [](double length) {
    return static_cast<float>(length / Acts::UnitConstants::mm * edm4eic::unit::mm);
  };
  for (const auto& actsSeed : actsSeeds) {
    const auto triplet = actsSeed.spacePointIndices();
    const auto seed    = makeSeed(m_cfg, hits[triplet[0]], hits[triplet[1]], hits[triplet[2]]);
    if (!seed) {
      continue;
    }

    const auto& par = seed->parameters;
    auto params     = trk_params->create();
    params.setType(-1); // seed
    params.setLoc({edmLength(par[Acts::eBoundLoc0]), edmLength(par[Acts::eBoundLoc1])});
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

    auto trkSeed = trk_seeds->create();
    trkSeed.setPerigee(
        {edmLength(seed->anchor.x()), edmLength(seed->anchor.y()), edmLength(seed->anchor.z())});
    trkSeed.setQuality(actsSeed.quality());
    trkSeed.setParams(params);
    for (const auto i : triplet) {
      trkSeed.addToHits((*trk_hits)[hits[i].index]);
    }
  }

  debug("{} hits, {} seed candidates, {} seeds", hits.size(), actsSeeds.size(), trk_seeds->size());
}

} // namespace eicrecon
