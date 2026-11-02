// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0TrackSeeding.h"

#include <Acts/Definitions/Direction.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/SeedContainer2.hpp>
#include <Acts/EventData/SeedProxy2.hpp>
#include <Acts/EventData/SpacePointColumns.hpp>
#include <Acts/EventData/SpacePointContainer2.hpp>
#include <Acts/EventData/SpacePointProxy2.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Seeding2/BroadTripletSeedFilter.hpp>
#include <Acts/Seeding2/DoubletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeedFinder.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/MathHelpers.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Utilities/UnitVectors.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

#include "ActsGeometryProvider.h"
#include "B0TrackSeedingConfig.h"

namespace eicrecon {

std::optional<edm4eic::MutableTrackSeed>
B0TrackSeeding::addSeed(const B0TrackSeedingConfig& cfg, const Acts::Vector3& first,
                        const Acts::Vector3& middle, const Acts::Vector3& last,
                        const Acts::Vector3& field, float time, edm4eic::TrackSeedCollection& seeds,
                        edm4eic::TrackParametersCollection& parameters) {
  const Acts::FreeVector helix = Acts::estimateTrackParamsFromSeed(first, 0., middle, last, field);
  const double qOverP = helix[Acts::eFreeQOverP] * Acts::UnitConstants::GeV / edm4eic::unit::GeV;
  if (!helix.allFinite() || qOverP == 0.) {
    return std::nullopt;
  }

  // CKFTracking rebuilds the perigee along global z from the seed's point, so the seed is at
  // loc0 = loc1 = 0. Forward tracks run nearly parallel to it, so the loc1 and phi errors are
  // scaled up.
  const Acts::Vector3 direction = helix.segment<3>(Acts::eFreeDir0);
  const double theta            = Acts::VectorHelpers::theta(direction);

  edm4eic::Cov6f cov;
  cov(0, 0) = static_cast<float>(Acts::square(cfg.positionError));
  cov(1, 1) = static_cast<float>(Acts::square(cfg.positionError / std::tan(theta)));
  cov(2, 2) = static_cast<float>(Acts::square(cfg.angleError / std::sin(theta)));
  cov(3, 3) = static_cast<float>(Acts::square(cfg.angleError));
  cov(4, 4) = static_cast<float>(Acts::square(cfg.qOverPRelativeError * qOverP));
  cov(5, 5) = static_cast<float>(Acts::square(cfg.timeError));

  auto params = parameters.create();
  params.setType(-1); // seed
  params.setLoc({0.F, 0.F});
  params.setPhi(static_cast<float>(Acts::VectorHelpers::phi(direction)));
  params.setTheta(static_cast<float>(theta));
  params.setQOverP(static_cast<float>(qOverP));
  params.setTime(time);
  params.setCovariance(cov);

  const Acts::Vector3 anchor =
      (first - cfg.anchorDistance * direction) / Acts::UnitConstants::mm * edm4eic::unit::mm;
  auto seed = seeds.create();
  seed.setPerigee({static_cast<float>(anchor.x()), static_cast<float>(anchor.y()),
                   static_cast<float>(anchor.z())});
  seed.setParams(params);
  return seed;
}

void B0TrackSeeding::process(const Input& input, const Output& output) const {
  const auto [hits]            = input;
  auto [trk_seeds, trk_params] = output;

  // Space point i is hit i, with its position in the frame of the current middle hit. All hits are
  // bottom and top candidates: the doublet finders keep those at least deltaRMin from the middle.
  const double toActs = Acts::UnitConstants::mm / edm4eic::unit::mm;
  std::vector<Acts::Vector3> positions;
  Acts::SpacePointContainer2 spacePoints(
      Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  for (const auto& hit : *hits) {
    const auto pos = hit.getPosition();
    const auto cov = hit.getPositionError();
    positions.push_back(toActs * Acts::Vector3(pos.x, pos.y, pos.z));
    auto sp = spacePoints.createSpacePoint();
    sp.varianceZ() =
        static_cast<float>(Acts::square(m_cfg.hitErrorScale * toActs) * std::max(cov.xx, cov.yy));
    sp.varianceR() = sp.varianceZ();
  }
  std::array candidates{spacePoints.range({0, spacePoints.size()}).asConst()};

  auto bottomConfig               = m_cfg.doublets;
  bottomConfig.candidateDirection = Acts::Direction::Backward();

  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  const Acts::BroadTripletSeedFilter filter(m_cfg.filter, filterState, filterCache,
                                            Acts::getDummyLogger());
  Acts::TripletSeeder::Cache seederCache;
  Acts::SeedContainer2 actsSeeds;

  const auto fieldProvider = m_geoSvc->getFieldProvider();
  auto fieldCache          = fieldProvider->makeCache(m_geoSvc->getActsMagneticFieldContext());

  for (std::uint32_t m = 0; m < hits->size(); ++m) {
    const Acts::Vector3 field = fieldProvider->getField(positions[m], fieldCache).value();
    const auto bField         = static_cast<float>(field.norm());
    if (!(bField > 0.F)) {
      continue;
    }

    const Acts::Vector3 fieldDirection = field.normalized();
    const auto [u, v]                  = Acts::createCurvilinearUnitVectors(fieldDirection);
    for (std::uint32_t i = 0; i < hits->size(); ++i) {
      const auto x = static_cast<float>(u.dot(positions[i]));
      const auto y = static_cast<float>(v.dot(positions[i]));
      auto sp      = spacePoints[i];
      sp.xy()      = {x, y};
      sp.zr()      = {static_cast<float>(fieldDirection.dot(positions[i])), std::hypot(x, y)};
    }

    actsSeeds.clear();
    // The middle radius range only applies to hits sorted in radius
    m_seeder.createSeedsFromGroups(
        seederCache, *Acts::DoubletSeedFinder::create({bottomConfig, bField}),
        *Acts::DoubletSeedFinder::create({m_cfg.doublets, bField}),
        *Acts::TripletSeedFinder::create({m_cfg.triplets, bField}), filter, spacePoints, candidates,
        spacePoints.range({m, m + 1}).asConst(), candidates, {0.F, 0.F}, actsSeeds);

    for (const auto& actsSeed : actsSeeds) {
      const auto triplet = actsSeed.spacePointIndices();
      auto seed =
          addSeed(m_cfg, positions[triplet[0]], positions[triplet[1]], positions[triplet[2]], field,
                  (*hits)[triplet[0]].getTime(), *trk_seeds, *trk_params);
      if (seed) {
        seed->setQuality(actsSeed.quality());
        for (const auto i : triplet) {
          seed->addToHits((*hits)[i]);
        }
      }
    }
  }
}

} // namespace eicrecon
