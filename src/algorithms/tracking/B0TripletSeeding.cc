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
#include <edm4eic/Cov6f.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstddef>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

#include "ActsGeometryProvider.h"
#include "B0TripletSeedingConfig.h"
#include "extensions/edm4eic/EDM4eicToActs.h"

namespace eicrecon {

void B0TripletSeeding::process(const Input& input, const Output& output) const {
  const auto [hits]                = input;
  auto [track_seeds, track_params] = output;

  if (hits->empty()) {
    return;
  }

  const auto& gctx = m_geoSvc->getActsGeometryContext();
  const auto& mctx = m_geoSvc->getActsMagneticFieldContext();
  const auto field = m_geoSvc->getFieldProvider();
  auto fieldCache  = field->makeCache(mctx);

  struct Point {
    Acts::Vector3 position;
    double time;
    std::size_t index;
  };

  // Group hits into stations: sort by z and start a new station at every gap.
  // Hits on the front and back sensors of a station stay separate points.
  std::vector<std::size_t> order(hits->size());
  std::iota(order.begin(), order.end(), 0);
  auto sortKey = [&](std::size_t i) {
    const auto& pos = (*hits)[i].getPosition();
    return std::tuple{pos.z, pos.x, pos.y, (*hits)[i].getCellID()};
  };
  std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return sortKey(a) < sortKey(b); });

  std::vector<std::vector<Point>> stations;
  double lastZ = -std::numeric_limits<double>::infinity();
  for (const auto i : order) {
    const auto& pos = (*hits)[i].getPosition();
    const Acts::Vector3 position{pos.x / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                 pos.y / edm4eic::unit::mm * Acts::UnitConstants::mm,
                                 pos.z / edm4eic::unit::mm * Acts::UnitConstants::mm};
    const double time = (*hits)[i].getTime() / edm4eic::unit::ns * Acts::UnitConstants::ns;
    if (position.z() - lastZ > m_cfg.stationGap) {
      stations.emplace_back();
    }
    lastZ = position.z();
    stations.back().push_back({position, time, i});
  }

  // Triplets of hits on three stations. A helix does not bend along the field,
  // so the middle hit must lie on the chord of the outer two in that direction.
  struct Candidate {
    double residual;
    std::array<const Point*, 3> points;
    Acts::Vector3 field;
  };
  std::vector<Candidate> candidates;
  for (std::size_t i = 0; i < stations.size(); ++i) {
    for (std::size_t j = i + 1; j < stations.size(); ++j) {
      for (const auto& b : stations[j]) {
        const auto bField = field->getField(b.position, fieldCache);
        if (!bField.ok() || bField->norm() == 0.) {
          continue;
        }
        const Acts::Vector3 fieldDirection = bField->normalized();
        for (std::size_t k = j + 1; k < stations.size(); ++k) {
          for (const auto& a : stations[i]) {
            for (const auto& c : stations[k]) {
              const Acts::Vector3 chord = c.position - a.position;
              const double fraction     = (b.position.z() - a.position.z()) / chord.z();
              const double residual =
                  std::abs((b.position - a.position - fraction * chord).dot(fieldDirection));
              if (residual < m_cfg.maxResidual) {
                candidates.push_back({residual, {&a, &b, &c}, *bField});
              }
            }
          }
        }
      }
    }
  }
  std::ranges::sort(candidates, [](const auto& lhs, const auto& rhs) {
    return std::tie(lhs.residual, lhs.points[0]->index, lhs.points[1]->index,
                    lhs.points[2]->index) <
           std::tie(rhs.residual, rhs.points[0]->index, rhs.points[1]->index, rhs.points[2]->index);
  });

  for (const auto& [residual, points, bField] : candidates) {
    if (track_seeds->size() >= m_cfg.maxSeeds) {
      warning("Keeping {} of {} seed candidates", m_cfg.maxSeeds, candidates.size());
      break;
    }
    const auto& [a, b, c] = points;

#if Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 2)
    Acts::FreeVector freeParams =
        Acts::estimateTrackParamsFromSeed(a->position, a->time, b->position, c->position, bField);
#else
    Acts::FreeVector freeParams =
        Acts::estimateTrackParamsFromSeed(a->position, b->position, c->position, bField);
    freeParams[Acts::eFreeTime] = a->time;
#endif
    if (!freeParams.allFinite() ||
        std::abs(freeParams[Acts::eFreeQOverP]) * m_cfg.minMomentum > 1.) {
      continue;
    }

    // Express the seed on a perigee surface just upstream of its first hit,
    // where it was measured, rather than transporting it to the origin
    const Acts::Vector3 anchor =
        a->position - m_cfg.anchorDistance * freeParams.segment<3>(Acts::eFreeDir0);
    freeParams.segment<3>(Acts::eFreePos0) = anchor;
    freeParams[Acts::eFreeTime] -= m_cfg.anchorDistance;
    const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(anchor);
    const auto local   = Acts::transformFreeToBoundParameters(freeParams, *perigee, gctx);
    if (!local.ok() || !local->allFinite()) {
      continue;
    }
    const auto& parameter = *local;

    // The azimuth and the position along the perigee line of a forward track
    // are poorly defined, so scale their errors
    const double theta = parameter[Acts::eBoundTheta];
    Acts::BoundVector errors;
    errors << m_cfg.positionError, m_cfg.positionError / std::tan(theta),
        m_cfg.angleError / std::sin(theta), m_cfg.angleError,
        m_cfg.qOverPRelativeError * std::abs(parameter[Acts::eBoundQOverP]), m_cfg.timeError;

    auto pars = track_params->create();
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
    edm4eic::Cov6f cov;
    for (std::size_t i = 0; const auto& [ia, x] : edm4eic_indexed_units) {
      for (std::size_t j = 0; const auto& [ib, y] : edm4eic_indexed_units) {
        cov(i, j) = (ia == ib ? errors[ia] * errors[ia] : 0.) / x / y;
        ++j;
      }
      ++i;
    }
    pars.setCovariance(cov);

    auto seed = track_seeds->create();
    seed.setPerigee({static_cast<float>(anchor.x() / Acts::UnitConstants::mm * edm4eic::unit::mm),
                     static_cast<float>(anchor.y() / Acts::UnitConstants::mm * edm4eic::unit::mm),
                     static_cast<float>(anchor.z() / Acts::UnitConstants::mm * edm4eic::unit::mm)});
    seed.setQuality(static_cast<float>(-residual / Acts::UnitConstants::mm));
    seed.setParams(pars);
    for (const auto* point : points) {
      seed.addToHits((*hits)[point->index]);
    }
  }

  debug("{} stations, {} triplet candidates, {} seeds", stations.size(), candidates.size(),
        track_seeds->size());
}

} // namespace eicrecon
