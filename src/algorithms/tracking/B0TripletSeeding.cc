// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0TripletSeeding.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Direction.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/GenericBoundTrackParameters.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#endif
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/TransformationHelpers.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Propagator/ActorList.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/PropagatorResult.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
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
#include <optional>
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

  // CKFTracking expects parameters on a perigee at the origin. The seed is
  // transported there through the field; this is a change of reference
  // surface, not a vertex constraint.
  const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  Acts::Propagator<Acts::EigenStepper<>> propagator{Acts::EigenStepper<>{field}};
  decltype(propagator)::Options<> options(gctx, mctx);
  options.direction = Acts::Direction::Backward();

  for (const auto& [residual, points, bField] : candidates) {
    if (track_seeds->size() >= m_cfg.maxSeeds) {
      warning("Keeping {} of {} seed candidates", m_cfg.maxSeeds, candidates.size());
      break;
    }
    const auto& [a, b, c] = points;

    auto transform          = Acts::Transform3::Identity();
    transform.translation() = a->position;
    const auto plane        = Acts::Surface::makeShared<Acts::PlaneSurface>(transform);
#if Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 2)
    const Acts::FreeVector freeParams =
        Acts::estimateTrackParamsFromSeed(a->position, a->time, b->position, c->position, bField);
#else
    Acts::FreeVector freeParams =
        Acts::estimateTrackParamsFromSeed(a->position, b->position, c->position, bField);
    freeParams[Acts::eFreeTime] = a->time;
#endif
    const auto local = Acts::transformFreeToBoundParameters(freeParams, *plane, gctx);
    if (!local.ok() || !local->allFinite() ||
        std::abs((*local)[Acts::eBoundQOverP]) * m_cfg.minMomentum > 1.) {
      continue;
    }

    // The azimuth of a forward track is poorly defined, so scale its error
    const double theta  = (*local)[Acts::eBoundTheta];
    const double qOverP = (*local)[Acts::eBoundQOverP];
    Acts::BoundVector errors;
    errors << m_cfg.positionError, m_cfg.positionError, m_cfg.angleError / std::sin(theta),
        m_cfg.angleError, m_cfg.qOverPRelativeError * std::abs(qOverP), m_cfg.timeError;
#if Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 1)
    const Acts::BoundMatrix covariance = errors.cwiseAbs2().asDiagonal();
#else
    const Acts::BoundSquareMatrix covariance = errors.cwiseAbs2().asDiagonal();
#endif

    const Acts::BoundTrackParameters start(plane, *local, covariance,
                                           Acts::ParticleHypothesis::pion());
    const auto result = propagator.propagate(start, *perigee, options);
    if (!result.ok() || !result->endParameters || !result->endParameters->covariance()) {
      debug("Seed transport to the perigee failed from z = {} mm",
            a->position.z() / Acts::UnitConstants::mm);
      continue;
    }
    const auto& parameter      = result->endParameters->parameters();
    const auto& seedCovariance = *result->endParameters->covariance();

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
        cov(i, j) = seedCovariance(ia, ib) / x / y;
        ++j;
      }
      ++i;
    }
    pars.setCovariance(cov);

    auto seed = track_seeds->create();
    seed.setPerigee({0.F, 0.F, 0.F});
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
