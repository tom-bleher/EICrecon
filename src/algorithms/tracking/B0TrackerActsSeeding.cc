// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0TrackerActsSeeding.h"

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Evaluator/DD4hepUnits.h>
#include <Eigen/Cholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/tracking/B0TrackerStubSeeder.h"

namespace eicrecon {
namespace b0acts {

  std::optional<Parameters> estimateParameters(const std::vector<Measurement>& measurements,
                                               const Acts::GeometryContext& gctx,
                                               const Acts::Vector3& field) {
    if (measurements.size() < 3 || !field.allFinite() || field.norm() == 0.0) {
      return std::nullopt;
    }
    std::vector<Acts::Vector3> positions;
    positions.reserve(measurements.size());
    for (const auto& measurement : measurements) {
      if (measurement.surface == nullptr || !measurement.local.allFinite()) {
        return std::nullopt;
      }
      const auto position =
          measurement.surface->localToGlobal(gctx, measurement.local, Acts::Vector3::UnitZ());
      if (!position.allFinite()) {
        return std::nullopt;
      }
      positions.push_back(position);
    }
    const auto free = Acts::estimateTrackParamsFromSpacePoints(positions, field);
    if (!free.ok() || !free->allFinite()) {
      return std::nullopt;
    }
    const Acts::Vector3 direction = free->segment<3>(Acts::eFreeDir0);
    const auto local              = measurements.front().surface->globalToLocal(
        gctx, free->segment<3>(Acts::eFreePos0), direction);
    if (!local.ok()) {
      return std::nullopt;
    }
    Parameters result;
    result << local->x(), local->y(), std::atan2(direction.y(), direction.x()),
        std::atan2(std::hypot(direction.x(), direction.y()), direction.z()),
        (*free)[Acts::eFreeQOverP];
    return result.allFinite() ? std::optional<Parameters>{result} : std::nullopt;
  }

  std::optional<Estimate> estimateWithCovariance(const std::vector<Measurement>& measurements,
                                                 const Acts::GeometryContext& gctx,
                                                 const Acts::Vector3& field) {
    const auto parameters = estimateParameters(measurements, gctx, field);
    if (!parameters) {
      return std::nullopt;
    }
    Estimate result{.parameters = *parameters, .covariance = Covariance::Zero()};
    for (std::size_t hit = 0; hit < measurements.size(); ++hit) {
      const auto& covariance = measurements[hit].covariance;
      if (!covariance.allFinite() || !covariance.isApprox(covariance.transpose()) ||
          Eigen::LLT<Acts::SquareMatrix2>(covariance).info() != Eigen::Success) {
        return std::nullopt;
      }
      Eigen::Matrix<double, 5, 2> jacobian;
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        auto plus         = measurements;
        auto minus        = measurements;
        const double step = std::max(1.0e-4 * Acts::UnitConstants::mm,
                                     1.0e-2 * std::sqrt(covariance(coordinate, coordinate)));
        plus[hit].local[coordinate] += step;
        minus[hit].local[coordinate] -= step;
        const auto upper = estimateParameters(plus, gctx, field);
        const auto lower = estimateParameters(minus, gctx, field);
        if (!upper || !lower) {
          return std::nullopt;
        }
        Parameters difference    = *upper - *lower;
        difference(2)            = std::remainder(difference(2), 2.0 * std::acos(-1.0));
        jacobian.col(coordinate) = difference / (2.0 * step);
      }
      result.covariance += jacobian * covariance * jacobian.transpose();
    }
    result.covariance = (0.5 * (result.covariance + result.covariance.transpose())).eval();
    return result.covariance.allFinite() ? std::optional<Estimate>{result} : std::nullopt;
  }
} // namespace b0acts

void B0TrackerActsSeeding::init() {
  if (!(m_cfg.covarianceInflation > 0.0) || !std::isfinite(m_cfg.covarianceInflation) ||
      !(m_cfg.timeVariance > 0.0) || !std::isfinite(m_cfg.timeVariance) ||
      !(m_cfg.minAbsField > 0.0) || !std::isfinite(m_cfg.minAbsField) ||
      !(m_cfg.scatteringScale >= 0.0) || !std::isfinite(m_cfg.scatteringScale) ||
      !(m_cfg.qOverPRelativeUncertainty >= 0.0) ||
      !std::isfinite(m_cfg.qOverPRelativeUncertainty) || !(m_cfg.minCurvatureSignificance >= 0.0) ||
      !std::isfinite(m_cfg.minCurvatureSignificance)) {
    throw std::invalid_argument("B0 ACTS seeding requires finite, valid covariance/field settings");
  }
  m_geometry = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (!m_geometry) {
    throw std::runtime_error("B0 ACTS seeding requires ACTS geometry");
  }
  // A geometry without B0 has no candidates, so the angle is not needed there.
  try {
    m_crossingAngle = m_geometry->dd4hepDetector()->constant<double>("CrossingAngle") / dd4hep::rad;
  } catch (const std::exception&) {
    m_crossingAngle = 0.0;
  }
}

void B0TrackerActsSeeding::process(const Input& input, const Output& output) const {
  const auto [candidates, measurements] = input;
  auto [seeds, parameters]              = output;
  using HitKey                          = std::pair<unsigned int, int>;
  const auto key                        = [](const auto& hit) -> HitKey {
    const auto id = hit.getObjectID();
    return {id.collectionID, id.index};
  };
  std::map<HitKey, std::vector<std::size_t>> byHit;
  for (std::size_t index = 0; index < measurements->size(); ++index) {
    for (const auto hit : (*measurements)[index].getHits()) {
      byHit[key(hit)].push_back(index);
    }
  }
  const auto& gctx         = m_geometry->getActsGeometryContext();
  const auto fieldProvider = m_geometry->getFieldProvider();
  auto fieldCache          = fieldProvider->makeCache(m_geometry->getActsMagneticFieldContext());
  const double ca          = std::cos(m_crossingAngle);
  const double sa          = std::sin(m_crossingAngle);
  std::set<std::vector<HitKey>> usedCandidates;
  for (const auto candidate : *candidates) {
    if (seeds->size() >= m_cfg.maxSeeds) {
      break;
    }
    std::vector<edm4eic::TrackerHit> hits(candidate.getHits().begin(), candidate.getHits().end());
    std::vector<HitKey> identity;
    for (const auto hit : hits) {
      identity.push_back(key(hit));
    }
    std::ranges::sort(identity);
    // Opposite-charge input hypotheses contain the same hits. Re-estimate
    // once, then determine charge ambiguity from the multipoint covariance.
    if (!usedCandidates.insert(identity).second || identity.size() < 3 ||
        std::ranges::adjacent_find(identity) != identity.end()) {
      continue;
    }
    std::ranges::stable_sort(hits, [&](const auto& left, const auto& right) {
      const auto a = left.getPosition();
      const auto b = right.getPosition();
      return a.x * sa + a.z * ca < b.x * sa + b.z * ca;
    });
    std::vector<b0acts::Measurement> points;
    std::vector<b0stub::Point3> ion;
    std::set<std::uint64_t> surfaces;
    for (const auto hit : hits) {
      const auto found = byHit.find(key(hit));
      if (found == byHit.end() || found->second.size() != 1) {
        break;
      }
      const auto measurement = (*measurements)[found->second.front()];
      const auto gid         = Acts::GeometryIdentifier{measurement.getSurface()};
      const auto* surface    = m_geometry->trackingGeometry()->findSurface(gid);
      if (surface == nullptr || gid.sensitive() == 0 || !surfaces.insert(gid.value()).second) {
        break;
      }
      b0acts::Measurement point;
      point.surface = surface;
      point.local << measurement.getLoc().a, measurement.getLoc().b;
      const auto covariance = measurement.getCovariance();
      point.covariance << covariance.xx, covariance.xy, covariance.xy, covariance.yy;
      points.push_back(point);
      const auto global = surface->localToGlobal(gctx, point.local, Acts::Vector3::UnitZ());
      ion.push_back(
          {global.x() * ca - global.z() * sa, global.y(), global.x() * sa + global.z() * ca});
    }
    if (points.size() != hits.size()) {
      debug("Skipping B0 candidate with missing, ambiguous or repeated sensor measurements");
      continue;
    }
    const auto path = b0stub::fitStub(ion);
    if (!path.valid) {
      continue;
    }
    const double z = 0.5 * (ion.front().z + ion.back().z);
    const Acts::Vector3 midpoint(path.x(z) * ca + z * sa, path.y(z), -path.x(z) * sa + z * ca);
    const auto field = fieldProvider->getField(midpoint, fieldCache);
    if (!field.ok() || !field->allFinite() ||
        field->norm() < m_cfg.minAbsField * Acts::UnitConstants::T) {
      debug("Skipping B0 candidate with an invalid or insufficient field");
      continue;
    }
    const auto estimate = b0acts::estimateWithCovariance(points, gctx, *field);
    if (!estimate || (*estimate).parameters(4) == 0.0) {
      debug("Skipping B0 candidate whose multipoint estimate/covariance failed");
      continue;
    }
    const auto& state   = estimate->parameters;
    const double qOverP = state(4);
    const auto charges =
        b0stub::chargeHypotheses(qOverP, estimate->covariance(4, 4), m_cfg.minCurvatureSignificance,
                                 qOverP < 0.0 ? -1 : 1, m_cfg.testBothCharges, m_cfg.charge);
    for (const int charge : charges) {
      if (seeds->size() >= m_cfg.maxSeeds) {
        break;
      }
      const double signedQOverP = charge * std::abs(qOverP);
      auto covariance           = estimate->covariance;
      if (signedQOverP * qOverP < 0.0) {
        covariance.row(4) *= -1.0;
        covariance.col(4) *= -1.0;
      }
      const auto model = b0stub::scatteringCovarianceAdditions(
          signedQOverP, state(3), m_cfg.scatteringScale, m_cfg.qOverPRelativeUncertainty);
      for (int i = 0; i < 3; ++i) {
        covariance(i + 2, i + 2) += model.at(i);
      }
      covariance *= m_cfg.covarianceInflation;
      if (Eigen::LLT<b0acts::Covariance>(covariance).info() != Eigen::Success) {
        debug("Skipping B0 seed with a non-positive covariance");
        continue;
      }
      auto parameter = parameters->create();
      parameter.setType(-1);
      parameter.setSurface(points.front().surface->geometryId().value());
      parameter.setLoc({static_cast<float>(state(0)), static_cast<float>(state(1))});
      parameter.setPhi(static_cast<float>(state(2)));
      parameter.setTheta(static_cast<float>(state(3)));
      parameter.setQOverP(static_cast<float>(signedQOverP));
      const auto first =
          points.front().surface->localToGlobal(gctx, points.front().local, Acts::Vector3::UnitZ());
      // Approximate relativistic time of flight; B0 timing is not fitted here.
      parameter.setTime(candidate.getParams().getTime() + first.norm() / Acts::UnitConstants::ns);
      edm4eic::Cov6f packed;
      for (int row = 0; row < 5; ++row) {
        for (int col = row; col < 5; ++col) {
          packed(row, col) = static_cast<float>(covariance(row, col));
        }
      }
      packed(5, 5) = m_cfg.timeVariance;
      parameter.setCovariance(packed);
      auto seed = seeds->create();
      seed.setParams(parameter);
      seed.setQuality(candidate.getQuality());
      // This metadata remains the eventual perigee target; parameter.surface
      // identifies the actual sensor-bound seed reference.
      seed.setPerigee(candidate.getPerigee());
      for (const auto hit : hits) {
        seed.addToHits(hit);
      }
    }
  }
}
} // namespace eicrecon
