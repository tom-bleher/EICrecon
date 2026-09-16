// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include "B0TelescopeSeeding.h"
#include "algorithms/interfaces/ActsSvc.h"
#include <Acts/Definitions/Units.hpp>
#include <edm4eic/Cov6f.h>
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <mutex>
#include <stdexcept>

namespace eicrecon {
namespace {
b0::State boundState(const b0::State& state, double z, const Acts::Surface& surface,
                     const Acts::GeometryContext& context, const Eigen::Matrix3d& rotation) {
  const Acts::Vector3 direction = (rotation.transpose() * Acts::Vector3(state(2), state(3), 1)).normalized();
  Acts::Vector3 point = rotation.transpose() * Acts::Vector3(state(0), state(1), z) * Acts::UnitConstants::mm;
  const Acts::Vector3 normal = surface.transform(context).linear().col(2);
  const double denominator = normal.dot(direction);
  if (std::abs(denominator) < 1e-8) throw std::runtime_error("B0 seed is tangent to its reference plane");
  point += direction * normal.dot(surface.center(context) - point) / denominator;
  const auto local = surface.globalToLocal(context, point, direction);
  if (!local.ok()) throw std::runtime_error("B0 seed cannot be bound to its reference plane");
  b0::State result;
  result << local.value().x() / Acts::UnitConstants::mm, local.value().y() / Acts::UnitConstants::mm,
            std::atan2(direction.y(), direction.x()),
            std::atan2(std::hypot(direction.x(), direction.y()), direction.z()), state(4);
  return result;
}
}
void b0::writeTelescopeSeed(const TelescopeFit& fit, const Acts::Surface& surface,
                            const Acts::GeometryContext& context, const Eigen::Matrix3d& rotation,
                            double time, double timeVariance,
                            edm4eic::MutableTrackParameters parameters) {
  const State nominal = boundState(fit.state, fit.zReference, surface, context, rotation);
  Covariance jacobian;
  for (int i = 0; i < 5; ++i) {
    const double step = std::max(1e-8, std::sqrt(std::max(0.0, fit.covariance(i, i))) * 1e-3);
    State plus = fit.state, minus = fit.state; plus(i) += step; minus(i) -= step;
    State difference = boundState(plus, fit.zReference, surface, context, rotation) -
                       boundState(minus, fit.zReference, surface, context, rotation);
    difference(2) = std::remainder(difference(2), 2 * std::acos(-1));
    jacobian.col(i) = difference / (2 * step);
  }
  Covariance cov = jacobian * fit.covariance * jacobian.transpose();
  cov = (0.5 * (cov + cov.transpose())).eval();
  if (!nominal.allFinite() || !cov.allFinite() || cov.llt().info() != Eigen::Success)
    throw std::runtime_error("Invalid B0 bound seed/covariance");
  parameters.setType(-1);
  parameters.setSurface(surface.geometryId().value());
  parameters.setLoc({static_cast<float>(nominal(0)), static_cast<float>(nominal(1))});
  parameters.setPhi(nominal(2)); parameters.setTheta(nominal(3)); parameters.setQOverP(nominal(4));
  parameters.setTime(std::isfinite(time) ? time : 0.0);
  edm4eic::Cov6f packed;
  for (int i = 0; i < 5; ++i) for (int j = i; j < 5; ++j) packed(i, j) = cov(i, j);
  packed(5, 5) = std::isfinite(timeVariance) && timeVariance > 0 ? timeVariance : 1e6;
  parameters.setCovariance(packed);
}
void B0TelescopeSeeding::init() {
  m_provider = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (!m_provider) throw std::runtime_error("B0: ACTS geometry unavailable");
  m_geometry = std::make_unique<b0::TelescopeGeometry>(*m_provider, m_cfg.stationZGap);
  m_geometry->writeManifest(m_cfg.stationMapFile);
  if (!m_cfg.diagnosticsFile.empty()) {
    std::ofstream file(m_cfg.diagnosticsFile);
    if (!(file << "{\"type\":\"configuration\",\"schema_version\":1,\"max_trials\":"
               << m_cfg.search.maxTrials << ",\"timing\":" << (m_cfg.search.useTiming ? "true" : "false")
               << ",\"field_integral\":" << (m_cfg.search.fit.fieldIntegral ? "true" : "false") << "}\n"))
      throw std::runtime_error("Cannot initialize B0 seed diagnostic file");
  }
  for (double sigma : {m_cfg.seedPositionSigma, m_cfg.seedAngularSigma, m_cfg.seedQOverPSigma})
    if (!(sigma > 0) || !std::isfinite(sigma)) throw std::invalid_argument("B0 seed floors must be finite and positive");
  b0::findTelescopeCandidates({}, [](const Eigen::Vector3d&) { return Eigen::Vector3d::Zero(); }, m_cfg.search);
  b0::chargeProposals({}, m_cfg.curvatureSignificance, m_cfg.proposalMaxMomentum, m_cfg.forceCharge, m_cfg.bothCharges);
}
void B0TelescopeSeeding::process(const Input& input, const Output& output) const {
  const auto [headers, measurements] = input;
  if (headers->empty()) throw std::runtime_error("Missing B0 EventHeader");
  auto [seeds, parameters] = output;
  const auto points = m_geometry->measurements(*measurements);
  const auto fieldProvider = m_provider->getFieldProvider();
  auto cache = fieldProvider->makeCache(m_provider->getActsMagneticFieldContext());
  const auto& rotation = m_geometry->labToIon();
  b0::Field field = [&](const Eigen::Vector3d& p) -> Eigen::Vector3d {
    const auto b = fieldProvider->getField(rotation.transpose() * p * Acts::UnitConstants::mm, cache);
    if (!b.ok()) return Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
    return rotation * b.value() / Acts::UnitConstants::T;
  };
  const auto search = b0::findTelescopeCandidates(points, field, m_cfg.search);
  for (const auto& candidate : search.candidates) {
    const auto first = candidate.hits.front();
    const auto measurement = (*measurements)[first];
    const auto* surface = m_geometry->surface(measurement.getSurface());
    for (auto proposal : b0::chargeProposals(candidate.fit, m_cfg.curvatureSignificance,
                                             m_cfg.proposalMaxMomentum, m_cfg.forceCharge, m_cfg.bothCharges)) {
      proposal.covariance(0, 0) += std::pow(m_cfg.seedPositionSigma, 2);
      proposal.covariance(1, 1) += std::pow(m_cfg.seedPositionSigma, 2);
      proposal.covariance(2, 2) += std::pow(m_cfg.seedAngularSigma, 2);
      proposal.covariance(3, 3) += std::pow(m_cfg.seedAngularSigma, 2);
      proposal.covariance(4, 4) += std::pow(m_cfg.seedQOverPSigma, 2);
      auto parameter = parameters->create();
      b0::writeTelescopeSeed(proposal, *surface, m_provider->getActsGeometryContext(), rotation,
                             measurement.getTime(), measurement.getCovariance().zz, parameter);
      auto seed = seeds->create(); seed.setParams(parameter);
      seed.setQuality(-candidate.fit.chi2 / std::max(1, candidate.fit.ndf));
      seed.setPerigee({0, 0, 0}); // unused when params.surface is nonzero
      for (auto index : candidate.hits) for (const auto& hit : (*measurements)[index].getHits()) seed.addToHits(hit);
    }
  }
  if (!m_cfg.diagnosticsFile.empty()) {
    std::lock_guard lock(m_fileMutex);
    std::ofstream file(m_cfg.diagnosticsFile, std::ios::app);
    file << "{\"type\":\"seed_event\",\"event\":" << (*headers)[0].getEventNumber()
         << ",\"hits\":" << points.size() << ",\"seeds\":" << seeds->size()
         << ",\"trials\":" << search.trials << ",\"fits\":" << search.fits
         << ",\"truncated\":" << (search.truncated ? "true" : "false")
         << ",\"timing_rejected\":" << search.timingRejected
         << ",\"candidate_cap_drops\":" << search.candidatesDropped << "}\n";
    if (!file) throw std::runtime_error("Cannot append B0 seed diagnostics");
  }
  debug("B0 local seeds: hits={} seeds={} trials={} fits={} truncated={} timingRejected={} candidateCapDrops={}",
        points.size(), seeds->size(), search.trials, search.fits, search.truncated,
        search.timingRejected, search.candidatesDropped);
}
} // namespace eicrecon
