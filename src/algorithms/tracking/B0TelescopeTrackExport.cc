// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include "B0TelescopeTrackExport.h"
#include "B0MeasurementTruth.h"
#include "algorithms/interfaces/ActsSvc.h"
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <edm4hep/MCParticle.h>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/ProxyAccessor.hpp>
#include <Acts/Utilities/UnitVectors.hpp>
#include <ActsExamples/EventData/IndexSourceLink.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <edm4eic/Cov6f.h>
#include <cmath>
#include <map>
#include <algorithm>
#include <array>
#include <vector>
#include <stdexcept>
#include <utility>
namespace eicrecon {
void B0TelescopeTrackExport::init() {
  m_provider = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (!m_provider) throw std::runtime_error("B0 export requires ACTS geometry context");
}
void B0TelescopeTrackExport::process(const Input& input, const Output& output) const {
  const auto [measurements, seeds, states, inputTracks, rawAssociations] = input;
  auto [trajectories, parameters, tracks, links, associations] = output;
  ActsExamples::ConstTrackContainer container(std::make_shared<Acts::ConstVectorTrackContainer>(*inputTracks),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(*states));
  const auto& context = m_provider->getActsGeometryContext();
  Acts::ConstProxyAccessor<unsigned int> seedColumn("seed");
  using Id = b0::TruthId;
  b0::RawTruthVotes rawTruth;
  std::map<Id, edm4hep::MCParticle> particles;
  for (const auto& association : *rawAssociations) {
    const double weight = association.getWeight();
    if (!std::isfinite(weight) || weight < 0) throw std::runtime_error("Invalid B0 raw association");
    const auto particle = association.getSimHit().getParticle();
    if (!particle.isAvailable() || weight == 0) continue;
    const auto id = b0::truthId(particle);
    particles[id] = particle;
    rawTruth[b0::truthId(association.getRawHit())][id] += weight;
  }
  const std::array<double, 6> units{Acts::UnitConstants::mm, Acts::UnitConstants::mm, 1, 1,
                                    1 / Acts::UnitConstants::GeV, Acts::UnitConstants::ns};
  for (const auto& track : container) {
    const auto summary = Acts::MultiTrajectoryHelpers::trajectoryState(container.trackStateContainer(), track.tipIndex());
    const auto& surface = track.referenceSurface();
    if (surface.type() != Acts::Surface::Plane) throw std::runtime_error("B0 local export expects a plane");
    auto trajectory = trajectories->create();
    trajectory.setNMeasurements(summary.nMeasurements); trajectory.setNStates(summary.nStates);
    trajectory.setNOutliers(summary.nOutliers); trajectory.setNHoles(summary.nHoles); trajectory.setNSharedHits(summary.nSharedHits);
    const auto seedIndex = seedColumn(track);
    if (seedIndex >= seeds->size()) throw std::runtime_error("B0 output lost its original seed identity");
    trajectory.setSeed((*seeds)[seedIndex]);
    for (double chi2 : summary.measurementChi2) trajectory.addToMeasurementChi2(chi2);
    for (double chi2 : summary.outlierChi2) trajectory.addToOutlierChi2(chi2);
    const auto& p = track.parameters();
    auto pars = parameters->create(); pars.setType(0); pars.setSurface(surface.geometryId().value());
    pars.setLoc({static_cast<float>(p(0) / units[0]), static_cast<float>(p(1) / units[1])});
    pars.setPhi(p(2)); pars.setTheta(p(3)); pars.setQOverP(p(4) / units[4]); pars.setTime(p(5) / units[5]);
    pars.setPdg(track.particleHypothesis().absolutePdg());
    Eigen::Matrix<double, 6, 6> boundCov;
    edm4eic::Cov6f packed;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
      boundCov(i, j) = track.covariance()(i, j) / units[i] / units[j];
      packed(i, j) = boundCov(i, j);
    }
    pars.setCovariance(packed); trajectory.addToTrackParameters(pars);
    auto result = tracks->create(); result.setTrajectory(trajectory);
    const Acts::Vector3 direction = Acts::makeDirectionFromPhiTheta(p(2), p(3));
    const Acts::Vector3 position = surface.localToGlobal(context, p.head<2>(), direction) / units[0];
    result.setPosition({static_cast<float>(position.x()), static_cast<float>(position.y()), static_cast<float>(position.z())});
    const double q = pars.getQOverP();
    if (!std::isfinite(q) || q == 0.0) {
      result.setType(-1);
      result.setMomentum({0, 0, 0});
      warning("B0 output has unresolved exact-zero/non-finite q/p; use bound parameters, not Cartesian momentum");
    } else {
      result.setType(0);
      const double momentum = std::abs(1 / q);
      const Acts::Vector3 vector = momentum * direction;
      result.setMomentum({static_cast<float>(vector.x()), static_cast<float>(vector.y()), static_cast<float>(vector.z())});
      result.setCharge(std::copysign(1.0, q));
      Eigen::Matrix<double, 6, 6> jacobian = Eigen::Matrix<double, 6, 6>::Zero();
      jacobian.block<3, 2>(0, 0) = surface.transform(context).linear().leftCols<2>();
      jacobian.block<3, 1>(3, 2) = momentum * Acts::Vector3(-std::sin(p(3)) * std::sin(p(2)), std::sin(p(3)) * std::cos(p(2)), 0);
      jacobian.block<3, 1>(3, 3) = momentum * Acts::Vector3(std::cos(p(3)) * std::cos(p(2)), std::cos(p(3)) * std::sin(p(2)), -std::sin(p(3)));
      jacobian.block<3, 1>(3, 4) = -momentum / q * direction;
      const Eigen::Matrix<double, 6, 6> cartesian = jacobian * boundCov * jacobian.transpose();
      edm4eic::Cov6f cartesianPacked;
      for (int i = 0; i < 6; ++i) for (int j = i; j < 6; ++j) cartesianPacked(i, j) = cartesian(i, j);
      result.setPositionMomentumCovariance(cartesianPacked);
    }
    result.setTime(pars.getTime()); result.setTimeError(std::sqrt(std::max(0.0, boundCov(5, 5))));
    result.setChi2(summary.chi2Sum); result.setNdf(summary.NDF); result.setPdg(pars.getPdg());
    b0::TruthVotes weights;
    for (const auto& state : track.trackStatesReversed()) {
      const auto flags = state.typeFlags();
#if Acts_VERSION_MAJOR >= 45
      if (!flags.isMeasurement() || flags.isHole() || flags.isOutlier()) continue;
#else
      if (!flags.test(Acts::TrackStateFlag::MeasurementFlag) || flags.test(Acts::TrackStateFlag::HoleFlag) || flags.test(Acts::TrackStateFlag::OutlierFlag)) continue;
#endif
      const auto index = state.getUncalibratedSourceLink().template get<ActsExamples::IndexSourceLink>().index();
      const auto measurement = (*measurements)[index];
      result.addToMeasurements(measurement);
      for (const auto& [id, weight] : b0::measurementTruth(measurement, rawTruth)) weights[id] += weight;
    }
    // Equal measurement weight, not equal raw-pad weight. Unassociated signal
    // remains in the denominator even when a cluster has a small true tail.
    if (summary.nMeasurements > 0) for (const auto& [id, sum] : weights) {
      const double weight = sum / summary.nMeasurements;
      auto link = links->create(); link.setFrom(result); link.setTo(particles.at(id)); link.setWeight(weight);
      auto association = associations->create(); association.setRec(result); association.setSim(particles.at(id)); association.setWeight(weight);
    }
  }
}
} // namespace eicrecon
