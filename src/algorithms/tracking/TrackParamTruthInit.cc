// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten, Dmitry Romanov

#include "TrackParamTruthInit.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Evaluator/DD4hepUnits.h>
#include <algorithms/logger.h>
#include <edm4eic/Cov6f.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3d.h>
#include <edm4hep/Vector3f.h>
#include <podio/ObjectID.h>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>

#include "extensions/spdlog/SpdlogFormatters.h" // IWYU pragma: keep

namespace eicrecon {

void TrackParamTruthInit::process(const Input& input, const Output& output) const {
  // MCParticles uses numerical values in its specified units,
  // while m_cfg is in the DD4hep unit system

  const auto [headers, mcparticles, simhit_collections] = input;
  auto [track_seeds, track_parameters]                  = output;

  // local random generator
  auto seed = m_uid.getUniqueID(*headers, "TrackParamTruthInit");
  std::default_random_engine generator(seed);
  std::normal_distribution<double> gaussian;

  // Earliest hit of each particle (by time), keyed by the particle's ObjectID.
  // Hits flagged as produced by a secondary carry the secondary's momentum and are skipped.
  std::map<podio::ObjectID, edm4hep::SimTrackerHit> first_hits;
  if (m_cfg.anchorAtFirstHit) {
    if (simhit_collections.empty()) {
      throw std::runtime_error(
          "TrackParamTruthInit: anchorAtFirstHit requires at least one sim hit collection");
    }
    for (const auto& simhits : simhit_collections) {
      for (const auto& hit : *simhits) {
        if (hit.isProducedBySecondary()) {
          continue;
        }
        const auto id = hit.getParticle().getObjectID();
        auto it       = first_hits.find(id);
        if (it == first_hits.end() || hit.getTime() < it->second.getTime()) {
          first_hits.insert_or_assign(id, hit);
        }
      }
    }
  }

  // Loop over input particles
  for (const auto& mcparticle : *mcparticles) {

    // Accept all generated particles that been transported in simulation.
    // Both primary / decayed signal or beam background particles can produce truth seeds.
    // SimulatorStatus will stay 0 for intermediate partons, resonances, or any particles
    // that are not transported in simulation hence have no chance to leave hits.
    if (mcparticle.getSimulatorStatus() == 0) {
      trace("ignoring particle with generatorStatus = {}, simulatorStatus = {}",
            mcparticle.getGeneratorStatus(), mcparticle.getSimulatorStatus());
      continue;
    }

    // require close to interaction vertex
    auto v = mcparticle.getVertex();
    if (std::abs(v.x) * dd4hep::mm > m_cfg.maxVertexX ||
        std::abs(v.y) * dd4hep::mm > m_cfg.maxVertexY ||
        std::abs(v.z) * dd4hep::mm > m_cfg.maxVertexZ) {
      trace("ignoring particle with vs = {} [mm]", v);
      continue;
    }

    // require minimum momentum
    const auto& p   = mcparticle.getMomentum();
    const auto pmag = std::hypot(p.x, p.y, p.z);
    if (pmag * dd4hep::GeV < m_cfg.minMomentum) {
      trace("ignoring particle with p = {} GeV ", pmag);
      continue;
    }

    // require minimum pseudorapidity
    const auto phi   = std::atan2(p.y, p.x);
    const auto theta = std::atan2(std::hypot(p.x, p.y), p.z);
    const auto eta   = -std::log(std::tan(theta / 2));
    if (eta > m_cfg.maxEtaForward || eta < -std::abs(m_cfg.maxEtaBackward)) {
      trace("ignoring particle with Eta = {}", eta);
      continue;
    }

    // get the particle charge
    // note that we cannot trust the mcparticles charge, as DD4hep
    // sets this value to zero! let's lookup by PDGID instead
    const auto pdg      = mcparticle.getPDG();
    const auto particle = m_particleSvc.particle(pdg);
    double charge       = std::copysign(1.0, particle.charge);
    if (std::abs(particle.charge) < std::numeric_limits<double>::epsilon()) {
      trace("ignoring neutral particle");
      continue;
    }

    // Kinematics at the seed start: the vertex, or anchorDistance before the particle's
    // earliest hit with the momentum it has at that hit
    double seed_pmag     = pmag;
    double seed_phi      = phi;
    double seed_theta    = theta;
    double seed_time     = mcparticle.getTime();
    Acts::Vector3 anchor = Acts::Vector3::Zero();
    if (m_cfg.anchorAtFirstHit) {
      const auto it = first_hits.find(mcparticle.getObjectID());
      if (it == first_hits.end()) {
        trace("ignoring particle without hits");
        continue;
      }
      const auto& hp = it->second.getMomentum();
      const auto& hx = it->second.getPosition();
      seed_pmag      = std::hypot(hp.x, hp.y, hp.z);
      if (seed_pmag * dd4hep::GeV < m_cfg.minMomentum) {
        trace("ignoring particle with p = {} GeV at its first hit", seed_pmag);
        continue;
      }
      seed_phi   = std::atan2(hp.y, hp.x);
      seed_theta = std::atan2(std::hypot(hp.x, hp.y), hp.z);
      seed_time  = it->second.getTime();
      anchor = Acts::Vector3(hx.x, hx.y, hx.z) -
               (m_cfg.anchorDistance / dd4hep::mm) * Acts::Vector3(hp.x, hp.y, hp.z) / seed_pmag;
    }

    // modify initial momentum to avoid bleeding truth to results when fit fails
    const auto pinit = seed_pmag * (1.0 + m_cfg.momentumSmear * gaussian(generator));

    Acts::Vector3 direction(sin(seed_theta) * cos(seed_phi), sin(seed_theta) * sin(seed_phi),
                            cos(seed_theta));

    // Choose the perigee reference point and compute local (d0, z0) coordinates.
    //
    // Default: back-extrapolate to the transverse PCA w.r.t. the beam axis at z=0.
    //
    // This is correct for IP-originating tracks but produces a reference point
    // that is far away from beam axis for displaced decay products (e.g. from
    // Lambda decay daughters at large z with 25 mrad crossing angle)
    //
    // As this reference point is determined based on a linear extrapolation
    // that is only valid for small extrapolation distances, this approach
    // leads to inconsistencies with the ACTS CKF that must propagate many
    // metres before reaching the detector.
    //
    // useVertexAsPerigee: anchor the perigee at the MCParticle vertex itself so
    // the CKF only needs to propagate forward from the actual production point.
    //
    // anchorAtFirstHit: anchor the perigee just before the particle's first hit, so the CKF
    // starts next to its first measurement (takes precedence over useVertexAsPerigee).
    Acts::Vector3 perigee_center;
    Acts::Vector3 global;

    if (m_cfg.anchorAtFirstHit) {
      perigee_center = anchor;
      global         = anchor; // d0 = z0 = 0
    } else if (m_cfg.useVertexAsPerigee) {
      // Place perigee at the vertex; the particle starts there so d0=z0=0 exactly.
      perigee_center = Acts::Vector3(v.x, v.y, v.z);
      global         = perigee_center; // PCA = vertex itself
    } else {
      // Standard: perigee at beam axis origin, back-extrapolate to transverse PCA.
      perigee_center             = Acts::Vector3(0, 0, 0);
      auto linesurface_parameter = -(v.x * p.x + v.y * p.y) / (p.x * p.x + p.y * p.y);
      global = Acts::Vector3(v.x + linesurface_parameter * p.x, v.y + linesurface_parameter * p.y,
                             v.z + linesurface_parameter * p.z);
    }

    // define line surface for local position values
    auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(perigee_center);

    // convert from global to local coordinates using the defined line surface
    auto local = perigee->globalToLocal(m_geoSvc->getActsGeometryContext(), global, direction);

    if (!local.ok()) {
      error("skipping the track because globaltoLocal function failed");
      continue;
    }

    Acts::Vector2 localpos = local.value();

    // Insert into edm4eic::TrackParameters, which uses numerical values in its specified units
    auto track_parameter = track_parameters->create();
    track_parameter.setType(-1); // type --> seed(-1)
    track_parameter.setLoc({static_cast<float>(localpos(0)),
                            static_cast<float>(localpos(1))}); // 2d location on surface [mm]
    track_parameter.setPhi(seed_phi);                          // phi [rad]
    track_parameter.setTheta(seed_theta);                      // theta [rad]
    track_parameter.setQOverP(charge / (pinit / dd4hep::GeV)); // Q/p [e/GeV]
    track_parameter.setTime(seed_time);                        // time [ns]
    edm4eic::Cov6f cov;
    cov(0, 0) = std::pow(m_cfg.locError / dd4hep::mm, 2);     // loc0 [mm^2]
    cov(1, 1) = std::pow(m_cfg.locError / dd4hep::mm, 2);     // loc1 [mm^2]
    cov(2, 2) = std::pow(m_cfg.phiError / dd4hep::rad, 2);    // phi [rad^2]
    cov(3, 3) = std::pow(m_cfg.thetaError / dd4hep::rad, 2);  // theta [rad^2]
    cov(4, 4) = std::pow(m_cfg.qOverPError * dd4hep::GeV, 2); // qOverP [(e/GeV)^2]
    cov(5, 5) = 10e9;                                         // time
    track_parameter.setCovariance(cov);

    // Insert into edm4eic::TrackSeeds
    auto track_seed = track_seeds->create();
    track_seed.setPerigee({static_cast<float>(perigee_center.x()),
                           static_cast<float>(perigee_center.y()),
                           static_cast<float>(perigee_center.z())});
    track_seed.setParams(track_parameter);
    // There are no hits to store to the seed

    // Debug output
    if (level() <= algorithms::LogLevel::kDebug) {
      debug("Invoke track finding seeded by truth particle with:");
      debug("   p     = {} GeV (smeared to {} GeV)", seed_pmag / dd4hep::GeV, pinit / dd4hep::GeV);
      debug("   q     = {}", charge);
      debug("   q/p   = {} e/GeV (smeared to {} e/GeV)", charge / (seed_pmag / dd4hep::GeV),
            charge / (pinit / dd4hep::GeV));
      debug("   theta = {}", seed_theta);
      debug("   phi   = {}", seed_phi);
      debug("   start = ({}, {}, {}) mm", perigee_center.x(), perigee_center.y(),
            perigee_center.z());
    }
  }
}

} // namespace eicrecon
