// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include "B0TelescopeTruthSeeding.h"
#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/interfaces/UniqueIDGenSvc.h"
#include "services/particle/ParticleSvc.h"
#include <algorithms/geo.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/VolumeManager.h>
#include <edm4hep/MCParticle.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace eicrecon {
namespace {
using Id = std::pair<int, int>;
template <typename T> Id objectId(const T& p) {
  const auto id = p.getObjectID(); return {id.collectionID, id.index};
}
bool supported(int pdg) {
  return pdg == 11 || pdg == -11 || pdg == 13 || pdg == -13 || pdg == 211 || pdg == -211 ||
         pdg == 321 || pdg == -321 || pdg == 2212 || pdg == -2212;
}
}
void B0TelescopeTruthSeeding::init() {
  m_provider = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (!m_provider) throw std::runtime_error("B0 truth reference requires ACTS geometry");
  m_geometry = std::make_unique<b0::TelescopeGeometry>(*m_provider, m_cfg.stationZGap);
  if (m_cfg.minStations < 3 || !(m_cfg.relativeMomentumSmear >= 0) ||
      !std::isfinite(m_cfg.relativeMomentumSmear) || !(m_cfg.measurementPurityMin >= 0.5) ||
      !(m_cfg.measurementPurityMin < 1)) throw std::invalid_argument("Invalid B0 truth configuration");
  for (double value : {m_cfg.minMomentumGeV, m_cfg.positionSigmaMm, m_cfg.slopeSigma,
                      m_cfg.qOverPSigma, m_cfg.timeSigmaNs, m_cfg.maxProjectionMm})
    if (!(value > 0) || !std::isfinite(value)) throw std::invalid_argument("B0 truth scales must be positive/finite");
  if (!m_cfg.referenceFile.empty()) {
    std::ofstream file(m_cfg.referenceFile);
    if (!(file << "{\"type\":\"configuration\",\"schema_version\":1,\"truth_binding\":\"sim_hit_tangent_to_sensor_plane\","
               << "\"min_stations\":" << m_cfg.minStations << ",\"momentum_smear\":" << m_cfg.relativeMomentumSmear
               << ",\"max_projection_mm\":" << m_cfg.maxProjectionMm << "}\n"))
      throw std::runtime_error("Cannot initialize B0 truth reference file");
  }
}
void B0TelescopeTruthSeeding::process(const Input& input, const Output& output) const {
  const auto [headers, simHits, measurements, associations] = input;
  auto [seeds, parameters] = output;
  if (headers->empty()) throw std::runtime_error("Missing B0 truth EventHeader");
  const auto event = (*headers)[0].getEventNumber();
  const auto& gctx = m_provider->getActsGeometryContext();
  const auto& rotation = m_geometry->labToIon();
  const auto* converter = algorithms::GeoSvc::instance().cellIDPositionConverter();
  if (!converter) throw std::runtime_error("Missing B0 cell-ID converter");
  struct Crossing { edm4hep::SimTrackerHit hit; std::uint64_t surface; double z; };
  std::map<Id, std::vector<Crossing>> crossings;
  std::size_t unmapped = 0, invalidMomentum = 0;
  for (const auto& h : *simHits) {
    const auto particle = h.getParticle();
    if (!particle.isAvailable()) continue;
    const Acts::Surface* surface = nullptr;
    try {
      const auto* context = converter->findContext(h.getCellID());
      surface = context ? m_geometry->volumeSurface(context->identifier) : nullptr;
    } catch (const std::exception&) { ++unmapped; continue; }
    if (!surface) { ++unmapped; continue; }
    const auto p = h.getPosition();
    const Eigen::Vector3d ion = rotation * Eigen::Vector3d(p.x, p.y, p.z);
    if (!ion.allFinite()) { ++unmapped; continue; }
    crossings[objectId(particle)].push_back({h, surface->geometryId().value(), ion.z()});
  }
  // Truth association is confined to this explicitly named validation seeder.
  std::map<Id, std::map<Id, double>> rawVotes;
  for (const auto& a : *associations) {
    const auto particle = a.getSimHit().getParticle();
    if (particle.isAvailable() && std::isfinite(a.getWeight()) && a.getWeight() > 0)
      rawVotes[objectId(a.getRawHit())][objectId(particle)] += a.getWeight();
  }
  std::map<Id, std::map<std::uint64_t, std::pair<double, std::size_t>>> matched;
  for (std::size_t i = 0; i < measurements->size(); ++i) {
    const auto m = (*measurements)[i]; std::map<Id, double> votes;
    for (const auto& h : m.getHits()) {
      const auto found = rawVotes.find(objectId(h.getRawHit()));
      if (found == rawVotes.end()) continue;
      double total = 0; for (const auto& [id, w] : found->second) { (void)id; total += w; }
      if (total <= 0) continue;
      for (const auto& [id, w] : found->second) votes[id] += w / total;
    }
    if (m.hits_size() == 0) continue;
    for (const auto& [id, weight] : votes) {
      const double purity = weight / m.hits_size();
      if (purity <= m_cfg.measurementPurityMin) continue;
      auto& bySurface = matched[id]; const auto previous = bySurface.find(m.getSurface());
      if (previous == bySurface.end() || purity > previous->second.first)
        bySurface[m.getSurface()] = {purity, i};
    }
  }
  auto seedNumber = algorithms::UniqueIDGenSvc::instance().getUniqueID(*headers, std::string(name()));
  std::mt19937_64 random(seedNumber); std::normal_distribution<double> normal;
  std::ostringstream records; records << std::setprecision(17);
  for (auto& [id, hits] : crossings) {
    std::sort(hits.begin(), hits.end(), [](const Crossing& a, const Crossing& b) {
      if (a.z != b.z) return a.z < b.z;
      if (a.surface != b.surface) return a.surface < b.surface;
      return a.hit.getTime() < b.hit.getTime();
    });
    const auto particle = hits.front().hit.getParticle(); const int pdg = particle.getPDG();
    std::set<std::size_t> stations;
    for (const auto& h : hits) stations.insert(*m_geometry->stations().stationForSurface(h.surface));
    const auto vertex = particle.getVertex(); const auto originalP = particle.getMomentum();
    const double originalMomentum = std::hypot(originalP.x, originalP.y, originalP.z);
    if (!std::isfinite(originalMomentum) || !std::isfinite(vertex.x) || !std::isfinite(vertex.y) ||
        !std::isfinite(vertex.z)) { ++invalidMomentum; continue; }
    records << "{\"type\":\"particle\",\"event\":" << event << ",\"particle_id\":[" << id.first << ',' << id.second
            << "],\"pdg\":" << pdg << ",\"generator_status\":" << particle.getGeneratorStatus()
            << ",\"n_stations\":" << stations.size() << ",\"supported\":" << (supported(pdg) ? "true" : "false")
            << ",\"vertex_mm\":[" << vertex.x << ',' << vertex.y << ',' << vertex.z
            << "],\"p_production_gev\":" << originalMomentum << "}\n";
    if (!supported(pdg)) continue;
    const double charge = algorithms::ParticleSvc::instance().particle(pdg).charge;
    if (!std::isfinite(charge) || std::abs(std::abs(charge) - 1) > 1e-6) continue;
    bool seeded = false;
    std::set<std::uint64_t> recorded;
    for (const auto& crossing : hits) {
      if (!recorded.insert(crossing.surface).second) continue;
      const auto pos = crossing.hit.getPosition();
      const auto mom = crossing.hit.getMomentum();
      const Eigen::Vector3d point = rotation * Eigen::Vector3d(pos.x, pos.y, pos.z);
      const Eigen::Vector3d momentum = rotation * Eigen::Vector3d(mom.x, mom.y, mom.z);
      const double p = momentum.norm();
      if (!momentum.allFinite() || !(momentum.z() > 0) || !(p > 0)) { ++invalidMomentum; continue; }
      const auto* surface = m_geometry->surface(crossing.surface);
      const Eigen::Vector3d direction = rotation.transpose() * momentum / p;
      const Eigen::Vector3d labPoint(pos.x, pos.y, pos.z);
      const Eigen::Vector3d normalVector = surface->transform(gctx).linear().col(2);
      const double projection = normalVector.dot(surface->center(gctx) / Acts::UnitConstants::mm - labPoint) /
                                 normalVector.dot(direction);
      if (!std::isfinite(projection) || std::abs(projection) > m_cfg.maxProjectionMm) { ++invalidMomentum; continue; }
      b0::TelescopeFit truth;
      truth.valid = true; truth.zReference = point.z();
      truth.state << point.x(), point.y(), momentum.x() / momentum.z(), momentum.y() / momentum.z(), charge / p;
      truth.covariance.diagonal() << std::pow(m_cfg.positionSigmaMm, 2), std::pow(m_cfg.positionSigmaMm, 2),
          std::pow(m_cfg.slopeSigma, 2), std::pow(m_cfg.slopeSigma, 2), std::pow(m_cfg.qOverPSigma, 2);
      edm4eic::TrackParametersCollection truthParameters; auto bound = truthParameters.create();
      const double time = crossing.hit.getTime();
      if (!std::isfinite(time)) { ++invalidMomentum; continue; }
      const double mass = particle.getMass();
      if (!(mass >= 0) || !std::isfinite(mass)) { ++invalidMomentum; continue; }
      const double planeTime = time + projection * std::sqrt(1 + mass * mass / (p * p)) / b0::speedOfLight;
      b0::writeTelescopeSeed(truth, *surface, gctx, rotation, planeTime, std::pow(m_cfg.timeSigmaNs, 2), bound);
      records << "{\"type\":\"truth_state\",\"event\":" << event << ",\"particle_id\":[" << id.first << ',' << id.second
              << "],\"surface\":\"" << crossing.surface << "\",\"projection_mm\":" << projection
              << ",\"state\":[" << bound.getLoc().a << ',' << bound.getLoc().b << ',' << bound.getPhi() << ','
              << bound.getTheta() << ',' << bound.getQOverP() << ',' << bound.getTime() << "]}\n";
      // No IP, eta, generator-status or production-z selection.
      if (seeded || stations.size() < m_cfg.minStations || p < m_cfg.minMomentumGeV) continue;
      const double smear = std::exp(m_cfg.relativeMomentumSmear * normal(random));
      if (!std::isfinite(smear) || smear <= 0) { ++invalidMomentum; continue; }
      truth.state(4) /= smear;
      auto parameter = parameters->create();
      b0::writeTelescopeSeed(truth, *surface, gctx, rotation, planeTime, std::pow(m_cfg.timeSigmaNs, 2), parameter);
      parameter.setPdg(std::abs(pdg));
      auto seed = seeds->create(); seed.setParams(parameter); seed.setPerigee({0, 0, 0});
      const auto accepted = matched.find(id);
      if (accepted != matched.end()) for (const auto& [surfaceId, match] : accepted->second) {
        (void)surfaceId;
        for (const auto& hit : (*measurements)[match.second].getHits()) seed.addToHits(hit);
      }
      records << "{\"type\":\"truth_seed\",\"event\":" << event << ",\"seed\":" << seeds->size() - 1
              << ",\"particle_id\":[" << id.first << ',' << id.second << "]}\n";
      seeded = true;
    }
  }
  records << "{\"type\":\"event\",\"event\":" << event << ",\"seeds\":" << seeds->size()
          << ",\"unmapped_hits\":" << unmapped << ",\"invalid_momentum_or_projection\":" << invalidMomentum << "}\n";
  if (!m_cfg.referenceFile.empty()) {
    std::lock_guard lock(m_fileMutex); std::ofstream out(m_cfg.referenceFile, std::ios::app);
    if (!(out << records.str())) throw std::runtime_error("Cannot append B0 truth reference");
  }
}
} // namespace eicrecon
