// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#include "B0TelescopeSeeding.h"

#if EICRECON_HAS_B0_TELESCOPE
#include "B0TelescopeMath.h"
#include "B0TelescopeSeedFinder.h"
#include "B0TrackParameterEstimator.h"
#include "algorithms/interfaces/ActsSvc.h"
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/SpacePointColumns.hpp>
#include <DD4hep/IDDescriptor.h>
#include <DD4hep/Readout.h>
#include <DDSegmentation/BitFieldCoder.h>
#include <Evaluator/DD4hepUnits.h>
#include <edm4eic/Cov6f.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace eicrecon {

void B0TelescopeSeeding::init() {
  if (m_cfg.layersPerStation != 1 && m_cfg.layersPerStation != 2) {
    throw std::invalid_argument("B0 layersPerStation must be 1 (ideal) or 2 (front/back)");
  }
  for (const double value : {m_cfg.minDeltaZ, m_cfg.maxSlope, m_cfg.maxResidualX,
                             m_cfg.maxResidualY, m_cfg.positionSigma, m_cfg.angleSigma,
                             m_cfg.relativeMomentumSigma, m_cfg.qOverPSigmaFloor,
                             m_cfg.timeSigma, m_cfg.propagationPathLimit}) {
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument("B0 telescope cuts and seed uncertainties must be finite and positive");
    }
  }
  if (m_cfg.maxCombinations == 0 || m_cfg.maxSeeds == 0 || m_cfg.maxSeedsPerMiddle == 0) {
    throw std::invalid_argument("B0 telescope candidate/seed limits must be positive");
  }
  m_geo = algorithms::ActsSvc::instance().acts_geometry_provider();
  const auto* detector = m_geo->dd4hepDetector();
  m_decoder = detector->readout("B0TrackerHits").idSpec().decoder();
  (void)m_decoder->index("layer");
  (void)m_decoder->index("system");
  m_system = detector->constant<long>("B0Tracker_Station_1_ID");
  m_rotation = detector->constant<double>("B0Tracker_rotation") / dd4hep::rad;
  if (!std::isfinite(m_rotation)) {
    throw std::invalid_argument("B0 detector rotation must be finite");
  }
}

void B0TelescopeSeeding::process(const Input& input, const Output& output) const {
  const auto [hits] = input;
  auto [outSeeds, outParameters] = output;
  using namespace b0telescope;
  if (hits->empty()) {
    return;
  }
  if (hits->size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::length_error("B0 hit collection exceeds ACTS index range");
  }
  std::array<std::vector<std::uint32_t>, 4> stations;
  std::size_t invalidHits = 0;
  for (std::uint32_t i = 0; i < hits->size(); ++i) {
    const auto hit = (*hits)[i];
    const auto p = hit.getPosition();
    const auto station = stationFromLayer(
        static_cast<int>(m_decoder->get(hit.getCellID(), "layer")), m_cfg.layersPerStation);
    if (!station || m_decoder->get(hit.getCellID(), "system") != m_system ||
        !finite({p.x, p.y, p.z}) || !std::isfinite(hit.getTime()) ||
        std::hypot(p.x, p.y) <= std::numeric_limits<float>::epsilon() * edm4eic::unit::mm) {
      ++invalidHits;
      continue;
    }
    stations[*station].push_back(i);
  }
  if (invalidHits != 0) {
    warning("B0 telescope rejected {} hits with invalid identity/position/time", invalidHits);
  }
  const std::array<std::size_t, 4> counts{
      stations[0].size(), stations[1].size(), stations[2].size(), stations[3].size()};
  if (!withinCombinationBudget(counts, m_cfg.maxCombinations)) {
    warning("B0 telescope combination budget exceeded: station hits [{}, {}, {}, {}]; "
            "event produces no seeds (no partial prefix kept)", counts[0], counts[1], counts[2], counts[3]);
    return;
  }
  // Sort by stable data, never addresses. Keep indices, not pointers to PODIO proxies.
  const auto key = [&](std::uint32_t i) {
    const auto hit = (*hits)[i];
    const auto p = hit.getPosition();
    const auto id = hit.getObjectID();
    return std::tuple{hit.getCellID(), p.x, p.y, p.z, hit.getTime(), id.collectionID, id.index};
  };
  SpacePoints points(Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR);
  points.reserve(hits->size());
  std::vector<std::uint32_t> hitIndices;
  hitIndices.reserve(hits->size());
  std::array<std::uint32_t, 5> offsets{};
  constexpr double mm = Acts::UnitConstants::mm / edm4eic::unit::mm;
  for (unsigned station = 0; station < 4; ++station) {
    std::sort(stations[station].begin(), stations[station].end(),
               [&](auto a, auto b) { return key(a) < key(b); });
    for (const auto index : stations[station]) {
      const auto p = (*hits)[index].getPosition();
      auto point = points.createSpacePoint();
      point.xy() = {static_cast<float>(p.x * mm), static_cast<float>(p.y * mm)};
      point.zr() = {static_cast<float>(p.z * mm), static_cast<float>(std::hypot(p.x, p.y) * mm)};
      hitIndices.push_back(index);
    }
    offsets[station + 1] = static_cast<std::uint32_t>(points.size());
  }
  Seeds seeds;
  seeds.assignSpacePointContainer(points);
  findTelescopeSeeds(points, offsets, m_rotation, m_cfg, seeds);

  // Global deterministic ranking; apply an event cap only after considering all stations.
  using Candidate = std::pair<float, std::array<std::uint32_t, 3>>;
  std::vector<Candidate> candidates;
  for (const auto seed : seeds) {
    const auto indices = seed.spacePointIndices();
    candidates.emplace_back(-seed.quality(), std::array<std::uint32_t, 3>{indices[0], indices[1], indices[2]});
  }
  std::sort(candidates.begin(), candidates.end());
  candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
  const auto candidateCount = candidates.size();
  if (candidateCount > m_cfg.maxSeeds) {
    warning("B0 telescope seed cap keeps {} of {} ranked triplets", m_cfg.maxSeeds, candidateCount);
    candidates.resize(m_cfg.maxSeeds);
  }
  std::size_t failedEstimates = 0;
  for (const auto& [score, indices] : candidates) {
    std::array<Acts::Vector3, 3> positions;
    for (std::size_t i = 0; i < 3; ++i) {
      const auto p = (*hits)[hitIndices[indices[i]]].getPosition();
      positions[i] = Acts::Vector3{p.x * mm, p.y * mm, p.z * mm};
    }
    const auto estimate = estimateTrackParameters(
        positions, (*hits)[hitIndices[indices[0]]].getTime() * Acts::UnitConstants::ns / edm4eic::unit::ns,
        m_geo->getActsGeometryContext(), m_geo->getActsMagneticFieldContext(),
        m_geo->getFieldProvider(), m_cfg);
    if (!estimate) {
      ++failedEstimates;
      continue;
    }
    const std::array<double, 6> toEdm{
        edm4eic::unit::mm / Acts::UnitConstants::mm, edm4eic::unit::mm / Acts::UnitConstants::mm,
        1.0 / Acts::UnitConstants::rad, 1.0 / Acts::UnitConstants::rad,
        Acts::UnitConstants::GeV / edm4eic::unit::GeV, edm4eic::unit::ns / Acts::UnitConstants::ns};
    const auto& p = estimate->parameters;
    auto parameters = outParameters->create();
    parameters.setType(-1);
    parameters.setLoc({static_cast<float>(p[0] * toEdm[0]), static_cast<float>(p[1] * toEdm[1])});
    parameters.setPhi(static_cast<float>(p[2] * toEdm[2]));
    parameters.setTheta(static_cast<float>(p[3] * toEdm[3]));
    parameters.setQOverP(static_cast<float>(p[4] * toEdm[4]));
    parameters.setTime(static_cast<float>(p[5] * toEdm[5]));
    edm4eic::Cov6f covariance{};
    for (std::size_t i = 0; i < 6; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        covariance(i, j) = static_cast<float>(estimate->covariance(i, j) * toEdm[i] * toEdm[j]);
      }
    }
    parameters.setCovariance(covariance);
    auto seed = outSeeds->create();
    seed.setPerigee({0.F, 0.F, 0.F});
    seed.setParams(parameters);
    seed.setQuality(-score);
    for (const auto index : indices) {
      seed.addToHits((*hits)[hitIndices[index]]);
    }
  }
  debug("B0 telescope: stations [{}, {}, {}, {}], {} ranked triplets, {} estimate/transport "
        "failures, {} output seeds", counts[0], counts[1], counts[2], counts[3], candidateCount,
        failedEstimates, outSeeds->size());
}

} // namespace eicrecon
#endif
