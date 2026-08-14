// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
//
// Stub seeder for the B0 tracker.
//
// Why this exists: the B0 sits inside the B0pf dipole (B along y). The
// orthogonal seeder's solenoid helix model returns invalid parameters there
// (charge sign is a coin flip, |p| spans 30 MeV - 77 GeV on 41 GeV protons),
// so the CKF can never attach hits. Here we instead exploit the spectrometer
// layout directly:
//
//   1. rotate hits into the ion frame (crossing angle about y), where the
//      dipole bend lies purely in the x-z plane;
//   2. per candidate (one hit per station): fit y(z) with a line and x(z)
//      with a parabola; the curvature gives |p| (sagitta), the tangent the
//      direction;
//   3. back-extrapolate analytically: parabola to the field entrance, then a
//      straight line to the origin; express the result on the origin perigee
//      surface, which is what CKFTracking hard-codes as seed reference.
//
// The truth-seeded chain reconstructs these tracks at 84.6% band efficiency
// from exactly such origin-perigee parameters, so approximate stub seeds are
// sufficient; the CKF does the rest.

#include "B0TrackerStubSeeder.h"

#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/tracking/ActsGeometryProvider.h"

#include <algorithms/geo.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <DD4hep/Detector.h>
#include <DD4hep/Segmentations.h>
#include <DDRec/CellIDPositionConverter.h>
#include <Eigen/Dense>
#include <edm4eic/Cov6f.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace eicrecon {

namespace {

  struct StubCandidate {
    std::vector<std::size_t> hitIndices;
    double chi2{0.0};
    // ion-frame fit results
    double c0{0}, c1{0}, c2{0}; // x'(z') = c0 + c1 z' + c2 z'^2
    double b0{0}, b1{0};        // y'(z') = b0 + b1 z'
  };

} // namespace

void B0TrackerStubSeeder::init() {
  const auto& geo = algorithms::GeoSvc::instance();
  m_converter     = geo.cellIDPositionConverter();
  // Fail early if this compact has no B0 tracker readout.
  (void)geo.detector()->readout(m_cfg.readout);
  m_acts_context = algorithms::ActsSvc::instance().acts_geometry_provider();

  if (m_converter == nullptr || m_acts_context == nullptr) {
    throw std::runtime_error("B0TrackerStubSeeder: required geometry service is unavailable");
  }
}

void B0TrackerStubSeeder::process(const Input& input, const Output& output) const {
  const auto [measurements]         = input;
  auto [seeds, track_params_output] = output;

  if (measurements->empty()) {
    return;
  }

  const double ca = std::cos(m_cfg.crossingAngle);
  const double sa = std::sin(m_cfg.crossingAngle);

  // ------------------------------------------------------------------
  // Transform fitted measurement positions to global coordinates, rotate them
  // into the ion frame, and group them by station. Constituent TrackerHits are
  // retained so TrackSeed truth relations remain complete.
  // ------------------------------------------------------------------
  struct IonHit {
    double x, y, z;
    std::vector<edm4eic::TrackerHit> constituents;
  };
  std::vector<IonHit> ion;
  ion.reserve(measurements->size());
  std::map<unsigned int, std::vector<std::size_t>> byStation;

  const auto& surfaceMap = m_acts_context->surfaceMap();
  const auto& gctx       = m_acts_context->getActsGeometryContext();

  for (const auto& measurement : *measurements) {
    const auto relatedHits = measurement.getHits();
    const auto weights     = measurement.getWeights();
    if (relatedHits.empty() || weights.size() != relatedHits.size()) {
      warning("Skipping B0 measurement with {} hits and {} weights", relatedHits.size(),
              weights.size());
      continue;
    }

    const auto maxWeight = std::max_element(weights.begin(), weights.end());
    const std::size_t representativeIndex =
        static_cast<std::size_t>(std::distance(weights.begin(), maxWeight));
    const auto representativeHit = relatedHits[representativeIndex];
    const auto* context          = m_converter->findContext(representativeHit.getCellID());
    if (context == nullptr) {
      warning("No DD4hep context for B0 measurement constituent {:#018x}",
              representativeHit.getCellID());
      continue;
    }

    bool mixedSensors = false;
    std::vector<edm4eic::TrackerHit> constituents;
    constituents.reserve(relatedHits.size());
    for (const auto& hit : relatedHits) {
      const auto* hitContext = m_converter->findContext(hit.getCellID());
      if (hitContext == nullptr || hitContext->identifier != context->identifier) {
        mixedSensors = true;
        break;
      }
      constituents.push_back(hit);
    }
    if (mixedSensors) {
      warning("Skipping B0 measurement whose constituent hits span multiple sensors");
      continue;
    }

    const auto surfaceIt = surfaceMap.find(context->identifier);
    if (surfaceIt == surfaceMap.end() ||
        surfaceIt->second->geometryId().value() != measurement.getSurface()) {
      warning("B0 measurement surface does not match constituent sensor {:#018x}",
              context->identifier);
      continue;
    }

    const auto loc = measurement.getLoc();
    const Acts::Vector2 local{loc.a, loc.b};
    const Acts::Vector3 global =
        surfaceIt->second->localToGlobal(gctx, local, Acts::Vector3::Zero());
    if (!global.allFinite()) {
      warning("Skipping B0 measurement with non-finite global position");
      continue;
    }

    ion.push_back({global.x() * ca - global.z() * sa, global.y(), global.x() * sa + global.z() * ca,
                   std::move(constituents)});
  }

  // Group by ion-frame z, not cellID layer. Official B0 uses layer 1-4 for
  // four disks; the realistic geometry uses layer 1-8 (front/back per disk).
  // (layer+1)/2 therefore collapses the official detector to two stations
  // and emits no seeds. Official disks are ~270 mm apart; realistic
  // front/back faces of one disk are ~7 mm.
  if (!ion.empty()) {
    std::vector<std::size_t> order(ion.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return ion[a].z < ion[b].z; });
    unsigned int station = 0;
    double last_z        = ion[order.front()].z;
    for (std::size_t idx : order) {
      if (ion[idx].z - last_z > m_cfg.stationZGap) {
        ++station;
      }
      last_z = ion[idx].z;
      byStation[station].push_back(idx);
    }
  }

  // The parabola fit needs >= 3 points, so 3 stations is the hard floor
  // regardless of configuration (a 2-station mode would need a line fit with
  // a momentum prior instead).
  const unsigned int minStations = std::max(m_cfg.minStations, 3u);
  if (byStation.size() < minStations) {
    return;
  }

  // ------------------------------------------------------------------
  // Enumerate candidates: one hit per station, capped combinatorics.
  // In addition to the full station set, try leave-one-station-out
  // subsets (>= minStations): a proton that scattered hard at station k
  // has a clean subset excluding the post-kink stations, whose fit is
  // unbiased, while the all-station fit averages over the kink.
  // ------------------------------------------------------------------
  std::vector<std::vector<std::size_t>> allStationHits;
  allStationHits.reserve(byStation.size());
  for (const auto& [station, idxs] : byStation) {
    allStationHits.push_back(idxs);
  }

  std::vector<std::vector<std::vector<std::size_t>>> subsets;
  subsets.push_back(allStationHits);
  if (allStationHits.size() > minStations) {
    for (std::size_t drop = 0; drop < allStationHits.size(); ++drop) {
      std::vector<std::vector<std::size_t>> sub;
      for (std::size_t s = 0; s < allStationHits.size(); ++s) {
        if (s != drop) {
          sub.push_back(allStationHits[s]);
        }
      }
      subsets.push_back(sub);
    }
  }

  std::vector<StubCandidate> candidates;
  // Split the combination budget across subsets so a busy event cannot
  // exhaust it inside the full-station set and starve the leave-one-out
  // subsets that exist to recover kinked tracks.
  const unsigned int perSubsetBudget =
      std::max(1u, m_cfg.maxCombinations / static_cast<unsigned int>(subsets.size()));

  auto fitCandidate = [&](const std::vector<std::size_t>& idxs) -> StubCandidate {
    StubCandidate cand;
    cand.hitIndices = idxs;

    const std::size_t n = idxs.size();
    // parabola fit x'(z') via 3x3 normal equations
    Eigen::Matrix3d A  = Eigen::Matrix3d::Zero();
    Eigen::Vector3d bx = Eigen::Vector3d::Zero();
    // line fit y'(z') via 2x2
    double s0 = n, s1 = 0, s2 = 0, sy = 0, syz = 0;

    for (std::size_t idx : idxs) {
      const double z  = ion[idx].z;
      const double x  = ion[idx].x;
      const double y  = ion[idx].y;
      const double z2 = z * z;
      A(0, 0) += 1.0;
      A(0, 1) += z;
      A(0, 2) += z2;
      A(1, 2) += z * z2;
      A(2, 2) += z2 * z2;
      bx(0) += x;
      bx(1) += x * z;
      bx(2) += x * z2;
      s1 += z;
      s2 += z2;
      sy += y;
      syz += y * z;
    }
    A(1, 0) = A(0, 1);
    A(1, 1) = A(0, 2);
    A(2, 0) = A(0, 2);
    A(2, 1) = A(1, 2);

    const Eigen::Vector3d c = A.fullPivLu().solve(bx);
    cand.c0                 = c(0);
    cand.c1                 = c(1);
    cand.c2                 = c(2);

    const double det = s0 * s2 - s1 * s1;
    cand.b1          = (s0 * syz - s1 * sy) / det;
    cand.b0          = (sy - cand.b1 * s1) / s0;

    for (std::size_t idx : idxs) {
      const double z  = ion[idx].z;
      const double rx = ion[idx].x - (cand.c0 + cand.c1 * z + cand.c2 * z * z);
      const double ry = ion[idx].y - (cand.b0 + cand.b1 * z);
      cand.chi2 += rx * rx + ry * ry;
    }
    cand.chi2 /= static_cast<double>(n);
    return cand;
  };

  for (const auto& stationHits : subsets) {
    std::vector<std::size_t> cursor(stationHits.size(), 0);
    unsigned int tried = 0;
    while (tried < perSubsetBudget) {
      std::vector<std::size_t> idxs;
      idxs.reserve(stationHits.size());
      for (std::size_t s = 0; s < stationHits.size(); ++s) {
        idxs.push_back(stationHits[s][cursor[s]]);
      }
      StubCandidate cand = fitCandidate(idxs);
      if (std::isfinite(cand.chi2) && std::isfinite(cand.c1) && std::isfinite(cand.c2) &&
          std::isfinite(cand.b0) && std::isfinite(cand.b1)) {
        candidates.push_back(std::move(cand));
      }
      ++tried;

      // advance the multi-index cursor
      std::size_t s = 0;
      for (; s < stationHits.size(); ++s) {
        if (++cursor[s] < stationHits[s].size()) {
          break;
        }
        cursor[s] = 0;
      }
      if (s == stationHits.size()) {
        break; // enumerated this subset
      }
    }
  }

  if (candidates.empty()) {
    return;
  }

  // ------------------------------------------------------------------
  // Rank by residual, deduplicate by shared hits, emit up to maxSeeds.
  // ------------------------------------------------------------------
  // More stations first (a 3-point parabola interpolates exactly, so raw
  // chi2 would always favour subsets over genuine full-station candidates),
  // then residual within equal station count.
  std::sort(candidates.begin(), candidates.end(),
            [](const StubCandidate& a, const StubCandidate& b) {
              if (a.hitIndices.size() != b.hitIndices.size()) {
                return a.hitIndices.size() > b.hitIndices.size();
              }
              return a.chi2 < b.chi2;
            });

  std::vector<const StubCandidate*> accepted;
  for (const auto& cand : candidates) {
    if (accepted.size() >= m_cfg.maxSeeds) {
      break;
    }
    bool overlaps = false;
    for (const auto* acc : accepted) {
      unsigned int shared = 0;
      for (std::size_t i : cand.hitIndices) {
        for (std::size_t j : acc->hitIndices) {
          if (i == j) {
            ++shared;
          }
        }
      }
      if (shared > m_cfg.maxSharedHits) {
        overlaps = true;
        break;
      }
    }
    if (!overlaps) {
      accepted.push_back(&cand);
    }
  }

  // ------------------------------------------------------------------
  // Convert accepted candidates to origin-perigee seed parameters.
  // ------------------------------------------------------------------
  for (const auto* cand : accepted) {
    // momentum from the parabola curvature: dtx/dz = kappa = q*0.2998e-3*B/p
    // (p in GeV, B in T, z in mm)
    const double kappa = 2.0 * cand->c2;
    double p           = m_cfg.momentumPrior;
    if (std::isfinite(kappa) && std::abs(kappa) > 0.0) {
      const double pFit = 2.998e-4 * m_cfg.bFieldY / std::abs(kappa);
      if (std::isfinite(pFit) && pFit > m_cfg.pMin && pFit < m_cfg.pMax) {
        p = pFit;
      }
    }

    // state at the field entrance (ion frame)
    const double zIn  = m_cfg.zFieldEntrance;
    const double txIn = cand->c1 + 2.0 * cand->c2 * zIn;
    const double xIn  = cand->c0 + cand->c1 * zIn + cand->c2 * zIn * zIn;
    const double yIn  = cand->b0 + cand->b1 * zIn;

    // Direction upstream of the field, in the ion frame.
    //
    // Beamline-constrained (default): the vertex is at the origin, so use the
    // chord origin -> field-entrance point. Position errors of O(0.1 mm) over
    // ~5.8 m give ~20 urad direction accuracy, whereas the fitted slope
    // back-extrapolated over the same lever arm amplifies mm-level fit errors
    // into meter-level perigee positions.
    double txi;
    double tyi;
    double x0;
    double y0;
    if (m_cfg.constrainToBeamline) {
      txi = xIn / zIn;
      tyi = yIn / zIn;
      x0  = 0.0;
      y0  = 0.0;
    } else {
      txi = txIn;
      tyi = cand->b1;
      x0  = xIn - txIn * zIn;
      y0  = cand->b0; // line intercept at z' = 0
    }

    // rotate direction and reference point back to the global frame
    double dx          = txi * ca + 1.0 * sa;
    double dy          = tyi;
    double dz          = -txi * sa + 1.0 * ca;
    const double dnorm = std::sqrt(dx * dx + dy * dy + dz * dz);
    dx /= dnorm;
    dy /= dnorm;
    dz /= dnorm;

    const double px0 = x0 * ca; // ion (x0, y0, 0) -> global
    const double py0 = y0;
    const double pz0 = -x0 * sa;

    // point of closest approach to the z axis
    const double denom = dx * dx + dy * dy;
    const double tPca  = denom > 0.0 ? -(px0 * dx + py0 * dy) / denom : 0.0;
    const double xPca  = px0 + tPca * dx;
    const double yPca  = py0 + tPca * dy;
    const double zPca  = pz0 + tPca * dz;

    const double phi   = std::atan2(dy, dx);
    const double theta = std::acos(std::clamp(dz, -1.0, 1.0));

    // ACTS perigee local frame: loc0 = signed transverse impact, loc1 = z
    const double loc0 = -xPca * std::sin(phi) + yPca * std::cos(phi);
    const double loc1 = zPca;

    auto trackparam = track_params_output->create();
    trackparam.setType(-1); // seed
    trackparam.setLoc({static_cast<float>(loc0), static_cast<float>(loc1)});
    trackparam.setPhi(static_cast<float>(phi));
    trackparam.setTheta(static_cast<float>(theta));
    trackparam.setQOverP(static_cast<float>(m_cfg.charge / p));
    trackparam.setTime(10);

    edm4eic::Cov6f cov;
    cov(0, 0) = m_cfg.locaError;
    cov(1, 1) = m_cfg.locbError;
    cov(2, 2) = m_cfg.phiError;
    cov(3, 3) = m_cfg.thetaError;
    cov(4, 4) = m_cfg.qOverPError;
    cov(5, 5) = m_cfg.timeError;
    trackparam.setCovariance(cov);

    auto seed = seeds->create();
    seed.setPerigee({0.f, 0.f, 0.f});
    // Larger is better for ACTS seed quality; use the negative mean residual.
    seed.setQuality(static_cast<float>(-cand->chi2));
    seed.setParams(trackparam);
    for (std::size_t idx : cand->hitIndices) {
      for (const auto& hit : ion[idx].constituents) {
        seed.addToHits(hit);
      }
    }

    trace("B0 stub seed: p={:.2f} GeV theta={:.4f} phi={:.3f} loc=({:.2f},{:.2f}) chi2={:.4f}", p,
          theta, phi, loc0, loc1, cand->chi2);
  }
}

} // namespace eicrecon
