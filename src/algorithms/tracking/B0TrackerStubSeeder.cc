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
// The seed is only a starting point: the CKF refits with the full field and
// material, so the stub momentum must not be read as a measurement. Because
// the seed is expressed at the origin it is a prompt-track solution; see
// constrainToBeamline in the config for why displaced decays need more.

#include "B0TrackerStubSeeder.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Eigen/Dense>
#include <edm4eic/Cov6f.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "algorithms/interfaces/ActsSvc.h"

namespace eicrecon {

namespace b0stub {

  StubFit fitStub(const std::vector<Point3>& pts) {
    StubFit fit;
    const std::size_t n = pts.size();
    if (n < 3) {
      return fit;
    }

    // Centre and scale z before fitting. B0 measurements sit near z = 6 m;
    // fitting powers through z^4 in normal equations is needlessly ill
    // conditioned. QR on u = (z - zRef)/zScale gives the same polynomial with
    // stable coefficients.
    double zRef = 0.0;
    for (const auto& p : pts) {
      zRef += p.z;
    }
    zRef /= static_cast<double>(n);
    double zScale = 0.0;
    for (const auto& p : pts) {
      zScale = std::max(zScale, std::abs(p.z - zRef));
    }
    if (!(zScale > 0.0)) {
      return fit;
    }

    Eigen::MatrixXd designX(n, 3);
    Eigen::MatrixXd designY(n, 2);
    Eigen::VectorXd valuesX(n);
    Eigen::VectorXd valuesY(n);
    for (std::size_t row = 0; row < n; ++row) {
      const double u  = (pts[row].z - zRef) / zScale;
      designX(row, 0) = 1.0;
      designX(row, 1) = u;
      designX(row, 2) = u * u;
      valuesX(row)    = pts[row].x;
      designY(row, 0) = 1.0;
      designY(row, 1) = u;
      valuesY(row)    = pts[row].y;
    }
    const Eigen::Vector3d ax = designX.colPivHouseholderQr().solve(valuesX);
    const Eigen::Vector2d ay = designY.colPivHouseholderQr().solve(valuesY);

    // Convert back from the shifted/scaled variable to plain powers of z.
    fit.c2 = ax(2) / (zScale * zScale);
    fit.c1 = ax(1) / zScale - 2.0 * zRef * fit.c2;
    fit.c0 = ax(0) - ax(1) * zRef / zScale + ax(2) * zRef * zRef / (zScale * zScale);
    fit.b1 = ay(1) / zScale;
    fit.b0 = ay(0) - ay(1) * zRef / zScale;

    double sumRx2 = 0.0;
    double sumRy2 = 0.0;
    for (const auto& p : pts) {
      const double rx = p.x - fit.x(p.z);
      const double ry = p.y - fit.y(p.z);
      sumRx2 += rx * rx;
      sumRy2 += ry * ry;
    }
    fit.rmsX = std::sqrt(sumRx2 / static_cast<double>(n));
    fit.rmsY = std::sqrt(sumRy2 / static_cast<double>(n));

    fit.valid = std::isfinite(fit.c0) && std::isfinite(fit.c1) && std::isfinite(fit.c2) &&
                std::isfinite(fit.b0) && std::isfinite(fit.b1) && std::isfinite(fit.rmsX) &&
                std::isfinite(fit.rmsY);
    return fit;
  }

  double momentumFromCurvature(double c2, double fieldY) {
    // For a track travelling in +z, d(tx)/dz = -q * 0.2998e-3 * By / p, with z
    // in mm, By in T and p in GeV.
    const double kappa = 2.0 * c2;
    if (!std::isfinite(kappa) || std::abs(kappa) < 1.0e-12) {
      return 0.0;
    }
    return 2.998e-4 * std::abs(fieldY) / std::abs(kappa);
  }

  int chargeFromCurvature(double c2, double fieldY) { return 2.0 * c2 * fieldY > 0.0 ? -1 : 1; }

  PerigeeParams perigeeFromRay(const Point3& ref, const Point3& dir, const Point3& perigee) {
    PerigeeParams out;

    const double dnorm = std::sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
    if (!(dnorm > 0.0)) {
      return out;
    }
    const double dx = dir.x / dnorm;
    const double dy = dir.y / dnorm;
    const double dz = dir.z / dnorm;

    // Point of closest approach to the perigee line (parallel to z through
    // `perigee`), measured from `ref`.
    const double rx    = ref.x - perigee.x;
    const double ry    = ref.y - perigee.y;
    const double denom = dx * dx + dy * dy;
    const double t     = denom > 0.0 ? -(rx * dx + ry * dy) / denom : 0.0;

    const double xPca = rx + t * dx;
    const double yPca = ry + t * dy;
    const double zPca = ref.z + t * dz - perigee.z;

    out.phi   = std::atan2(dy, dx);
    out.theta = std::acos(std::clamp(dz, -1.0, 1.0));
    // ACTS perigee local frame: loc0 = signed transverse impact, loc1 = z
    out.loc0 = -xPca * std::sin(out.phi) + yPca * std::cos(out.phi);
    out.loc1 = zPca;
    return out;
  }

} // namespace b0stub

namespace {

  struct StubCandidate {
    std::vector<std::size_t> hitIndices;
    b0stub::StubFit fit;
    double fieldY{0.0};
    double momentum{0.0};
    int inferredCharge{0};
  };

} // namespace

void B0TrackerStubSeeder::init() {
  m_acts_context = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (m_acts_context == nullptr) {
    throw std::runtime_error("B0TrackerStubSeeder: ACTS geometry service is unavailable");
  }
}

void B0TrackerStubSeeder::process(const Input& input, const Output& output) const {
  const auto [hits]                 = input;
  auto [seeds, track_params_output] = output;

  if (hits->empty()) {
    return;
  }

  const double ca = std::cos(m_cfg.crossingAngle);
  const double sa = std::sin(m_cfg.crossingAngle);

  // ------------------------------------------------------------------
  // Rotate hit positions into the ion frame and group them by station.
  // ------------------------------------------------------------------
  std::vector<b0stub::Point3> ion;
  ion.reserve(hits->size());
  std::vector<edm4eic::TrackerHit> ionHits;
  ionHits.reserve(hits->size());

  for (const auto& hit : *hits) {
    const auto& p = hit.getPosition();
    const b0stub::Point3 q{p.x * ca - p.z * sa, p.y, p.x * sa + p.z * ca};
    if (!(std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z))) {
      debug("Skipping B0 hit with non-finite position");
      continue;
    }
    ion.push_back(q);
    ionHits.push_back(hit);
  }
  if (ion.empty()) {
    return;
  }

  // Group by ion-frame z, not cellID layer. Official B0 uses layer 1-4 for
  // four disks; the realistic geometry uses layer 1-8 (front/back per disk).
  // (layer+1)/2 therefore collapses the official detector to two stations
  // and emits no seeds. Official disks are ~270 mm apart; realistic
  // front/back faces of one disk are ~7 mm.
  std::map<unsigned int, std::vector<std::size_t>> byStation;
  {
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

  // Sample the same ACTS/DD4hep field provider used by CKFTracking. The B0
  // field has a dipole component and a gradient, so an analytic nominal B is
  // not portable across beam-energy configurations or hit positions.
  const auto fieldProvider = m_acts_context->getFieldProvider();
  auto fieldCache = fieldProvider->makeCache(m_acts_context->getActsMagneticFieldContext());

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

  auto makeCompatible = [&](StubCandidate& cand) -> bool {
    if (!cand.fit.valid || cand.fit.rmsX > m_cfg.maxXResidual ||
        cand.fit.rmsY > m_cfg.maxYResidual) {
      return false;
    }

    const auto zMinMax =
        std::minmax_element(cand.hitIndices.begin(), cand.hitIndices.end(),
                            [&](std::size_t a, std::size_t b) { return ion[a].z < ion[b].z; });
    const double zFirst = ion[*zMinMax.first].z;
    const double zLast  = ion[*zMinMax.second].z;
    if (std::abs(cand.fit.tx(zFirst)) > m_cfg.maxAbsTransverseSlope ||
        std::abs(cand.fit.b1) > m_cfg.maxAbsTransverseSlope) {
      return false;
    }

    const unsigned int nSamples = std::max(1u, m_cfg.fieldSamples);
    double sumBy                = 0.0;
    unsigned int nField         = 0;
    for (unsigned int i = 0; i < nSamples; ++i) {
      const double f = nSamples == 1 ? 0.5 : static_cast<double>(i) / (nSamples - 1);
      const double z = zFirst + f * (zLast - zFirst);
      const double x = cand.fit.x(z);
      const double y = cand.fit.y(z);
      const Acts::Vector3 global{x * ca + z * sa, y, -x * sa + z * ca};
      const auto field = fieldProvider->getField(global, fieldCache);
      if (!field.ok() || !field.value().allFinite()) {
        continue;
      }
      sumBy += field.value().y() / Acts::UnitConstants::T;
      ++nField;
    }
    if (nField == 0) {
      return false;
    }
    cand.fieldY = sumBy / static_cast<double>(nField);
    if (!std::isfinite(cand.fieldY) || std::abs(cand.fieldY) < m_cfg.minAbsFieldY) {
      return false;
    }

    cand.momentum = b0stub::momentumFromCurvature(cand.fit.c2, cand.fieldY);
    if (!std::isfinite(cand.momentum) || cand.momentum < m_cfg.pMin || cand.momentum > m_cfg.pMax) {
      return false;
    }
    cand.inferredCharge = b0stub::chargeFromCurvature(cand.fit.c2, cand.fieldY);
    return true;
  };

  unsigned int truncated = 0;
  for (const auto& stationHits : subsets) {
    // The dipole does not bend in y, so the hits of a single track lie on a
    // straight line in y(z). Anchor on the outermost two stations and keep
    // only the intermediate hits that sit inside a road around that chord.
    // This removes most cross-track combinations before any fit, so the
    // budget below binds only in genuinely dense events rather than
    // truncating the enumeration in lexicographic order.
    const std::size_t nStations = stationHits.size();
    const auto& firstStation    = stationHits.front();
    const auto& lastStation     = stationHits.back();

    unsigned int tried = 0;
    bool budgetHit     = false;
    for (std::size_t iFirst = 0; iFirst < firstStation.size() && !budgetHit; ++iFirst) {
      for (std::size_t iLast = 0; iLast < lastStation.size() && !budgetHit; ++iLast) {
        const std::size_t hFirst = firstStation[iFirst];
        const std::size_t hLast  = lastStation[iLast];
        const double zA          = ion[hFirst].z;
        const double zB          = ion[hLast].z;
        if (!(zB > zA)) {
          continue;
        }
        const double slopeY = (ion[hLast].y - ion[hFirst].y) / (zB - zA);

        // Candidate hits per intermediate station, inside the y road.
        std::vector<std::vector<std::size_t>> middle;
        middle.reserve(nStations - 2);
        bool empty = false;
        for (std::size_t s = 1; s + 1 < nStations; ++s) {
          std::vector<std::size_t> keep;
          for (std::size_t idx : stationHits[s]) {
            const double yPred = ion[hFirst].y + slopeY * (ion[idx].z - zA);
            if (std::abs(ion[idx].y - yPred) <= m_cfg.yRoadWidth) {
              keep.push_back(idx);
            }
          }
          if (keep.empty()) {
            empty = true;
            break;
          }
          middle.push_back(std::move(keep));
        }
        if (empty) {
          continue;
        }

        // Enumerate the surviving middle combinations.
        std::vector<std::size_t> cursor(middle.size(), 0);
        while (true) {
          if (tried >= perSubsetBudget) {
            budgetHit = true;
            ++truncated;
            break;
          }
          std::vector<std::size_t> idxs;
          idxs.reserve(nStations);
          idxs.push_back(hFirst);
          for (std::size_t s = 0; s < middle.size(); ++s) {
            idxs.push_back(middle[s][cursor[s]]);
          }
          idxs.push_back(hLast);

          std::vector<b0stub::Point3> pts;
          pts.reserve(idxs.size());
          for (std::size_t idx : idxs) {
            pts.push_back(ion[idx]);
          }
          StubCandidate cand;
          cand.hitIndices = std::move(idxs);
          cand.fit        = b0stub::fitStub(pts);
          if (makeCompatible(cand)) {
            candidates.push_back(std::move(cand));
          }
          ++tried;

          std::size_t s = 0;
          for (; s < middle.size(); ++s) {
            if (++cursor[s] < middle[s].size()) {
              break;
            }
            cursor[s] = 0;
          }
          if (s == middle.size()) {
            break; // enumerated this pair
          }
        }
      }
    }
  }
  if (truncated > 0) {
    debug("B0 stub seeding hit the combination budget in {} of {} station subsets; some "
          "combinations were not tried",
          truncated, subsets.size());
  }

  if (candidates.empty()) {
    return;
  }

  // ------------------------------------------------------------------
  // Rank by residual, deduplicate by shared hits, emit up to maxSeeds.
  // ------------------------------------------------------------------
  // More stations first (a 3-point parabola interpolates exactly, so raw
  // residuals would always favour subsets over genuine full-station
  // candidates), then residual within equal station count.
  std::sort(candidates.begin(), candidates.end(),
            [](const StubCandidate& a, const StubCandidate& b) {
              if (a.hitIndices.size() != b.hitIndices.size()) {
                return a.hitIndices.size() > b.hitIndices.size();
              }
              return a.fit.rmsX * a.fit.rmsX + a.fit.rmsY * a.fit.rmsY <
                     b.fit.rmsX * b.fit.rmsX + b.fit.rmsY * b.fit.rmsY;
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
  unsigned int emitted = 0;
  for (const auto* cand : accepted) {
    if (emitted >= m_cfg.maxSeeds) {
      break;
    }
    const double p = cand->momentum;

    // state at the field entrance (ion frame)
    const double zIn  = m_cfg.zFieldEntrance;
    const double txIn = cand->fit.tx(zIn);
    const double xIn  = cand->fit.x(zIn);
    const double yIn  = cand->fit.y(zIn);

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
      tyi = cand->fit.b1;
      x0  = xIn - txIn * zIn;
      y0  = cand->fit.b0; // line intercept at z' = 0
    }

    // rotate direction and reference point back to the global frame
    const b0stub::Point3 dir{txi * ca + sa, tyi, -txi * sa + ca};
    const b0stub::Point3 ref{x0 * ca, y0, -x0 * sa};
    const auto perigee = b0stub::perigeeFromRay(ref, dir, {0.0, 0.0, 0.0});

    std::vector<int> charges;
    if (m_cfg.testBothCharges) {
      charges = {-1, 1};
    } else if (m_cfg.charge == -1 || m_cfg.charge == 1) {
      charges = {m_cfg.charge};
    } else {
      charges = {cand->inferredCharge};
    }

    for (const int charge : charges) {
      if (emitted >= m_cfg.maxSeeds) {
        break;
      }
      auto trackparam = track_params_output->create();
      trackparam.setType(-1); // seed
      trackparam.setLoc({static_cast<float>(perigee.loc0), static_cast<float>(perigee.loc1)});
      trackparam.setPhi(static_cast<float>(perigee.phi));
      trackparam.setTheta(static_cast<float>(perigee.theta));
      trackparam.setQOverP(static_cast<float>(charge / p));
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
      seed.setQuality(
          static_cast<float>(-(cand->fit.rmsX * cand->fit.rmsX + cand->fit.rmsY * cand->fit.rmsY)));
      seed.setParams(trackparam);
      for (std::size_t idx : cand->hitIndices) {
        seed.addToHits(ionHits[idx]);
      }
      ++emitted;

      trace("B0 stub seed: q={} p={:.2f} GeV By={:.3f} T theta={:.4f} phi={:.3f} "
            "loc=({:.2f},{:.2f}) rms=({:.3f},{:.3f})",
            charge, p, cand->fieldY, perigee.theta, perigee.phi, perigee.loc0, perigee.loc1,
            cand->fit.rmsX, cand->fit.rmsY);
    }
  }
}

} // namespace eicrecon
