// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include "B0TelescopeMath.h"

#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>

namespace eicrecon::b0 {
namespace {
using Transport = Eigen::Matrix<double, 8, 1>;
// Lorentz equation with z as independent coordinate; full slope factors.
Eigen::Vector2d acceleration(double tx, double ty, const Eigen::Vector3d& b) {
  const double scale = bendConstant * std::sqrt(1.0 + tx * tx + ty * ty);
  return scale * Eigen::Vector2d(tx * ty * b.x() - (1.0 + tx * tx) * b.y() + ty * b.z(),
                                (1.0 + ty * ty) * b.x() - tx * ty * b.y() - tx * b.z());
}
Eigen::Vector2d response(const State& ref, double z0, double z, const Field& field,
                         const FitOptions& options) {
  const double dz = z - z0;
  if (!options.fieldIntegral) {
    const auto b = field({ref(0) + 0.5 * dz * ref(2), ref(1) + 0.5 * dz * ref(3), z0 + 0.5 * dz});
    return (0.5 * dz * dz * acceleration(ref(2), ref(3), b)).eval();
  }
  // Integrate the reference trajectory AND the two field response integrals.
  // Do not divide trajectory displacement by q/p: that is singular at q/p=0.
  Transport s = Transport::Zero();
  s.head<4>() = ref.head<4>();
  const double h = dz / options.integrationSteps;
  auto rhs = [&](double at, const Transport& p) -> Transport {
    const auto f = acceleration(p(2), p(3), field({p(0), p(1), at}));
    Transport d;
    d << p(2), p(3), ref(4) * f(0), ref(4) * f(1), p(6), p(7), f(0), f(1);
    return d;
  };
  for (unsigned i = 0; i < options.integrationSteps; ++i) {
    const double at = z0 + i * h;
    const Transport k1 = rhs(at, s);
    const Transport k2 = rhs(at + h / 2, s + h * k1 / 2);
    const Transport k3 = rhs(at + h / 2, s + h * k2 / 2);
    const Transport k4 = rhs(at + h, s + h * k3);
    s += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
  }
  return s.segment<2>(4);
}
Eigen::Vector2d predict(const TelescopeFit& fit, double z, const Field& field,
                         const FitOptions& options) {
  return (fit.state.head<2>() + (z - fit.zReference) * fit.state.segment<2>(2) +
          fit.state(4) * response(fit.state, fit.zReference, z, field, options)).eval();
}
bool goodHit(const TelescopeHit& h) {
  return h.position.allFinite() && h.covariance.allFinite() &&
         h.covariance(0, 0) > 0 && h.covariance(1, 1) > 0 &&
         h.covariance.isApprox(h.covariance.transpose(), 1e-10);
}
std::size_t stationCount(const std::vector<TelescopeHit>& hits) {
  std::set<std::size_t> stations;
  for (const auto& hit : hits) stations.insert(hit.station);
  return stations.size();
}
} // namespace

TelescopeFit fitTelescope(const std::vector<TelescopeHit>& hits, const Field& field,
                          const FitOptions& options) {
  TelescopeFit out;
  if (!field || options.iterations == 0 || options.iterations > 20 ||
      options.integrationSteps == 0 || options.integrationSteps > 1024 || hits.size() < 3)
    return out;
  if (!std::all_of(hits.begin(), hits.end(), goodHit) || stationCount(hits) < 3) return out;
  const auto endpoints = std::minmax_element(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
    return a.position.z() < b.position.z();
  });
  out.zReference = endpoints.first->position.z();
  const double scale = endpoints.second->position.z() - out.zReference;
  if (!(scale > 0)) return out;
  State reference = State::Zero();
  reference.head<2>() = endpoints.first->position.head<2>();
  reference.segment<2>(2) = (endpoints.second->position - endpoints.first->position).head<2>() / scale;
  Eigen::MatrixXd design(2 * hits.size(), 5);
  Eigen::VectorXd values(2 * hits.size());
  Eigen::Matrix<double, 5, 5> unscale = Covariance::Identity();
  unscale(2, 2) = unscale(3, 3) = 1.0 / scale;
  for (unsigned iteration = 0; iteration < options.iterations; ++iteration) {
    for (std::size_t i = 0; i < hits.size(); ++i) {
      const auto& hit = hits[i];
      const double dz = hit.position.z() - out.zReference;
      const auto integral = response(reference, out.zReference, hit.position.z(), field, options);
      Eigen::Matrix<double, 2, 5> row;
      row << 1, 0, dz / scale, 0, integral.x(), 0, 1, 0, dz / scale, integral.y();
      Eigen::Matrix<double, 2, 3> projection;
      projection << 1, 0, -reference(2), 0, 1, -reference(3);
      const Eigen::Matrix2d residualCov = projection * hit.covariance * projection.transpose();
      const Eigen::LLT<Eigen::Matrix2d> llt(residualCov);
      if (llt.info() != Eigen::Success) return out;
      design.middleRows<2>(2 * i) = llt.matrixL().solve(row);
      values.segment<2>(2 * i) = llt.matrixL().solve(hit.position.head<2>());
    }
    if (!design.allFinite()) return out;
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
    if (qr.rank() != 5) return out; // no measured curvature in zero field
    const State scaled = qr.solve(values);
    reference = unscale * scaled;
    if (!reference.allFinite()) return out;
    const Eigen::MatrixXd inverse = qr.solve(Eigen::MatrixXd::Identity(design.rows(), design.rows()));
    out.covariance = unscale * inverse * inverse.transpose() * unscale.transpose();
    out.chi2 = (design * scaled - values).squaredNorm();
  }
  out.state = reference;
  out.ndf = static_cast<int>(2 * hits.size()) - 5;
  out.covariance = (0.5 * (out.covariance + out.covariance.transpose())).eval();
  out.valid = out.covariance.allFinite() && out.covariance.llt().info() == Eigen::Success &&
              std::isfinite(out.chi2);
  return out;
}

std::vector<TelescopeFit> chargeProposals(const TelescopeFit& fit, double significance,
                                         double maxMomentum, int forceCharge, bool bothCharges) {
  if (!(significance >= 0) || !std::isfinite(significance) || !(maxMomentum > 0) ||
      !std::isfinite(maxMomentum) || (forceCharge != 0 && forceCharge != -1 && forceCharge != 1))
    throw std::invalid_argument("Invalid B0 charge-proposal configuration");
  if (!fit.valid || !fit.state.allFinite() || !fit.covariance.allFinite() ||
      !(fit.covariance(4, 4) > 0)) return {};
  const double q = fit.state(4), sigma = std::sqrt(fit.covariance(4, 4));
  const bool unresolved = q == 0.0 || std::abs(q) / sigma < significance;
  std::vector<int> charges = bothCharges ? std::vector<int>{-1, 1} :
      forceCharge != 0 ? std::vector<int>{forceCharge} :
      unresolved ? std::vector<int>{-1, 1} : std::vector<int>{q < 0 ? -1 : 1};
  std::vector<TelescopeFit> result;
  for (int charge : charges) {
    auto proposal = fit;
    const double magnitude = unresolved ? std::max({std::abs(q), sigma, 1.0 / maxMomentum}) :
                                         std::max(std::abs(q), 1.0 / maxMomentum);
    proposal.state(4) = charge * magnitude;
    proposal.state.head<4>() += fit.covariance.block<4, 1>(0, 4) / fit.covariance(4, 4) *
                                (proposal.state(4) - q);
    result.push_back(proposal);
  }
  return result;
}

bool timeCompatible(const TelescopeHit& a, const TelescopeHit& b, double qOverP,
                    double mass, double nSigma, double modelSigmaNs) {
  if (!(nSigma > 0) || !(mass >= 0) || !(modelSigmaNs >= 0) || !std::isfinite(qOverP) ||
      !std::isfinite(mass) || !std::isfinite(nSigma) || !std::isfinite(modelSigmaNs))
    throw std::invalid_argument("Invalid B0 timing configuration");
  if (!std::isfinite(a.time) || !std::isfinite(b.time) || !(a.timeVariance > 0) ||
      !(b.timeVariance > 0) || !std::isfinite(a.timeVariance + b.timeVariance)) return true;
  const double inverseBeta = std::sqrt(1.0 + mass * mass * qOverP * qOverP); // |charge|=1
  const double distance = (b.position - a.position).norm();
  const double expected = std::copysign(distance, b.position.z() - a.position.z()) * inverseBeta / speedOfLight;
  const double sigma = std::sqrt(a.timeVariance + b.timeVariance + modelSigmaNs * modelSigmaNs);
  return std::abs((b.time - a.time) - expected) <= nSigma * sigma;
}

SearchResult findTelescopeCandidates(std::vector<TelescopeHit> hits, const Field& field,
                                     const SearchOptions& cfg) {
  if (cfg.minStations < 3 || cfg.maxTrials == 0 || cfg.maxCandidates == 0 ||
      !(cfg.minMomentum > 0) || !(cfg.maxSlope > 0) || !(cfg.fieldBoundTesla > 0) ||
      !(cfg.roadWidthMm > 0) || !(cfg.maxChi2PerDof > 0) || !(cfg.extensionChi2 > 0) ||
      cfg.fit.iterations == 0 || cfg.fit.iterations > 20 || cfg.fit.integrationSteps == 0 ||
      cfg.fit.integrationSteps > 1024 || !field ||
      !std::isfinite(cfg.minMomentum) || !std::isfinite(cfg.maxSlope) ||
      !std::isfinite(cfg.fieldBoundTesla) || !std::isfinite(cfg.roadWidthMm) ||
      !std::isfinite(cfg.maxChi2PerDof) || !std::isfinite(cfg.extensionChi2) ||
      !(cfg.timingNSigma > 0) || !std::isfinite(cfg.timingNSigma) ||
      !(cfg.timingModelSigmaNs >= 0) || !std::isfinite(cfg.timingModelSigmaNs) ||
      !(cfg.massGeV >= 0) || !std::isfinite(cfg.massGeV))
    throw std::invalid_argument("Invalid B0 candidate-search configuration");
  hits.erase(std::remove_if(hits.begin(), hits.end(), [](const auto& h) { return !goodHit(h); }), hits.end());
  std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
    return std::tie(a.station, a.position.y(), a.position.x(), a.surface, a.index) <
           std::tie(b.station, b.position.y(), b.position.x(), b.surface, b.index);
  });
  std::map<std::size_t, std::vector<const TelescopeHit*>> groups;
  for (const auto& hit : hits) groups[hit.station].push_back(&hit);
  std::vector<std::vector<const TelescopeHit*>> stations;
  for (const auto& [id, group] : groups) { (void)id; stations.push_back(group); }
  SearchResult result;
  if (stations.size() < cfg.minStations) return result;
  auto spend = [&]() {
    if (result.trials >= cfg.maxTrials) { result.truncated = true; return false; }
    ++result.trials;
    return true;
  };
  auto admissible = [&](const TelescopeFit& f) {
    return f.valid && std::abs(f.state(4)) <= 1.0 / cfg.minMomentum &&
           f.state.segment<2>(2).cwiseAbs().maxCoeff() <= cfg.maxSlope &&
           f.chi2 <= cfg.maxChi2PerDof * std::max(1, f.ndf);
  };
  auto rank = [](const Candidate& a, const Candidate& b) {
    if (a.stations != b.stations) return a.stations > b.stations;
    if (a.hits.size() != b.hits.size()) return a.hits.size() > b.hits.size();
    if (a.fit.chi2 != b.fit.chi2) return a.fit.chi2 < b.fit.chi2;
    return a.hits < b.hits;
  };
  std::set<std::vector<std::size_t>> unique;
  auto retain = [&](std::vector<TelescopeHit> selected, const TelescopeFit& fit) {
    const auto nStations = stationCount(selected);
    if (nStations < cfg.minStations) return;
    std::sort(selected.begin(), selected.end(), [](const auto& a, const auto& b) {
      return std::tie(a.position.z(), a.surface, a.index) < std::tie(b.position.z(), b.surface, b.index);
    });
    Candidate c; c.fit = fit; c.stations = nStations;
    for (const auto& h : selected) c.hits.push_back(h.index);
    auto key = c.hits; std::sort(key.begin(), key.end());
    if (!unique.insert(key).second) return;
    result.candidates.push_back(std::move(c));
    std::sort(result.candidates.begin(), result.candidates.end(), rank);
    if (result.candidates.size() > cfg.maxCandidates) {
      auto dropped = result.candidates.back().hits;
      std::sort(dropped.begin(), dropped.end()); unique.erase(dropped);
      result.candidates.pop_back(); ++result.candidatesDropped;
    }
  };
  const double slope2 = cfg.maxSlope * cfg.maxSlope;
  const double curvature = bendConstant / cfg.minMomentum * cfg.fieldBoundTesla *
                           std::sqrt(1 + 2 * slope2) * (1 + 2 * slope2 + cfg.maxSlope);
  for (std::size_t a = 0; a + 2 < stations.size() && !result.truncated; ++a)
    for (std::size_t c = a + 2; c < stations.size() && !result.truncated; ++c)
      for (std::size_t b = a + 1; b < c && !result.truncated; ++b)
        for (const auto* first : stations[a]) {
          if (result.truncated) break;
          for (const auto* last : stations[c]) {
            if (!spend()) break; // endpoint enumeration itself is bounded
            const double dz = last->position.z() - first->position.z();
            if (!(dz > 0)) continue;
            const Eigen::Vector2d slope = (last->position - first->position).head<2>() / dz;
            if (slope.cwiseAbs().maxCoeff() > cfg.maxSlope + curvature * dz) continue;
            for (const auto* middle : stations[b]) {
              if (!spend()) break;
              const double z = middle->position.z() - first->position.z();
              const Eigen::Vector2d residual = middle->position.head<2>() - first->position.head<2>() - z * slope;
              if (residual.cwiseAbs().maxCoeff() > cfg.roadWidthMm + 0.5 * curvature * z * (dz - z)) continue;
              std::vector<TelescopeHit> selected{*first, *middle, *last};
              ++result.fits;
              auto fit = fitTelescope(selected, field, cfg.fit);
              if (!admissible(fit)) continue;
              if (cfg.useTiming && (!timeCompatible(*first, *middle, fit.state(4), cfg.massGeV, cfg.timingNSigma, cfg.timingModelSigmaNs) ||
                                    !timeCompatible(*middle, *last, fit.state(4), cfg.massGeV, cfg.timingNSigma, cfg.timingModelSigmaNs))) {
                ++result.timingRejected; continue;
              }
              retain(selected, fit); // keep genuine three-station recovery candidates
              std::set<std::uint64_t> used{first->surface, middle->surface, last->surface};
              std::map<std::uint64_t, std::pair<double, const TelescopeHit*>> extensions;
              for (const auto& h : hits) {
                if (!spend()) break;
                if (used.count(h.surface)) continue;
                const auto d = (h.position.head<2>() - predict(fit, h.position.z(), field, cfg.fit)).eval();
                // A deliberately conservative search covariance; final fitting uses
                // the original measurement covariance, not this road allowance.
                Eigen::Matrix2d cov = h.covariance.topLeftCorner<2, 2>();
                cov.diagonal().array() += cfg.roadWidthMm * cfg.roadWidthMm / cfg.extensionChi2;
                const double chi2 = d.dot(cov.ldlt().solve(d));
                if (chi2 > cfg.extensionChi2) continue;
                if (cfg.useTiming && !timeCompatible(*first, h, fit.state(4), cfg.massGeV, cfg.timingNSigma, cfg.timingModelSigmaNs)) {
                  ++result.timingRejected; continue;
                }
                const auto old = extensions.find(h.surface);
                if (old == extensions.end() || chi2 < old->second.first) extensions[h.surface] = {chi2, &h};
              }
              if (result.truncated) break; // never emit a silently partial extension
              for (const auto& [surface, choice] : extensions) { (void)surface; selected.push_back(*choice.second); }
              if (!extensions.empty()) {
                ++result.fits;
                fit = fitTelescope(selected, field, cfg.fit);
                if (admissible(fit)) retain(selected, fit);
              }
            }
          }
        }
  return result;
}
} // namespace eicrecon::b0
