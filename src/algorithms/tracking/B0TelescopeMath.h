// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace eicrecon::b0 {
using State = Eigen::Matrix<double, 5, 1>; // x, y [mm], dx/dz, dy/dz, q/p [GeV^-1]
using Covariance = Eigen::Matrix<double, 5, 5>;
using Field = std::function<Eigen::Vector3d(const Eigen::Vector3d&)>; // ion frame, tesla
inline constexpr double bendConstant = 2.99792458e-4; // mm, tesla, GeV
inline constexpr double speedOfLight = 299.792458;   // mm/ns

struct TelescopeHit {
  std::size_t index{}; // original Measurement2D index, never a proximity surrogate
  std::uint64_t surface{};
  std::size_t station{};
  Eigen::Vector3d position = Eigen::Vector3d::Zero();
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  double time = std::numeric_limits<double>::quiet_NaN();
  double timeVariance = 0.0;
};
struct TelescopeFit {
  State state = State::Zero();
  Covariance covariance = Covariance::Zero();
  double zReference{};
  double chi2{};
  int ndf{};
  bool valid{false};
};
struct FitOptions {
  bool fieldIntegral{false};
  unsigned iterations{3};
  unsigned integrationSteps{16};
};

/// A small-angle, five-parameter seed fit, not a final material-aware fit.
/// Field-integral mode iterates a frozen-path linearization using RK4 transport.
/// Covariance includes hit correlations and tilted-plane residual projection,
/// but not material, field-map uncertainty or nonlinear linearization error.
TelescopeFit fitTelescope(const std::vector<TelescopeHit>& hits, const Field& field,
                          const FitOptions& options = {});

/// Finite, self-consistent proposals for an unresolved curvature. Each changed
/// q/p also updates the correlated position/slopes by Gaussian conditioning.
/// This does not turn a charge hypothesis into a separate measurement.
std::vector<TelescopeFit> chargeProposals(const TelescopeFit& fit, double significance = 3.0,
                                         double maxMomentum = 10000.0, int forceCharge = 0,
                                         bool bothCharges = false);

bool timeCompatible(const TelescopeHit& a, const TelescopeHit& b, double qOverP,
                    double mass, double nSigma, double modelSigmaNs);

struct SearchOptions {
  std::size_t minStations{3};
  std::size_t maxTrials{4096}; // includes endpoint tests AND complete candidate tests
  std::size_t maxCandidates{64};
  double minMomentum{3.0};
  double maxSlope{0.15};
  double fieldBoundTesla{5.0}; // explicit bound for preselection, not a measured field
  double roadWidthMm{5.0};
  double maxChi2PerDof{50.0};
  double extensionChi2{50.0};
  bool useTiming{false};
  double timingNSigma{5.0};
  double timingModelSigmaNs{0.1};
  double massGeV{0.93827208816};
  FitOptions fit;
};
struct Candidate {
  std::vector<std::size_t> hits; // original measurement indices, ordered in z
  TelescopeFit fit;
  std::size_t stations{};
};
struct SearchResult {
  std::vector<Candidate> candidates;
  std::size_t trials{};
  std::size_t fits{};
  std::size_t timingRejected{};
  std::size_t candidatesDropped{};
  bool truncated{false};
};
/// Ordered station triplets, bounded work, both-projection roads, then sensor-
/// aware extension. One measurement per actual surface; single-sided stations
/// remain valid. Exact hit-set deduplication never merges close physical tracks.
SearchResult findTelescopeCandidates(std::vector<TelescopeHit> hits, const Field& field,
                                     const SearchOptions& options = {});
} // namespace eicrecon::b0
