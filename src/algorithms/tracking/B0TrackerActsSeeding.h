// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <algorithms/algorithm.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"

namespace eicrecon {

namespace b0acts {
  using Parameters = Eigen::Matrix<double, 5, 1>;
  using Covariance = Eigen::Matrix<double, 5, 5>;

  /// Independent measured coordinates on a fixed detector surface, in ACTS units.
  struct Measurement {
    const Acts::Surface* surface{};
    Acts::Vector2 local            = Acts::Vector2::Zero();
    Acts::SquareMatrix2 covariance = Acts::SquareMatrix2::Zero();
  };

  struct Estimate {
    Parameters parameters = Parameters::Zero();
    Covariance covariance = Covariance::Zero();
  };

  /// The state is bound to the first measurement surface. The field is held
  /// fixed across the candidate and during numerical covariance propagation.
  std::optional<Parameters> estimateParameters(const std::vector<Measurement>& measurements,
                                               const Acts::GeometryContext& gctx,
                                               const Acts::Vector3& field);
  std::optional<Estimate> estimateWithCovariance(const std::vector<Measurement>& measurements,
                                                 const Acts::GeometryContext& gctx,
                                                 const Acts::Vector3& field);
} // namespace b0acts

struct B0TrackerActsSeedingConfig {
  /// A seed already uses the measurements that CKF will update with again.
  /// Inflate the measurement/model covariance to weaken that correlated prior.
  double covarianceInflation = 100.0;
  /// Provisional first-surface scattering/model additions; not a calibration.
  double scatteringScale           = 3.0e-4; // GeV rad; angular sigma = scale * |q/p|
  double qOverPRelativeUncertainty = 0.025;
  double timeVariance              = 100.0; // ns^2
  double minAbsField               = 0.05;  // T
  double minCurvatureSignificance  = 2.0;
  int charge                       = 0;
  bool testBothCharges             = false;
  unsigned int maxSeeds            = 20;
};

using B0TrackerActsSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackSeedCollection, edm4eic::Measurement2DCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Re-estimate the B0 stub candidates with all their assigned sensor measurements.
class B0TrackerActsSeeding : public B0TrackerActsSeedingAlgorithm,
                             public WithPodConfig<B0TrackerActsSeedingConfig> {
public:
  explicit B0TrackerActsSeeding(std::string_view name)
      : B0TrackerActsSeedingAlgorithm{name,
                                      {"candidates", "measurements"},
                                      {"seeds", "parameters"},
                                      "ACTS multipoint B0 seed parameters at the first sensor"} {}
  void init() final;
  void process(const Input& input, const Output& output) const final;

private:
  std::shared_ptr<const ActsGeometryProvider> m_geometry;
  double m_crossingAngle{};
};
} // namespace eicrecon
