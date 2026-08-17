// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <array>
#include <memory>
#include <string_view>
#include <vector>

#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"
#include "algorithms/tracking/B0TrackerStubSeederConfig.h"

namespace eicrecon {

/// Geometry-free pieces of the B0 stub seeder, separated out so that the
/// trajectory model can be exercised without a detector or a field service.
namespace b0stub {

  struct Point3 {
    double x{};
    double y{};
    double z{};
    /// Measurement variances along the ion-frame bend/non-bend sensor axes.
    double varianceX{};
    double varianceY{};
  };

  /// Auxiliary fit in the ion frame: a parabola in the bend plane and a line
  /// in the non-bend plane. It supplies compatibility residuals and the
  /// initial field-sampling path; momentum comes from FieldIntegralFit.
  struct StubFit {
    double c0{}; ///< x(z) = c0 + c1 z + c2 z^2
    double c1{};
    double c2{};
    double b0{}; ///< y(z) = b0 + b1 z
    double b1{};
    double rmsX{};
    double rmsY{};
    double zRef{};
    double zScale{};
    /// Coefficients and covariances in the numerically stable basis
    /// u=(z-zRef)/zScale, in row-major order. They are populated when every
    /// input point has positive finite x/y variances.
    std::array<double, 3> basisX{};
    std::array<double, 2> basisY{};
    std::array<double, 9> covarianceX{};
    std::array<double, 4> covarianceY{};
    bool covarianceValid{false};
    bool valid{false};

    double x(double z) const { return c0 + c1 * z + c2 * z * z; }
    double y(double z) const { return b0 + b1 * z; }
    /// bend-plane slope dx/dz
    double tx(double z) const { return c1 + 2.0 * c2 * z; }
  };

  /// Least-squares fit of `pts` (ion frame, mm). Needs at least three points;
  /// z is centred and scaled internally so the equations stay well conditioned.
  StubFit fitStub(const std::vector<Point3>& pts);

  struct FieldSample {
    double z{};      ///< ion-frame position [mm]
    double fieldY{}; ///< dipole field [T]
  };

  struct FieldIntegral {
    double first{};  ///< integral B_y ds [T mm]
    double second{}; ///< integral (z-s) B_y ds [T mm^2]
  };

  /// Integrate a sorted field mesh with a piecewise-linear B_y model. The
  /// returned moments have the same length and are referenced to samples[0].
  std::vector<FieldIntegral> integrateFieldSamples(const std::vector<FieldSample>& samples);

  /// Weighted linear bend-plane fit
  /// x(z)=x_ref+tx_ref*(z-z_ref)-K*(q/p)*I(z).
  struct FieldIntegralFit {
    double zReference{};
    double xReference{};
    double txReference{};
    double qOverP{};
    double rmsX{};
    /// Covariance of (xReference,txReference,qOverP), row-major.
    std::array<double, 9> covariance{};
    bool covarianceValid{false};
    bool valid{false};
  };
  FieldIntegralFit fitFieldIntegral(const std::vector<Point3>& pts,
                                    const std::vector<double>& secondIntegrals, double zReference);

  struct EndpointCompatibility {
    double slopeY{};
    double beamResidual{};
    double score{};
    bool valid{false};
  };
  /// Fast non-bending-plane compatibility for an outer-station hit pair.
  EndpointCompatibility endpointCompatibility(const Point3& first, const Point3& last,
                                              double maxAbsSlope, double maxBeamResidual,
                                              bool constrainToBeamline);

  /// Diagonal (phi, theta, q/p) variance additions for calibrated residual and
  /// external field-integral uncertainties.
  std::array<double, 3> covarianceModelAdditions(double qOverP, double phiModelVariance,
                                                 double phiQOverPScale, double thetaModelVariance,
                                                 double thetaQOverPScale,
                                                 double qOverPModelVariance,
                                                 double qOverPRelativeUncertainty,
                                                 double fieldRelativeUncertainty);

  /// Propagate the fitted coefficient covariance to
  /// (loc0,loc1,phi,theta,q/p), returned as a row-major 5x5 matrix.
  std::array<double, 25> seedCovarianceFromFit(const FieldIntegralFit& bendFit,
                                               const StubFit& nonBendFit, double crossingAngle,
                                               bool constrainToBeamline);

  /// Perigee parameters of the straight ray leaving `ref` along `dir`, expressed
  /// on a perigee surface centred at `perigee`. All positions in mm.
  struct PerigeeParams {
    double loc0{};
    double loc1{};
    double phi{};
    double theta{};
  };
  PerigeeParams perigeeFromRay(const Point3& ref, const Point3& dir, const Point3& perigee);

} // namespace b0stub

using B0TrackerStubSeederAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Dedicated seeder for the B0 tracker (dipole spectrometer at z ~ 6 m).
///
/// Fits one hit per station with a straight line (non-bend plane) and a
/// sampled field-integral basis (bend plane) in the ion-rotated frame, extracts direction,
/// signed q/p (from the sampled dipole field integral), back-extrapolates
/// analytically to the origin, and emits seed parameters on the origin
/// perigee surface that CKFTracking expects.
class B0TrackerStubSeeder : public B0TrackerStubSeederAlgorithm,
                            public WithPodConfig<B0TrackerStubSeederConfig> {
public:
  B0TrackerStubSeeder(std::string_view name)
      : B0TrackerStubSeederAlgorithm{name,
                                     {"inputHits"},
                                     {"outputSeeds", "outputTrackParameters"},
                                     "stub-based track seeds for the B0 dipole tracker"} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  std::shared_ptr<const ActsGeometryProvider> m_acts_context;
  double m_crossing_angle{};
};

} // namespace eicrecon
