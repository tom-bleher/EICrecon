// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
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
  };

  /// Trajectory of one stub in the ion frame: a parabola in the bend plane and a
  /// straight line in the non-bend plane.
  struct StubFit {
    double c0{}; ///< x(z) = c0 + c1 z + c2 z^2
    double c1{};
    double c2{};
    double b0{}; ///< y(z) = b0 + b1 z
    double b1{};
    double rmsX{};
    double rmsY{};
    bool valid{false};

    double x(double z) const { return c0 + c1 * z + c2 * z * z; }
    double y(double z) const { return b0 + b1 * z; }
    /// bend-plane slope dx/dz
    double tx(double z) const { return c1 + 2.0 * c2 * z; }
  };

  /// Least-squares fit of `pts` (ion frame, mm). Needs at least three points for
  /// the parabola; z is centred and scaled internally so that the normal
  /// equations stay well conditioned six metres from the origin.
  StubFit fitStub(const std::vector<Point3>& pts);

  /// Momentum [GeV] implied by the bend-plane curvature in a field `fieldY` [T],
  /// for positions in mm. Returns 0 for a straight (infinite momentum) fit.
  double momentumFromCurvature(double c2, double fieldY);

  /// Charge sign implied by that curvature for a track travelling towards +z.
  /// A positive particle in a +y field bends towards -x, i.e. d(tx)/dz < 0.
  int chargeFromCurvature(double c2, double fieldY);

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
/// parabola (bend plane) in the ion-rotated frame, extracts direction,
/// momentum (from the dipole sagitta) and charge, back-extrapolates
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
};

} // namespace eicrecon
