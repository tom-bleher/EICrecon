// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string_view>
#include <unordered_map>
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

  /// Least-squares fit in the ion frame: a parabola in the bend plane and a
  /// line in the non-bend plane. In the uniform dipole field of the B0pf the
  /// parabola is the small-angle trajectory, so its curvature is the momentum
  /// measurement (see bendFitFromParabola).
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
    /// Response of basisY to the q/p used to remove non-bend curvature.
    /// covarianceY remains the conditional covariance at fixed q/p.
    std::array<double, 2> derivativeYQOverP{};
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

  /// Bend-plane state at a reference z (the B0pf entrance), derived from the
  /// parabola in a uniform field: x(z) = x_ref + tx_ref (z-z_ref) - K (q/p) B_y
  /// (z-z_ref)^2 / 2.
  struct BendFit {
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
  /// Signed q/p = -2 c2 / (K B_y) and the entrance state from the parabola
  /// fit, with the coefficient covariance transformed linearly. `fieldY` is the
  /// ion-frame dipole component [T] along the hits; it must be non-zero. The
  /// parabola is used from `zFirst` (first hit) on; the gap from `zReference`
  /// to `zFirst` is stepped back analytically with `fieldYGap`, since the
  /// unrotated magnet puts a lower field there than along the stations.
  BendFit bendFitFromParabola(const StubFit& fit, double fieldY, double zReference, double zFirst,
                              double fieldYGap);
  /// Uniform-field convenience: the parabola itself is evaluated at zReference.
  BendFit bendFitFromParabola(const StubFit& fit, double fieldY, double zReference);

  struct EndpointCompatibility {
    double slopeY{};
    double beamResidual{};
    double score{};
    bool valid{false};
  };
  /// Bound |d^2 y/dz^2| [1/mm] from ion-frame field components [T], a bound
  /// on both transverse slopes, and the minimum momentum [GeV]. A missing
  /// momentum/field bound returns infinity, leaving rejection to the fit.
  double nonBendCurvatureBound(double fieldX, double fieldY, double fieldZ, double maxAbsSlope,
                               double pMin);

  /// Fast non-bending-plane compatibility for an outer-station hit pair.
  /// Allow for |y''| <= maxAbsCurvatureY downstream of zFieldEntrance [mm].
  /// The unchanged beamline tolerance is reapplied after the signed-curvature
  /// correction; this envelope only prevents premature straight-road rejection.
  EndpointCompatibility endpointCompatibility(const Point3& first, const Point3& last,
                                              double maxAbsSlope, double maxBeamResidual,
                                              bool constrainToBeamline,
                                              double maxAbsCurvatureY = 0.0,
                                              double zFieldEntrance   = 0.0);

  /// Remove the quadrupole bending y'' = kappa (q/p) Bx (uniform `fieldX` [T]
  /// from `zReference` on) from the non-bend coordinates so that a straight
  /// line describes the field-free upstream trajectory.
  std::vector<Point3> removeNonBendCurvature(const std::vector<Point3>& pts, double fieldX,
                                             double qOverP, double zReference);

  /// Fit the corrected coordinates, retaining their shared dependence on q/p
  /// for joint bend/non-bend covariance propagation (field held fixed).
  StubFit fitWithNonBendCorrection(const std::vector<Point3>& pts, double fieldX, double qOverP,
                                   double zReference);

  /// Diagonal (phi, theta, q/p) variance additions for multiple scattering
  /// and the trajectory model, both scaling with |q/p|.
  std::array<double, 3> scatteringCovarianceAdditions(double qOverP, double theta,
                                                      double scatteringScale,
                                                      double qOverPRelativeUncertainty);

  /// Seed parameters (loc0, loc1, phi, theta, q/p) on the origin perigee from
  /// the bend and non-bend fits, in the lab frame.
  std::array<double, 5> seedParametersFromFit(const BendFit& bendFit, const StubFit& nonBendFit,
                                              double crossingAngle, bool constrainToBeamline);

  /// Propagate the fitted coefficient covariance to
  /// (loc0,loc1,phi,theta,q/p), returned as a row-major 5x5 matrix.
  std::array<double, 25> seedCovarianceFromFit(const BendFit& bendFit, const StubFit& nonBendFit,
                                               double crossingAngle, bool constrainToBeamline);

  /// Combine row-major fit and beam-spot covariances. Missing fit angular or
  /// momentum variances require fallbacks even if the beam spot contributes;
  /// position fallbacks are needed only when neither source constrains them.
  std::array<double, 25>
  seedCovarianceWithFallbacks(const std::array<double, 25>& fitCovariance,
                              const std::array<double, 25>& beamSpotCovariance,
                              const std::array<double, 5>& fallbackVariances);

  /// Seed covariance contribution of the interaction-vertex spread, as a
  /// row-major 5x5 matrix, for the beamline-constrained mode only.
  ///
  /// That mode reports the parameters of a track leaving the origin. A vertex
  /// displaced by v gives different parameters, above all a different polar
  /// angle: theta moves by about theta * v_z / L over the lever arm L of 6 m
  /// to the B0pf entrance, which is two orders of magnitude larger than the
  /// hit-resolution term for a bunch tens of mm long. Returns
  /// J diag(sigma^2) J^T with J the derivative of the seed state with respect
  /// to the lab-frame vertex, so the loc0-phi and loc1-theta correlations the
  /// spread induces are carried too. Zero in the unconstrained mode, where the
  /// direction is measured from the stub and the vertex never enters.
  std::array<double, 25> beamSpotCovarianceAdditions(const BendFit& bendFit,
                                                     const StubFit& nonBendFit,
                                                     double crossingAngle, bool constrainToBeamline,
                                                     double sigmaX, double sigmaY, double sigmaZ);

  /// Perigee parameters of the straight ray leaving `ref` along `dir`, expressed
  /// on a perigee surface centred at `perigee`. All positions in mm.
  struct PerigeeParams {
    double loc0{};
    double loc1{};
    double phi{};
    double theta{};
  };
  PerigeeParams perigeeFromRay(const Point3& ref, const Point3& dir, const Point3& perigee);

  /// One physical B0 station after z-gap clustering of surface or hit z.
  struct StationInterval {
    double zMin{};
    double zMax{};
    double zMean{};
  };

  /// Single-linkage clustering of ion-frame z values. Values farther apart
  /// than `gap` start a new station. Official disks (~270 mm) split;
  /// realistic front/back faces (~7 mm) stay one station.
  std::vector<StationInterval> clusterStations(const std::vector<double>& zValues, double gap);

  /// Index of the cached station containing `z`, or -1 if none is within `gap`
  /// of that station's [zMin, zMax] interval.
  int assignStation(double z, const std::vector<StationInterval>& stations, double gap);

  /// Group hit indices by station. When `stations` is non-empty, assign each
  /// hit to the cached geometry cluster; otherwise cluster the hit z values.
  std::map<unsigned int, std::vector<std::size_t>>
  groupHitsByStation(const std::vector<double>& hitZ, const std::vector<StationInterval>& stations,
                     double gap);

  /// Charge hypotheses to emit for one stub candidate.
  ///
  /// The sign of the charge is only known if the fitted curvature is
  /// significantly non-zero. Exactly zero curvature and a non-positive or
  /// non-finite variance are maximally ambiguous and must yield both
  /// hypotheses -- inheriting the sign convention of a ternary on zero would
  /// silently drop every negative candidate at the ambiguous limit.
  ///
  /// `configuredCharge` of -1 or +1 pins the charge; any other value leaves it
  /// to the curvature. `testBothCharges` forces both regardless.
  std::vector<int> chargeHypotheses(double qOverP, double qOverPVariance,
                                    double minCurvatureSignificance, int inferredCharge,
                                    bool testBothCharges, int configuredCharge);

} // namespace b0stub

using B0TrackerStubSeederAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Dedicated seeder for the B0 tracker (dipole spectrometer at z ~ 6 m).
///
/// Fits one hit per station in the ion-rotated frame: a parabola in the bend
/// plane gives signed q/p and direction from the local dipole field, the
/// quadrupole Bx bending is removed from the non-bend plane, then the state is
/// back-extrapolated analytically to the origin and emitted on the origin
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
  double m_z_field_entrance{};
  /// Lab z of the B0pf entrance face when derived from the magnet geometry
  /// (0 when zFieldEntrance is configured explicitly). The magnet is not
  /// rotated with the beam, so its face is a plane of constant lab z whose
  /// ion-frame z depends on the track's x.
  double m_z_face_lab{0.0};
  std::vector<b0stub::StationInterval> m_stations;
  std::unordered_map<std::uint64_t, unsigned int> m_volume_to_station;
  /// False when the geometry has no B0 to seed; process() then emits nothing.
  bool m_enabled{true};
};

} // namespace eicrecon
