// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/Cov4f.h>
#include <edm4eic/Cov6f.h>
#include <edm4eic/MCRecoTrackParticleAssociationCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrajectoryCollection.h>
#include <edm4eic/VertexCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/Vector4f.h>
#include <podio/RelationRange.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <string>

#include "algorithms/tracking/TrackVertexConstraint.h"
#include "algorithms/tracking/TrackVertexConstraintConfig.h"

namespace {

// Like a B0 track at the perigee of the origin: d0 known to 2 mm and anticorrelated with phi, z0
// to 60 mm and correlated with theta
constexpr double sigmaD0    = 2.;
constexpr double sigmaZ0    = 60.;
constexpr double sigmaPhi   = 1.e-2;
constexpr double sigmaTheta = 3.e-4;
constexpr double rho        = 0.99;

/// The 5x41 GeV luminous region, given explicitly: the tests have no geometry to take it from
eicrecon::TrackVertexConstraintConfig config() {
  eicrecon::TrackVertexConstraintConfig cfg;
  cfg.beamSpotX = cfg.beamSpotY = cfg.beamSpotZ = 0.;
  cfg.beamSpotSigmaX                            = 0.17;
  cfg.beamSpotSigmaY                            = 0.019;
  cfg.beamSpotSigmaZ                            = 38.;
  return cfg;
}

struct Output {
  edm4eic::TrackParametersCollection params;
  edm4eic::TrajectoryCollection trajectories;
  edm4eic::TrackCollection tracks;
  edm4eic::MCRecoTrackParticleAssociationCollection associations;
};

struct Event {
  edm4eic::TrackParametersCollection params;
  edm4hep::MCParticleCollection particles;
  edm4eic::MCRecoTrackParticleAssociationCollection associations;
  edm4eic::TrajectoryCollection trajectories;
  edm4eic::TrackCollection tracks;
  edm4eic::TrackCollection vertexTracks;
  edm4eic::VertexCollection vertices;

  /// A track with parameters (d0, z0, phi, theta) and the given point of closest approach
  void addTrack(edm4eic::TrackCollection& to, double d0, double z0, double phi, double theta,
                const std::array<double, 4>& sigma, double correlation,
                edm4hep::Vector3f reference = {}) {
    auto p = params.create();
    p.setLoc({static_cast<float>(d0), static_cast<float>(z0)});
    p.setPhi(static_cast<float>(phi));
    p.setTheta(static_cast<float>(theta));
    p.setQOverP(0.05F);
    edm4eic::Cov6f cov;
    for (unsigned i = 0; const double si : sigma) {
      cov(i, i) = static_cast<float>(si * si);
      ++i;
    }
    cov(0, 2) = static_cast<float>(-correlation * sigma[0] * sigma[2]);
    cov(1, 3) = static_cast<float>(correlation * sigma[1] * sigma[3]);
    cov(4, 4) = 1.e-6F;
    cov(5, 5) = 1.F;
    p.setCovariance(cov);
    auto trajectory = trajectories.create();
    trajectory.addToTrackParameters(p);
    auto track = to.create();
    track.setTrajectory(trajectory);
    track.setPosition(reference);
    track.setChi2(4.F);
    track.setNdf(5);
    if (&to == &tracks) {
      auto a = associations.create();
      a.setRec(track);
      a.setSim(particles.create());
      a.setWeight(1.F);
    }
  }

  Output run(const eicrecon::TrackVertexConstraintConfig& cfg = config()) {
    eicrecon::TrackVertexConstraint algo("TrackVertexConstraint");
    algo.applyConfig(cfg);
    algo.init();
    Output out;
    algo.process({&tracks, &associations, &vertices, &vertexTracks},
                 {&out.params, &out.trajectories, &out.tracks, &out.associations});
    return out;
  }
};

constexpr std::array<double, 4> b0Sigma{sigmaD0, sigmaZ0, sigmaPhi, sigmaTheta};

} // namespace

TEST_CASE("A track is pulled to the vertex, and its correlated angles with it",
          "[TrackVertexConstraint]") {
  Event event;
  const double phi   = std::numbers::pi;
  const double theta = 0.03;
  event.addTrack(event.tracks, 1., 20., phi, theta, b0Sigma, rho);
  // A point-like beam spot at the origin
  auto cfg           = config();
  cfg.beamSpotSigmaX = cfg.beamSpotSigmaY = cfg.beamSpotSigmaZ = 1.e-3;
  const auto result                                            = event.run(cfg);
  const auto& out                                              = result.params;
  REQUIRE(out.size() == 1);
  // The track now passes through the origin: d0 and z0 near 0, with small errors
  CHECK(out[0].getLoc().a == Catch::Approx(0.).margin(1.e-2));
  CHECK(out[0].getLoc().b == Catch::Approx(0.).margin(1.e-1));
  CHECK(out[0].getCovariance()(0, 0) < 1.e-3F);
  // phi moves by cov(phi, d0) / var(d0) times the d0 residual: -rho * sigmaPhi / sigmaD0 * (0 - 1)
  CHECK(out[0].getPhi() == Catch::Approx(phi + rho * sigmaPhi / sigmaD0 * 1.).epsilon(1.e-4));
  CHECK(out[0].getTheta() ==
        Catch::Approx(theta - rho * sigmaTheta / sigmaZ0 * 20.).epsilon(1.e-4));
  CHECK(out[0].getCovariance()(2, 2) < sigmaPhi * sigmaPhi * (1. - rho * rho) * 1.01);

  // The track carries the constrained momentum, the vertex in its chi2, and the input track
  REQUIRE(result.tracks.size() == 1);
  const auto track = result.tracks[0];
  const double p   = 1. / 0.05;
  CHECK(track.getMomentum().x ==
        Catch::Approx(p * std::sin(out[0].getTheta()) * std::cos(out[0].getPhi())).epsilon(1.e-5));
  CHECK(track.getMomentum().z == Catch::Approx(p * std::cos(out[0].getTheta())).epsilon(1.e-5));
  CHECK(track.getNdf() == 7);
  CHECK(track.getChi2() > 4.F);
  // whose trajectory carries the constrained parameters, as read downstream
  REQUIRE(track.getTrajectory().getTrackParameters().size() == 1);
  CHECK(track.getTrajectory().getTrackParameters()[0] == out[0]);
  REQUIRE(track.getTracks().size() == 1);
  CHECK(track.getTracks()[0] == event.tracks[0]);
  // and its association is that of the input track
  REQUIRE(result.associations.size() == 1);
  CHECK(result.associations[0].getRec() == track);
  CHECK(result.associations[0].getSim() == event.associations[0].getSim());
}

TEST_CASE("Displaced tracks and tracks away from the origin are left alone",
          "[TrackVertexConstraint]") {
  Event event;
  // d0 = 50 mm, 25 sigma from the beam spot
  event.addTrack(event.tracks, 50., 20., std::numbers::pi, 0.03, b0Sigma, rho);
  // expressed on a perigee near the first B0 hit
  event.addTrack(event.tracks, 0., 0., std::numbers::pi, 0.03, b0Sigma, rho, {-150.F, 0.F, 5840.F});
  const auto result = event.run();
  REQUIRE(result.params.size() == 2);
  REQUIRE(result.tracks.size() == 2);
  for (std::size_t i = 0; i < 2; ++i) {
    CHECK(result.params[i].getLoc().a == event.params[i].getLoc().a);
    CHECK(result.params[i].getPhi() == event.params[i].getPhi());
    CHECK(result.tracks[i].getNdf() == 5);
  }
}

TEST_CASE("The vertex z comes from a central track, or from a reconstructed vertex",
          "[TrackVertexConstraint]") {
  eicrecon::TrackVertexConstraint algo("TrackVertexConstraint");
  algo.applyConfig(config());
  algo.init();
  Event event;
  // A central track at 90 degrees measures z0 = z of the vertex
  constexpr std::array<double, 4> centralSigma{0.02, 0.05, 1.e-3, 1.e-3};
  event.addTrack(event.vertexTracks, 0., 12., 0., std::numbers::pi / 2., centralSigma, 0.);
  // and an outlier far from the beam spot, rejected
  event.addTrack(event.vertexTracks, 5., -300., 0., std::numbers::pi / 2., centralSigma, 0.);
  auto vertex = algo.estimateVertex(event.vertices, event.vertexTracks);
  CHECK(vertex.position.z() == Catch::Approx(12.).margin(1.e-3));
  CHECK(std::sqrt(vertex.covariance(2, 2)) == Catch::Approx(0.05).epsilon(1.e-2));

  // A reconstructed vertex takes precedence over the tracks
  auto v = event.vertices.create();
  v.setPosition({0.F, 0.F, -7.F, 0.F});
  v.setPositionError(edm4eic::Cov4f{0.01F, 0.01F, 0.01F, 1.F});
  vertex = algo.estimateVertex(event.vertices, event.vertexTracks);
  CHECK(vertex.position.z() == Catch::Approx(-7.).margin(1.e-3));
}

TEST_CASE("The luminous region follows the beam configuration", "[TrackVertexConstraint]") {
  using eicrecon::TrackVertexConstraint;
  const auto hiDiv = TrackVertexConstraint::beamSpot(18., 275., "hiDiv");
  const auto hiAcc = TrackVertexConstraint::beamSpot(18., 275., "hiAcc");
  REQUIRE(hiDiv);
  REQUIRE(hiAcc);
  CHECK(hiDiv->sigma.x() < hiAcc->sigma.x());
  CHECK(TrackVertexConstraint::beamSpot(5., 41., "hiDiv")->sigma.z() == Catch::Approx(37.5));
  CHECK_FALSE(TrackVertexConstraint::beamSpot(10., 130., "hiDiv"));
  CHECK_FALSE(TrackVertexConstraint::beamSpot(5., 41., "other"));
}
