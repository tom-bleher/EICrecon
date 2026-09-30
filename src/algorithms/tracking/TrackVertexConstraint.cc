// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "TrackVertexConstraint.h"

#include <DD4hep/Detector.h>
#include <Evaluator/DD4hepUnits.h>
#include <algorithms/geo.h>
#include <edm4eic/Cov4f.h>
#include <edm4eic/Cov6f.h>
#include <edm4eic/TrackSeed.h>
#include <edm4eic/unit_system.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/Vector4f.h>
#include <podio/ObjectID.h>
#include <podio/RelationRange.h>
#include <Eigen/Cholesky>
#include <Eigen/LU> // IWYU pragma: keep
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <gsl/pointers>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "TrackVertexConstraintConfig.h"

namespace eicrecon {

namespace {

  // Order of the edm4eic track parameters and their covariance (the ACTS order)
  constexpr std::size_t kD0 = 0, kZ0 = 1, kPhi = 2, kTheta = 3, kQOverP = 4, kTime = 5;
  // A track whose point of closest approach is farther than this from the beam spot is not at
  // the perigee of the origin (e.g. kept on the perigee of its seed)
  constexpr double kMaxReferenceDistance = 1000. * edm4eic::unit::mm;

  using Vector6 = Eigen::Matrix<double, 6, 1>;
  using Matrix6 = Eigen::Matrix<double, 6, 6>;

  /// Luminous region spread (x, y, z) [mm] of the IP6 afterburner configurations (afterburner
  /// 0.2.1, 5000 events each); the centres are at the origin within their statistical error
  struct BeamConfiguration {
    double electronEnergy;
    double ionEnergy;
    std::array<double, 3> hiDiv;
    std::array<double, 3> hiAcc;
  };
  constexpr std::array<BeamConfiguration, 5> kBeamConfigurations{{
      {.electronEnergy = 5.,
       .ionEnergy      = 41.,
       .hiDiv          = {0.166, 0.0187, 37.5},
       .hiAcc          = {0.166, 0.0187, 37.5}},
      {.electronEnergy = 5.,
       .ionEnergy      = 100.,
       .hiDiv          = {0.125, 0.0079, 35.0},
       .hiAcc          = {0.135, 0.0090, 35.0}},
      {.electronEnergy = 10.,
       .ionEnergy      = 100.,
       .hiDiv          = {0.132, 0.0087, 35.0},
       .hiAcc          = {0.148, 0.0106, 35.0}},
      {.electronEnergy = 10.,
       .ionEnergy      = 275.,
       .hiDiv          = {0.111, 0.0059, 30.1},
       .hiAcc          = {0.150, 0.0109, 30.1}},
      {.electronEnergy = 18.,
       .ionEnergy      = 275.,
       .hiDiv          = {0.141, 0.0074, 30.2},
       .hiAcc          = {0.223, 0.0171, 30.2}},
  }};

  /// (d0, z0) at the perigee of the origin of a straight line with direction (phi, theta)
  /// through a point (x, y, z); linear in the point
  Eigen::Matrix<double, 2, 3> projection(double phi, double theta) {
    Eigen::Matrix<double, 2, 3> h;
    h << -std::sin(phi), std::cos(phi), 0., -std::cos(phi) / std::tan(theta),
        -std::sin(phi) / std::tan(theta), 1.;
    return h;
  }

  Vector6 state(const edm4eic::TrackParameters& p) {
    Vector6 x;
    x << p.getLoc().a, p.getLoc().b, p.getPhi(), p.getTheta(), p.getQOverP(), p.getTime();
    return x;
  }

  Matrix6 covariance(const edm4eic::TrackParameters& p) {
    Matrix6 c;
    for (unsigned i = 0; i < 6; ++i) {
      for (unsigned j = 0; j < 6; ++j) {
        c(i, j) = p.getCovariance()(i, j);
      }
    }
    return c;
  }

  /// The first parameters of the track's trajectory, which are those at its reference surface
  std::optional<edm4eic::TrackParameters> parameters(const edm4eic::Track& track) {
    const auto trajectory = track.getTrajectory();
    if (!trajectory.isAvailable() || trajectory.getTrackParameters().empty()) {
      return std::nullopt;
    }
    return trajectory.getTrackParameters()[0];
  }

} // namespace

std::optional<TrackVertexConstraint::BeamSpot>
TrackVertexConstraint::beamSpot(double electronEnergy, double ionEnergy,
                                const std::string& setting) {
  if (setting != "hiDiv" && setting != "hiAcc") {
    return std::nullopt;
  }
  for (const auto& c : kBeamConfigurations) {
    if (std::abs(c.electronEnergy - electronEnergy) < 0.5 &&
        std::abs(c.ionEnergy - ionEnergy) < 0.5) {
      const auto& s = setting == "hiDiv" ? c.hiDiv : c.hiAcc;
      return BeamSpot{.sigma = {s[0], s[1], s[2]}};
    }
  }
  return std::nullopt;
}

void TrackVertexConstraint::init() {
  m_spot      = {m_cfg.beamSpotX, m_cfg.beamSpotY, m_cfg.beamSpotZ};
  m_spotSigma = {m_cfg.beamSpotSigmaX, m_cfg.beamSpotSigmaY, m_cfg.beamSpotSigmaZ};
  if ((m_spotSigma.array() > 0.).all()) {
    return;
  }
  const auto detector   = algorithms::GeoSvc::instance().detector();
  double electronEnergy = 0.;
  double ionEnergy      = 0.;
  try {
    electronEnergy = detector->constant<double>("ElectronBeamEnergy") / dd4hep::GeV;
    ionEnergy      = detector->constant<double>("IonBeamEnergy") / dd4hep::GeV;
  } catch (const std::exception& e) {
    error("The geometry has no beam energies ({}): set the beam spot spreads", e.what());
    throw;
  }
  const auto spot = beamSpot(electronEnergy, ionEnergy, m_cfg.beamSetting);
  if (!spot) {
    error("No luminous region for {:.0f}x{:.0f} GeV with beam setting '{}': set the beam spot "
          "spreads",
          electronEnergy, ionEnergy, m_cfg.beamSetting);
    throw std::runtime_error("unknown beam configuration for the vertex constraint");
  }
  for (int i = 0; i < 3; ++i) {
    if (!(m_spotSigma(i) > 0.)) {
      m_spotSigma(i) = spot->sigma(i);
    }
  }
  info("Luminous region of {:.0f}x{:.0f} GeV, {}: sigma ({}, {}, {}) mm", electronEnergy, ionEnergy,
       m_cfg.beamSetting, m_spotSigma.x(), m_spotSigma.y(), m_spotSigma.z());
}

TrackVertexConstraint::Vertex
TrackVertexConstraint::estimateVertex(const edm4eic::VertexCollection& vertices,
                                      const edm4eic::TrackCollection& vertexTracks) const {
  Vertex v{.position = m_spot, .covariance = m_spotSigma.cwiseAbs2().asDiagonal()};

  if (!vertices.empty()) {
    const auto p = vertices[0].getPosition();
    const auto e = vertices[0].getPositionError();
    Eigen::Matrix3d r;
    r << e.xx, e.xy, e.xz, e.xy, e.yy, e.yz, e.xz, e.yz, e.zz;
    const Eigen::Matrix3d k = v.covariance * (v.covariance + r).inverse();
    v.position += k * (Eigen::Vector3d{p.x, p.y, p.z} - v.position);
    v.covariance = (Eigen::Matrix3d::Identity() - k) * v.covariance;
    return v;
  }

  // The vertex tracks at the perigee of the origin, best measured z0 first
  std::vector<edm4eic::TrackParameters> tracks;
  for (const auto& track : vertexTracks) {
    const auto p   = parameters(track);
    const auto ref = track.getPosition();
    if (p && (Eigen::Vector3d{ref.x, ref.y, ref.z} - m_spot).norm() < kMaxReferenceDistance) {
      tracks.push_back(*p);
    }
  }
  std::ranges::stable_sort(tracks, {}, [](const auto& p) { return p.getCovariance()(kZ0, kZ0); });

  for (const auto& p : tracks) {
    const Matrix6 c                = covariance(p);
    const auto h                   = projection(p.getPhi(), p.getTheta());
    const Eigen::Vector2d residual = Eigen::Vector2d{p.getLoc().a, p.getLoc().b} - h * v.position;
    const Eigen::Matrix2d s        = h * v.covariance * h.transpose() + c.topLeftCorner<2, 2>();
    const double chi2              = residual.dot(s.ldlt().solve(residual));
    if (!(chi2 < m_cfg.maxChi2)) {
      trace("vertex track with d0 {} z0 {} incompatible with the vertex, chi2 {}", p.getLoc().a,
            p.getLoc().b, chi2);
      continue;
    }
    const Eigen::Matrix<double, 3, 2> k = v.covariance * h.transpose() * s.inverse();
    v.position += k * residual;
    v.covariance = (Eigen::Matrix3d::Identity() - k * h) * v.covariance;
  }
  return v;
}

std::optional<std::pair<edm4eic::MutableTrackParameters, double>>
TrackVertexConstraint::constrain(const edm4eic::Track& track, const Vertex& vertex) const {
  const auto p = parameters(track);
  if (!p) {
    return std::nullopt;
  }
  const auto ref = track.getPosition();
  if (!((Eigen::Vector3d{ref.x, ref.y, ref.z} - m_spot).norm() < kMaxReferenceDistance)) {
    return std::nullopt;
  }

  Vector6 x                      = state(*p);
  Matrix6 c                      = covariance(*p);
  const auto h                   = projection(x(kPhi), x(kTheta));
  const Eigen::Vector2d residual = h * vertex.position - x.head<2>();
  const Eigen::Matrix2d s        = c.topLeftCorner<2, 2>() + h * vertex.covariance * h.transpose();
  const double chi2              = residual.dot(s.ldlt().solve(residual));
  if (!(chi2 < m_cfg.maxChi2)) {
    return std::nullopt;
  }
  const Eigen::Matrix<double, 6, 2> k = c.leftCols<2>() * s.inverse();
  x += k * residual;
  c -= k * c.topRows<2>();

  auto out = p->clone();
  out.setLoc({static_cast<float>(x(kD0)), static_cast<float>(x(kZ0))});
  out.setPhi(static_cast<float>(x(kPhi)));
  out.setTheta(static_cast<float>(x(kTheta)));
  out.setQOverP(static_cast<float>(x(kQOverP)));
  out.setTime(static_cast<float>(x(kTime)));
  edm4eic::Cov6f cov;
  for (unsigned i = 0; i < 6; ++i) {
    for (unsigned j = i; j < 6; ++j) {
      cov(i, j) = static_cast<float>(0.5 * (c(i, j) + c(j, i)));
    }
  }
  out.setCovariance(cov);
  return std::pair{out, chi2};
}

void TrackVertexConstraint::process(const Input& input, const Output& output) const {
  const auto [tracks, associations, vertices, vertexTracks]         = input;
  auto [outParameters, outTrajectories, outTracks, outAssociations] = output;
  if (tracks->empty()) {
    return;
  }
  const auto vertex = estimateVertex(*vertices, *vertexTracks);
  debug("vertex ({}, {}, {}) mm, sigma ({}, {}, {}) mm", vertex.position.x(), vertex.position.y(),
        vertex.position.z(), std::sqrt(vertex.covariance(0, 0)), std::sqrt(vertex.covariance(1, 1)),
        std::sqrt(vertex.covariance(2, 2)));

  for (const auto& track : *tracks) {
    // The same measurements, and the input track as its segment
    auto out = track.clone();
    out.addToTracks(track);

    auto constrained = constrain(track, vertex);
    if (constrained) {
      const auto& [pars, chi2] = *constrained;
      const double d0          = pars.getLoc().a;
      const double phi         = pars.getPhi();
      const double theta       = pars.getTheta();
      const double p           = std::abs(1. / pars.getQOverP());
      out.setPosition({static_cast<float>(-d0 * std::sin(phi)),
                       static_cast<float>(d0 * std::cos(phi)), pars.getLoc().b});
      out.setMomentum({static_cast<float>(p * std::sin(theta) * std::cos(phi)),
                       static_cast<float>(p * std::sin(theta) * std::sin(phi)),
                       static_cast<float>(p * std::cos(theta))});
      out.setTime(pars.getTime());
      out.setChi2(track.getChi2() + static_cast<float>(chi2));
      out.setNdf(track.getNdf() + 2);
      outParameters->push_back(pars);
    } else if (const auto p = parameters(track)) {
      outParameters->push_back(p->clone());
    } else {
      outParameters->create();
    }

    // A trajectory carrying the new parameters
    auto trajectory = outTrajectories->create();
    if (const auto in = track.getTrajectory(); in.isAvailable()) {
      trajectory.setType(in.getType());
      trajectory.setNStates(in.getNStates());
      trajectory.setNMeasurements(in.getNMeasurements());
      trajectory.setNOutliers(in.getNOutliers());
      trajectory.setNHoles(in.getNHoles());
      trajectory.setNSharedHits(in.getNSharedHits());
      for (const auto chi2 : in.getMeasurementChi2()) {
        trajectory.addToMeasurementChi2(chi2);
      }
      for (const auto chi2 : in.getOutlierChi2()) {
        trajectory.addToOutlierChi2(chi2);
      }
      for (const auto& m : in.getMeasurements_deprecated()) {
        trajectory.addToMeasurements_deprecated(m);
      }
      for (const auto& m : in.getOutliers_deprecated()) {
        trajectory.addToOutliers_deprecated(m);
      }
      trajectory.setSeed(in.getSeed());
    }
    if (constrained || parameters(track)) {
      trajectory.addToTrackParameters((*outParameters)[outParameters->size() - 1]);
    }
    out.setTrajectory(trajectory);
    outTracks->push_back(out);
  }

  // The associations of the input tracks, pointing to the output tracks
  for (const auto& a : *associations) {
    const auto rec = a.getRec();
    if (rec.getObjectID().collectionID != tracks->getID()) {
      continue;
    }
    auto copy = outAssociations->create();
    copy.setWeight(a.getWeight());
    copy.setRec((*outTracks)[rec.getObjectID().index]);
    copy.setSim(a.getSim());
  }
}

} // namespace eicrecon
