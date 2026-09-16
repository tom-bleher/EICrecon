// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "ActsGeometryProvider.h"
#include "B0StationMap.h"
#include "B0TelescopeMath.h"
#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <DD4hep/Detector.h>
#include <Evaluator/DD4hepUnits.h>
#include <edm4eic/Measurement2DCollection.h>
#include <Eigen/Geometry>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>

namespace eicrecon::b0 {
/// A single geometry-ID/physical-station definition used by the local seeder,
/// truth selector and final station-quality check. No hit-based renumbering.
class TelescopeGeometry {
public:
  TelescopeGeometry(const ActsGeometryProvider& provider, double gap)
      : m_provider(provider), m_stations({}) {
    const auto* detector = provider.dd4hepDetector();
    if (!detector) throw std::runtime_error("B0: missing DD4hep geometry");
    m_angle = detector->constant<double>("CrossingAngle") / dd4hep::rad;
    if (!std::isfinite(m_angle)) throw std::runtime_error("B0 crossing angle is not finite");
    m_rotation = Eigen::AngleAxisd(-m_angle, Eigen::Vector3d::UnitY()).toRotationMatrix();
    const auto detectorId = detector->constant<unsigned long>("B0Tracker_Station_1_ID") & 0xff;
    std::vector<StationMap::Surface> centres;
    for (const auto& [volume, surface] : provider.surfaceMap()) {
      if (!surface || surface->geometryId().sensitive() == 0 ||
          surface->geometryId().extra() != detectorId) continue;
      if (surface->type() != Acts::Surface::Plane)
        throw std::runtime_error("B0 local seeding requires planar sensitive surfaces");
      const auto id = surface->geometryId().value();
      m_surfaces.emplace(id, surface);
      m_volumeToSurface.emplace(volume, id);
      centres.push_back({id, (m_rotation * surface->center(provider.getActsGeometryContext())).z() /
                               Acts::UnitConstants::mm});
    }
    m_stations = StationMap(std::move(centres), gap);
    if (m_stations.stations().empty()) throw std::runtime_error("B0: no sensitive surfaces found");
  }
  const StationMap& stations() const { return m_stations; }
  const Eigen::Matrix3d& labToIon() const { return m_rotation; }
  const Acts::Surface* surface(std::uint64_t id) const {
    const auto it = m_surfaces.find(id);
    return it == m_surfaces.end() ? nullptr : it->second;
  }
  const Acts::Surface* volumeSurface(std::uint64_t volume) const {
    const auto it = m_volumeToSurface.find(volume);
    return it == m_volumeToSurface.end() ? nullptr : surface(it->second);
  }
  std::vector<TelescopeHit> measurements(const edm4eic::Measurement2DCollection& input) const {
    std::vector<TelescopeHit> result;
    const auto& context = m_provider.getActsGeometryContext();
    for (std::size_t i = 0; i < input.size(); ++i) {
      const auto m = input[i];
      const auto* s = surface(m.getSurface());
      if (!s) throw std::runtime_error("B0 measurement has an unmapped sensitive surface");
      const auto local = m.getLoc();
      const auto global = s->localToGlobal(context,
          Acts::Vector2(local.a, local.b) * Acts::UnitConstants::mm, Acts::Vector3::UnitZ());
      Eigen::Matrix2d covariance;
      const auto cov = m.getCovariance();
      covariance << cov.xx, cov.xy, cov.xy, cov.yy;
      const Eigen::Matrix<double, 3, 2> axes = m_rotation * s->transform(context).linear().leftCols<2>();
      TelescopeHit h;
      h.index = i; h.surface = m.getSurface();
      h.station = m_stations.stationForSurface(h.surface).value();
      h.position = m_rotation * global / Acts::UnitConstants::mm;
      h.covariance = axes * covariance * axes.transpose();
      h.time = m.getTime(); h.timeVariance = cov.zz;
      result.push_back(h);
    }
    return result;
  }
  /// IDs are strings: JSON double precision cannot represent a uint64 ID.
  void writeManifest(const std::string& path) const {
    if (path.empty()) return;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write B0 station manifest: " + path);
    out << std::setprecision(17) << "{\"schema_version\":1,\"crossing_angle_rad\":" << m_angle
        << ",\"station_gap_mm\":" << m_stations.gap() << ",\"surfaces\":[";
    bool first = true;
    for (const auto& [id, s] : m_surfaces) {
      const auto& context = m_provider.getActsGeometryContext();
      const Acts::Vector3 p = s->center(context) / Acts::UnitConstants::mm;
      const auto n = s->transform(context).linear().col(2).eval();
      if (!first) out << ',';
      first = false;
      out << "{\"id\":\"" << id << "\",\"station\":" << m_stations.stationForSurface(id).value()
          << ",\"layer\":" << s->geometryId().layer() << ",\"z_ion_mm\":" << (m_rotation * p).z()
          << ",\"position_mm\":[" << p.x() << ',' << p.y() << ',' << p.z()
          << "],\"normal\":[" << n.x() << ',' << n.y() << ',' << n.z() << "]}";
    }
    out << "]}\n";
    if (!out) throw std::runtime_error("Failed while writing B0 station manifest: " + path);
  }
private:
  const ActsGeometryProvider& m_provider;
  double m_angle{};
  Eigen::Matrix3d m_rotation = Eigen::Matrix3d::Identity();
  StationMap m_stations;
  std::map<std::uint64_t, const Acts::Surface*> m_surfaces;
  std::map<std::uint64_t, std::uint64_t> m_volumeToSurface;
};
} // namespace eicrecon::b0
