// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include "B0ACLGAD_factory.h"
#include "algorithms/interfaces/ActsSvc.h"
#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Readout.h>
#include <DD4hep/detail/SegmentationsInterna.h>
#include <DDSegmentation/CartesianGridXY.h>
#include <TGeoBBox.h>
#include <edm4eic/Cov3f.h>
#include <edm4eic/CovDiag3f.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace eicrecon {
namespace {
std::uint64_t identity(const edm4hep::SimTrackerHit& hit) {
  const auto id = hit.getObjectID();
  if (!hit.isAvailable() || id.index < 0) throw std::runtime_error("Unavailable B0 sim hit");
  return (std::uint64_t(static_cast<std::uint32_t>(id.collectionID)) << 32) |
         static_cast<std::uint32_t>(id.index);
}
const dd4hep::DDSegmentation::CartesianGridXY& grid(const dd4hep::Detector& detector) {
  const auto segmentation = detector.readout("B0TrackerHits").segmentation();
  const auto* result = dynamic_cast<const dd4hep::DDSegmentation::CartesianGridXY*>(segmentation->segmentation);
  if (!result) throw std::runtime_error("B0 AC-LGAD requires CartesianGridXY");
  return *result;
}
bool close(double a, double b) { return std::abs(a - b) < 1e-8; }
}
void B0ACLGAD_factory::Configure() {
  if (config().calibrationFile.empty()) throw std::invalid_argument("B0 AC-LGAD requires calibrationFile");
  std::ifstream calibration(config().calibrationFile);
  if (!calibration) throw std::runtime_error("Cannot open B0 AC-LGAD calibration: " + config().calibrationFile);
  m_calibration = b0::aclgad::readCalibration(calibration, config().allowExperimental);
  m_provider = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (!m_provider || !m_provider->dd4hepDetector()) throw std::runtime_error("B0 AC-LGAD requires ACTS/DD4hep geometry");
  const auto& detector = *m_provider->dd4hepDetector();
  // No compatibility fallback: a 70 um file must not be reinterpreted as 500 um channels.
  if (detector.constant<int>("B0ACLGADReadoutVersion") != 1)
    throw std::runtime_error("Wrong B0 readout profile; regenerate simulation with the physical-pitch profile");
  const auto& segmentation = grid(detector);
  const auto& c = m_calibration.config;
  if (c.timeStep < .001 || !close(c.timeStep * 1000, std::round(c.timeStep * 1000)))
    throw std::runtime_error("B0 raw timestamps require a TDC step that is an integer number of picoseconds");
  if (!close(segmentation.gridSizeX() / dd4hep::mm, c.pitchX) ||
      !close(segmentation.gridSizeY() / dd4hep::mm, c.pitchY) ||
      !close(segmentation.offsetX() / dd4hep::mm, c.offsetX) ||
      !close(segmentation.offsetY() / dd4hep::mm, c.offsetY) ||
      segmentation.fieldNameX() != "x" || segmentation.fieldNameY() != "y")
    throw std::runtime_error("B0 response pitch/offset/field names do not match geometry");
  const auto system = detector.constant<unsigned long>("B0Tracker_Station_1_ID") & 0xff;
  const auto manager = detector.volumeManager();
  const auto& geometry = m_provider->getActsGeometryContext();
  for (const auto& [volume, surface] : m_provider->surfaceMap()) {
    if (!surface || surface->geometryId().sensitive() == 0 || surface->geometryId().extra() != system) continue;
    if (surface->type() != Acts::Surface::Plane) throw std::runtime_error("B0 response requires sensitive planes");
    const auto* context = manager.lookupContext(volume);
    if (!context) throw std::runtime_error("B0 sensitive surface has no DD4hep context");
    const auto* box = dynamic_cast<const TGeoBBox*>(context->volumePlacement().volume().solid().ptr());
    if (!box) throw std::runtime_error("B0 response requires rectangular sensors");
    Binding binding{context, surface, box->GetDX() / dd4hep::mm, box->GetDY() / dd4hep::mm,
                    box->GetDZ() / dd4hep::mm, Eigen::Matrix2d::Zero(), Acts::Vector2::Zero()};
    b0::aclgad::channelRange(binding.halfX, c.pitchX, c.offsetX);
    b0::aclgad::channelRange(binding.halfY, c.pitchY, c.offsetY);
    auto local = [&](double x, double y) -> Acts::Vector2 {
      const auto world = context->localToWorld(dd4hep::Position{x * dd4hep::mm, y * dd4hep::mm, 0});
      const Acts::Vector3 global = Acts::Vector3(world.x(), world.y(), world.z()) * (Acts::UnitConstants::mm / dd4hep::mm);
      const auto result = surface->globalToLocal(geometry, global, surface->transform(geometry).linear().col(2));
      if (!result.ok()) throw std::runtime_error("B0 sensor/ACTS plane correspondence failed");
      return result.value() / Acts::UnitConstants::mm;
    };
    binding.origin = local(0, 0);
    binding.jacobian.col(0) = local(1, 0) - binding.origin;
    binding.jacobian.col(1) = local(0, 1) - binding.origin;
    if (!binding.jacobian.allFinite() ||
        !(binding.jacobian.transpose() * binding.jacobian).isApprox(Eigen::Matrix2d::Identity(), 1e-7))
      throw std::runtime_error("B0 readout axes are not an orthonormal ACTS plane mapping");
    if (!m_bindings.emplace(context->identifier, binding).second)
      throw std::runtime_error("Duplicate B0 response sensor identity");
    m_sensors.push_back({context->identifier, binding.halfX, binding.halfY});
  }
  if (m_sensors.empty()) throw std::runtime_error("B0 response found no sensor surfaces");
  logger()->info("B0 AC-LGAD status={} reference={} sensors={} pitch=({}, {}) mm",
      m_calibration.status, m_calibration.reference, m_sensors.size(), c.pitchX, c.pitchY);
}
void B0ACLGAD_factory::Process(int32_t, uint64_t) {
  const auto* headers = m_headers();
  const auto* simHits = m_simHits();
  if (headers->size() != 1) throw std::runtime_error("B0 response requires exactly one EventHeader");
  const auto manager = m_provider->dd4hepDetector()->volumeManager();
  const auto& segmentation = grid(*m_provider->dd4hepDetector());
  const auto& c = m_calibration.config;
  std::vector<b0::aclgad::Deposit> deposits;
  std::map<std::uint64_t, edm4hep::SimTrackerHit> truth;
  for (const auto& hit : *simHits) {
    const auto* context = manager.lookupContext(hit.getCellID());
    if (!context || !m_bindings.count(context->identifier)) throw std::runtime_error("Unmapped B0 simulation hit");
    const auto& binding = m_bindings.at(context->identifier);
    const auto p = hit.getPosition();
    const auto local = context->worldToLocal(dd4hep::Position{p.x * dd4hep::mm, p.y * dd4hep::mm, p.z * dd4hep::mm});
    if (!std::isfinite(local.z()) || std::abs(local.z() / dd4hep::mm) > binding.halfZ + 1e-3)
      throw std::runtime_error("B0 simulation hit is outside its sensitive thickness");
    const auto key = identity(hit);
    if (!truth.emplace(key, hit).second) throw std::runtime_error("Duplicate B0 simulation hit");
    deposits.push_back({key, context->identifier, local.x() / dd4hep::mm, local.y() / dd4hep::mm,
                        hit.getEDep() * 1e9, hit.getTime()}); // EDM energy GeV -> eV
  }
  const auto response = b0::aclgad::simulate(m_sensors, std::move(deposits), c,
                                            m_uid.getUniqueID(*headers, "B0ACLGADResponse"));
  auto* raw = m_raw().get(); auto* channelHits = m_channelHits().get();
  for (const auto& channel : response.channels) {
    const auto& binding = m_bindings.at(channel.sensor);
    const dd4hep::Position local{channel.x * dd4hep::mm, channel.y * dd4hep::mm, 0};
    const auto global = binding.context->localToWorld(local);
    const dd4hep::DDSegmentation::Vector3D lp(local.x(), local.y(), local.z());
    const dd4hep::DDSegmentation::Vector3D gp(global.x(), global.y(), global.z());
    const auto cell = segmentation.cellID(lp, gp, channel.sensor);
    if ((cell & binding.context->mask) != channel.sensor)
      throw std::runtime_error("B0 channel encoding corrupted sensor volume bits");
    const auto roundtrip = segmentation.position(cell);
    if (!close(roundtrip.x() / dd4hep::mm, channel.x) || !close(roundtrip.y() / dd4hep::mm, channel.y))
      throw std::runtime_error("B0 physical channel ID round trip failed");
    const double timestamp = std::round(channel.time * 1000);
    if (!std::isfinite(timestamp) || timestamp < std::numeric_limits<std::int32_t>::min() ||
        timestamp > std::numeric_limits<std::int32_t>::max()) throw std::overflow_error("B0 TDC timestamp overflow");
    auto rawHit = raw->create(); rawHit.setCellID(cell); rawHit.setCharge(channel.adc);
    rawHit.setTimeStamp(static_cast<std::int32_t>(timestamp));
    auto rec = channelHits->create(); rec.setCellID(cell); rec.setRawHit(rawHit);
    rec.setPosition({static_cast<float>(global.x() / dd4hep::mm), static_cast<float>(global.y() / dd4hep::mm),
                     static_cast<float>(global.z() / dd4hep::mm)});
    // These are pad-center diagnostics. Only cluster Measurement2D covariance is used for tracking.
    rec.setPositionError({c.pitchX * c.pitchX / 12, c.pitchY * c.pitchY / 12, 0});
    rec.setTime(channel.time); rec.setTimeError(std::sqrt(c.sensorTimeSigma*c.sensorTimeSigma +
        c.clockTimeSigma*c.clockTimeSigma + channel.independentTimeVariance));
    rec.setEdep(channel.charge / c.gain * c.pairEnergy * 1e-9); rec.setEdepError(0);
    for (const auto& [key, weight] : channel.truth) {
      auto link = m_links()->create(); link.setFrom(rawHit); link.setTo(truth.at(key)); link.setWeight(weight);
      auto association = m_associations()->create(); association.setRawHit(rawHit);
      association.setSimHit(truth.at(key)); association.setWeight(weight);
    }
  }
  for (const auto& cluster : response.clusters) {
    const auto& binding = m_bindings.at(cluster.sensor);
    const auto local = binding.origin + binding.jacobian * Acts::Vector2(cluster.x, cluster.y);
    Eigen::Matrix2d covariance;
    covariance << cluster.xx, cluster.xy, cluster.xy, cluster.yy;
    covariance = (binding.jacobian * covariance * binding.jacobian.transpose()).eval();
    const auto cross = binding.jacobian * Acts::Vector2(cluster.xt, cluster.yt);
    auto measurement = m_measurements()->create(); measurement.setSurface(binding.surface->geometryId().value());
    measurement.setLoc({static_cast<float>(local.x()), static_cast<float>(local.y())});
    measurement.setTime(cluster.time);
    measurement.setCovariance({covariance(0,0), covariance(1,1), cluster.timeVariance,
                               covariance(0,1), cross.x(), cross.y()});
    for (auto index : cluster.channels) {
      measurement.addToHits((*channelHits)[index]);
      measurement.addToWeights(response.channels[index].charge / cluster.charge);
    }
  }
  logger()->debug("B0 AC-LGAD deposits={} channels={} clusters={} gateRejected={} belowThreshold={} saturated={} rejectedSaturatedClusters={} work={} edgeLoss_e={}",
      simHits->size(), response.channels.size(), response.clusters.size(), response.outsideGate,
      response.belowThreshold, response.saturatedChannels, response.rejectedSaturatedClusters,
      response.work, response.inputSignal - response.sharedSignal);
}
} // namespace eicrecon
