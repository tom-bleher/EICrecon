// SPDX-License-Identifier: LGPL-3.0-or-later
// Map finite Geant4 material steps to the B0 tracking geometry.
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/ApproachDescriptor.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/VolumeBounds.hpp>
#include <Acts/Material/BinnedSurfaceMaterialAccumulater.hpp>
#include <Acts/Material/IntersectionMaterialAssigner.hpp>
#include <Acts/Material/MaterialMapper.hpp>
#include <Acts/Material/ProtoSurfaceMaterial.hpp>
#include "B0WindowTrackingGeometry.h"
#include <ActsPlugins/Json/MaterialMapJsonConverter.hpp>
#include <ActsPlugins/Root/RootMaterialTrackIo.hpp>
#include <DD4hep/Detector.h>
#include <DD4hep/VolumeManager.h>
#include <TChain.h>
#include <TGeoMatrix.h>
#include <TGeoShape.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Json = nlohmann::json;

// The navigation carrier has circular outer bounds; the DD4hep Boolean
// window has two apertures. Do not offer the carrier for mapping through holes.
class WindowAssignmentFinder final : public Acts::IAssignmentFinder {
public:
  WindowAssignmentFinder(Acts::IntersectionMaterialAssigner::Config config,
                         const Acts::Surface* windowSurface, dd4hep::DetElement window)
      : m_native(config, Acts::getDefaultLogger("Assignment", Acts::Logging::WARNING))
      , m_windowSurface(windowSurface)
      , m_window(std::move(window)) {}
  std::pair<std::vector<SurfaceAssignment>, std::vector<VolumeAssignment>>
  assignmentCandidates(const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
                       const Acts::Vector3& position,
                       const Acts::Vector3& direction) const override {
    auto candidates = m_native.assignmentCandidates(gctx, mctx, position, direction);
    std::erase_if(candidates.first, [&](const SurfaceAssignment& candidate) {
      if (candidate.surface != m_windowSurface) {
        return false;
      }
      const std::array global{candidate.position.x() / Acts::UnitConstants::cm,
                              candidate.position.y() / Acts::UnitConstants::cm,
                              candidate.position.z() / Acts::UnitConstants::cm};
      std::array<double, 3> local{};
      m_window.nominal().worldTransformation().MasterToLocal(global.data(), local.data());
      return !m_window.volume().solid().ptr()->Contains(local.data());
    });
    return candidates;
  }

private:
  Acts::IntersectionMaterialAssigner m_native;
  const Acts::Surface* m_windowSurface;
  dd4hep::DetElement m_window;
};

double trackingExit(const Acts::TrackingVolume& envelope, const Acts::GeometryContext& gctx,
                    const Acts::Vector3& origin, const Acts::Vector3& direction) {
  if (!envelope.inside(origin)) {
    throw std::runtime_error("Recorded ray starts outside tracking envelope");
  }
  double exit = std::numeric_limits<double>::infinity();
  for (const auto& boundary : envelope.boundarySurfaces()) {
    const auto hit = boundary->surfaceRepresentation()
                         .intersect(gctx, origin, direction, Acts::BoundaryTolerance::None())
                         .closestForward();
    if (hit.isValid() && hit.pathLength() > 0. &&
        !envelope.inside(hit.position() + 1e-4 * direction)) {
      exit = std::min(exit, hit.pathLength());
    }
  }
  if (!std::isfinite(exit)) {
    throw std::runtime_error("Cannot find first tracking-envelope exit");
  }
  return exit;
}

} // namespace

int main(int argc, char** argv) try {
  if (argc != 7) {
    std::cerr << "Usage: eicrecon-b0-material-map XML INPUT.root OUTPUT.cbor FIRST COUNT PAD_MM\n"
                 "Input must contain straight Geant4 geantino pre-step material records.\n";
    return 2;
  }
  const std::span arguments(argv, static_cast<std::size_t>(argc));
  const long long first = std::stoll(arguments[4]);
  const long long count = std::stoll(arguments[5]);
  const double padding  = std::stod(arguments[6]);
  if (!std::isfinite(padding) || padding <= 0.) {
    throw std::invalid_argument("Invalid layer padding");
  }
  if (std::filesystem::exists(arguments[3])) {
    throw std::runtime_error("Output already exists");
  }
  auto detector = dd4hep::Detector::make_unique("");
  detector->fromCompact(arguments[1]);
  detector->volumeManager();
  detector->apply("DD4hepVolumeManager", 0, nullptr);
  const Acts::GeometryContext gctx{};
  const Acts::MagneticFieldContext mctx{};
  auto logger         = Acts::getDefaultLogger("B0MaterialMapping", Acts::Logging::WARNING);
  const auto geometry = b0window::convertDD4hepDetector(
      detector->world(), *logger, Acts::equidistant, Acts::equidistant, Acts::equidistant, 1.,
      padding, Acts::UnitConstants::fm, ActsPlugins::sortDetElementsByID, gctx);
  std::set<const Acts::Surface*> navigationSurfaces;
  geometry->visitVolumes([&](const Acts::TrackingVolume* volume) {
    if (volume->confinedLayers()) {
      for (const auto& layer : volume->confinedLayers()->arrayObjects()) {
        if (layer->layerType() == Acts::navigation) {
          navigationSurfaces.insert(&layer->surfaceRepresentation());
        }
      }
    }
  });
  std::map<Acts::GeometryIdentifier, const Acts::Surface*> selected;
  auto select = [&](const Acts::Surface* surface) {
    if (surface->surfaceMaterial() && !navigationSurfaces.contains(surface)) {
      selected.emplace(surface->geometryId(), surface);
    }
  };
  geometry->visitSurfaces(select, false);
  geometry->visitVolumes([&](const Acts::TrackingVolume* volume) {
    for (const auto& boundary : volume->boundarySurfaces()) {
      select(&boundary->surfaceRepresentation());
    }
    if (volume->confinedLayers()) {
      for (const auto& layer : volume->confinedLayers()->arrayObjects()) {
        select(&layer->surfaceRepresentation());
        if (layer->approachDescriptor()) {
          for (const auto* surface : layer->approachDescriptor()->containedSurfaces()) {
            select(surface);
          }
        }
      }
    }
  });
  std::vector<const Acts::Surface*> surfaces;
  for (const auto& [id, surface] : selected) {
    surfaces.push_back(surface);
  }
  if (surfaces.empty()) {
    throw std::runtime_error("No material surfaces selected");
  }
  Acts::IntersectionMaterialAssigner::Config finderConfig;
  finderConfig.surfaces = surfaces;
  Acts::BinnedSurfaceMaterialAccumulater::Config accumulatorConfig;
  accumulatorConfig.geoContext         = gctx;
  accumulatorConfig.materialSurfaces   = surfaces;
  accumulatorConfig.emptyBinCorrection = true;
  const auto window                    = detector->detector("B0Window");
  std::vector<const Acts::Surface*> windowSurfaces;
  for (const auto* surface : surfaces) {
    if (b0window::isWindowSurface(*surface, window, gctx)) {
      windowSurfaces.push_back(surface);
    }
  }
  if (windowSurfaces.size() != 1) {
    throw std::runtime_error("Expected exactly one representing window material surface");
  }
  Acts::MaterialMapper::Config config;
  config.assignmentFinder =
      std::make_shared<WindowAssignmentFinder>(finderConfig, windowSurfaces.front(), window);
  config.surfaceMaterialAccumulater = std::make_shared<Acts::BinnedSurfaceMaterialAccumulater>(
      accumulatorConfig, Acts::getDefaultLogger("Accumulator", Acts::Logging::WARNING));
  Acts::MaterialMapper mapper(config, Acts::getDefaultLogger("Mapper", Acts::Logging::WARNING));
  auto state = mapper.createState();
  // The carrier is centered on the global axis, while the physical-window
  // binning is centered on its DD4hep placement. Native proto adjustment would
  // replace these bin bounds/transform with the carrier's annulus. Preserve the
  // explicit physical binning using the accumulator's public state interface.
  const auto* windowProto =
      dynamic_cast<const Acts::ProtoSurfaceMaterial*>(windowSurfaces.front()->surfaceMaterial());
  if (windowProto == nullptr) {
    throw std::runtime_error("Window mapping requires its explicit proto binning");
  }
  auto& accumulatorState = dynamic_cast<Acts::BinnedSurfaceMaterialAccumulater::State&>(
      *state->surfaceMaterialAccumulaterState);
  accumulatorState.accumulatedMaterial[windowSurfaces.front()->geometryId()] =
      Acts::AccumulatedSurfaceMaterial(windowProto->binning());
  TChain tree("material-tracks");
  if (tree.Add(arguments[2]) != 1) {
    throw std::runtime_error("Cannot open input tree");
  }
  const auto entries = tree.GetEntries();
  if (first < 0 || count <= 0 || first >= entries || count > entries - first) {
    throw std::invalid_argument("Entry range outside recorded tree");
  }
  for (const auto* name : {"event_id", "v_x", "v_y", "v_z", "v_px", "v_py", "v_pz", "mat_x",
                           "mat_y", "mat_z", "mat_dx", "mat_dy", "mat_dz", "mat_step_length",
                           "mat_X0", "mat_L0", "mat_A", "mat_Z", "mat_rho"}) {
    if (tree.GetBranch(name) == nullptr) {
      throw std::runtime_error(std::string("Missing branch ") + name);
    }
  }
  ActsPlugins::RootMaterialTrackIo reader({});
  reader.connectForRead(tree);
  double keptX0       = 0.;
  double assignedX0   = 0.;
  double unassignedX0 = 0.;
  for (long long entry = first; entry < first + count; ++entry) {
    if (tree.GetEntry(entry) <= 0) {
      throw std::runtime_error("Failed reading entry");
    }
    auto ray           = reader.read();
    const auto& origin = ray.first.first;
    if (!origin.allFinite() || !ray.first.second.allFinite() || ray.first.second.norm() <= 0.) {
      throw std::runtime_error("Invalid ray");
    }
    const Acts::Vector3 direction = ray.first.second.normalized();
    ray.first.second              = direction;
    const double exit = trackingExit(*geometry->highestTrackingVolume(), gctx, origin, direction);
    auto original     = std::move(ray.second.materialInteractions);
    ray.second.materialInteractions.clear();
    ray.second.materialInX0 = 0.;
    ray.second.materialInL0 = 0.;
    for (auto step : original) {
      const auto& material = step.materialSlab.material();
      const double length  = step.materialSlab.thickness();
      if (!step.position.allFinite() || !step.direction.allFinite() || !std::isfinite(length) ||
          length < 0. || step.direction.norm() <= 0. || !std::isfinite(material.X0()) ||
          material.X0() <= 0. || !std::isfinite(material.L0()) || material.L0() <= 0. ||
          !std::isfinite(material.Ar()) || material.Ar() <= 0. || !std::isfinite(material.Z()) ||
          material.Z() <= 0. || !std::isfinite(material.massDensity()) ||
          material.massDensity() <= 0.) {
        throw std::runtime_error("Nonfinite or invalid recorded material step");
      }
      const Acts::Vector3 stepDirection = step.direction.normalized();
      if (stepDirection.dot(direction) < 1. - 1e-6) {
        throw std::runtime_error("Input is not a straight geantino trajectory");
      }
      const double start = (step.position - origin).dot(direction);
      // ROOT stores single-precision coordinates. Permit their roundoff,
      // but reject displaced or curved trajectories instead of projecting them.
      if ((step.position - origin - start * direction).norm() >
          1e-6 * std::max(1., std::abs(start))) {
        throw std::runtime_error("Recorded step is not on the initial straight ray");
      }
      const double begin    = std::max(0., start);
      const double end      = std::min(exit, start + length);
      const double kept     = std::max(0., end - begin);
      const double fraction = length > 0. ? kept / length : 0.;
      keptX0 += fraction * step.materialSlab.thicknessInX0();
      if (kept <= 0.) {
        continue;
      }
      step.position += (begin - start + 0.5 * kept) * stepDirection;
      step.materialSlab = Acts::MaterialSlab(material, kept);
      ray.second.materialInX0 += step.materialSlab.thicknessInX0();
      ray.second.materialInL0 += step.materialSlab.thicknessInL0();
      ray.second.materialInteractions.push_back(std::move(step));
    }
    const auto [assigned, unassigned] = mapper.mapMaterial(*state, gctx, mctx, ray);
    assignedX0 += assigned.second.materialInX0;
    unassignedX0 += unassigned.second.materialInX0;
    if ((entry - first + 1) % 100000 == 0) {
      std::cout << "processed=" << entry - first + 1 << std::endl;
    }
  }
  const double residual = keptX0 - assignedX0 - unassignedX0;
  if (!std::isfinite(residual) || std::abs(residual) > 1e-6 * std::max(1., keptX0)) {
    throw std::runtime_error("Retained material not conserved in assignment");
  }
  Acts::MaterialMapJsonConverter converter({}, Acts::Logging::WARNING);
  const auto materialMap = converter.materialMapsToJson(mapper.finalizeMaps(*state));
  std::ofstream output(arguments[3], std::ios::binary);
  output.exceptions(std::ios::failbit | std::ios::badbit);
  Json::to_cbor(materialMap, output);
  output.close();
  std::cout << "mapped_entries=" << count << " surfaces=" << surfaces.size()
            << " retained_X0=" << keptX0 << " unassigned_X0=" << unassignedX0
            << " conservation_residual=" << residual << std::endl;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}
