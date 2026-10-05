// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/ApproachDescriptor.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Material/ProtoSurfaceMaterial.hpp>
#include <Acts/Surfaces/PlanarBounds.hpp>
#include <DD4hep/Alignments.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Readout.h>
#include <DD4hep/VolumeManager.h>
#include <DDRec/CellIDPositionConverter.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "algorithms/tracking/ActsGeometryProvider.h"
#include "algorithms/tracking/B0TripletSeedingConfig.h"
#include "algorithms/tracking/CKFTrackingConfig.h"

// Geometry/mapping consistency check for the B0 tracker (sensor surfaces,
// station grouping, readout round trips, material-map coverage). This is not
// a material-budget or end-to-end navigation validation. Without
// B0_GEOMETRY_TEST_COMPACT the program exits 77 so plain `ctest` skips it;
// the B0 validation workflow supplies all inputs and asserts pass/fail.
namespace {
void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

Acts::Vector3 actsPosition(const dd4hep::Position& point) {
  constexpr double scale = Acts::UnitConstants::mm / dd4hep::mm;
  return {point.x() * scale, point.y() * scale, point.z() * scale};
}

dd4hep::Position ddPosition(const Acts::Vector3& point) {
  constexpr double scale = dd4hep::mm / Acts::UnitConstants::mm;
  return {point.x() * scale, point.y() * scale, point.z() * scale};
}

bool volumeContains(const Acts::TrackingVolume& volume, const Acts::GeometryContext& gctx,
                    const Acts::Vector3& point, double tolerance) {
#if Acts_VERSION_MAJOR >= 45
  return volume.inside(gctx, point, tolerance);
#else
  (void)gctx;
  return volume.inside(point, tolerance);
#endif
}
} // namespace

int main(int argc, char** argv) {
  // Each argument may also come from the environment, so that CTest can run
  // without arguments. Without a compact file the test skips (exit 77).
  const char* envCompact  = std::getenv("B0_GEOMETRY_TEST_COMPACT");
  const char* envExpected = std::getenv("B0_GEOMETRY_TEST_EXPECTED");
  const char* envMaterial = std::getenv("B0_GEOMETRY_TEST_MATERIAL");
  if (argc < 1 || argc > 4) {
    std::cerr
        << "Usage: b0_geometry_test [compact.xml] [expected-sensors (0=any)] [material-map]\n";
    return 2;
  }
  const std::string compact = argc >= 2 ? argv[1] : (envCompact != nullptr ? envCompact : "");
  if (compact.empty()) {
    std::cout << "b0_geometry_test: skipped, no compact file "
                 "(pass an argument or set B0_GEOMETRY_TEST_COMPACT)\n";
    return 77;
  }
  try {
    const std::string expectedArg =
        argc >= 3 ? argv[2] : (envExpected != nullptr ? envExpected : "0");
    const std::size_t expected = std::stoul(expectedArg);
    const std::string material = argc == 4 ? argv[3] : (envMaterial != nullptr ? envMaterial : "");
    auto detector              = dd4hep::Detector::make_unique("");
    detector->fromCompact(compact);
    detector->apply("DD4hepVolumeManager", 0, nullptr);
    ActsGeometryProvider provider;
    auto logger = spdlog::default_logger();
    logger->set_level(spdlog::level::warn);
    provider.initialize(detector.get(), material, logger, logger);
    const auto& gctx        = provider.getActsGeometryContext();
    const auto readout      = detector->readout("B0TrackerHits");
    const auto segmentation = readout.segmentation();
    const dd4hep::rec::CellIDPositionConverter converter(*detector);
    std::set<dd4hep::VolumeID> sensitiveIDs;
    std::set<const Acts::Layer*> layers;
    std::vector<std::pair<const Acts::Layer*, double>> sensorLayers;
    constexpr double tolerance = 1.e-6 * Acts::UnitConstants::mm;

    auto checkElement = [&](auto&& self, dd4hep::DetElement element) -> void {
      if (element.volume().isSensitive()) {
        const auto volumeID     = element.volumeID();
        const std::string label = element.path();
        require(sensitiveIDs.insert(volumeID).second, label + ": duplicate sensitive volume ID");
        const auto found = provider.surfaceMap().find(volumeID);
        require(found != provider.surfaceMap().end(), label + ": missing ACTS surface");
        const auto* surface = found->second;
        const auto* layer   = surface->associatedLayer();
        require(layer != nullptr && layer->trackingVolume() != nullptr,
                label + ": ACTS surface has no owning tracking volume");
        layers.insert(layer);
        const auto* volume = layer->trackingVolume();
        const auto center  = surface->center(gctx);
        sensorLayers.emplace_back(layer, center.z());
        const auto alignment = element.nominal();
        require((actsPosition(alignment.localToWorld(dd4hep::Position{})) - center).norm() <
                    tolerance,
                label + ": DD4hep and ACTS sensor centers disagree");
        require(volumeContains(*volume, gctx, center, tolerance),
                label + ": sensor center outside ACTS volume");
        const auto* bounds = dynamic_cast<const Acts::PlanarBounds*>(&surface->bounds());
        require(bounds != nullptr, label + ": expected planar B0 sensor bounds");
        const auto normal = surface->normal(gctx, center, Acts::Vector3::UnitZ());
        auto vertices     = bounds->vertices();
        require(!vertices.empty(), label + ": sensor bounds have no vertices");
        vertices.push_back(Acts::Vector2::Zero());
        for (const auto& vertex : vertices) {
          const auto corner = surface->localToGlobal(gctx, vertex, normal);
          require(volumeContains(*volume, gctx, corner, tolerance),
                  label + ": sensor corner outside ACTS volume");

          // Sample inside every quadrant, avoiding edge pixels whose centers may
          // round outside the sensor. Use the configured segmentation axes.
          const auto sample   = surface->localToGlobal(gctx, 0.5 * vertex, normal);
          const auto local    = alignment.worldToLocal(ddPosition(sample));
          const auto cellID   = segmentation.cellID(local, ddPosition(sample), volumeID);
          const auto* context = converter.findContext(cellID);
          require(context != nullptr && context->identifier == volumeID,
                  label + ": pixel cellID does not resolve to its sensor");
          const auto pixel     = actsPosition(converter.position(cellID));
          const auto onSurface = surface->globalToLocal(gctx, pixel, normal, tolerance);
          require(onSurface.ok(), label + ": decoded pixel position is off its ACTS sensor plane");
          require((surface->localToGlobal(gctx, *onSurface, normal) - pixel).norm() < tolerance,
                  label + ": ACTS pixel coordinate roundtrip failed");
          const auto decodedLocal = segmentation.position(cellID);
          require((actsPosition(alignment.localToWorld(decodedLocal)) - pixel).norm() < tolerance,
                  label + ": DD4hep pixel coordinate roundtrip failed");
        }
      }
      for (const auto& [name, child] : element.children()) {
        self(self, child);
      }
    };
    checkElement(checkElement, detector->detector("B0Tracker"));
    require(!sensitiveIDs.empty(), "B0Tracker has no sensitive sensors");
    require(expected == 0 || sensitiveIDs.size() == expected, "Unexpected B0 sensor count");

    // The seeder groups hit positions and the CKF groups surface centers; both
    // must agree with the physical stations. Faces pair into stations with
    // the configured gap: within-station spread < gap < between-station gap.
    const double zGap = eicrecon::CKFTrackingConfig{}.measurementGroupZGap;
    require(eicrecon::B0TripletSeedingConfig{}.stationGap == zGap,
            "seeder and CKF station gaps differ");
    std::sort(sensorLayers.begin(), sensorLayers.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    std::vector<std::vector<std::pair<const Acts::Layer*, double>>> stations(1);
    for (const auto& entry : sensorLayers) {
      auto& current = stations.back();
      if (!current.empty() && entry.second - current.back().second > zGap) {
        stations.emplace_back();
      }
      stations.back().push_back(entry);
    }
    require(stations.size() == 4, "B0 sensors do not form four station groups");
    double maxSpread = 0.;
    double minGap    = std::numeric_limits<double>::max();
    for (std::size_t s = 0; s < stations.size(); ++s) {
      const auto& station = stations[s];
      std::set<const Acts::Layer*> stationLayers;
      for (const auto& [layer, z] : station) {
        static_cast<void>(z);
        stationLayers.insert(layer);
      }
      std::ostringstream label;
      label << "B0 station " << s;
      require(stationLayers.size() == 2, label.str() + ": expected a face pair");
      const double spread = station.back().second - station.front().second;
      require(spread < zGap, label.str() + ": within-station spread exceeds the gap");
      maxSpread = std::max(maxSpread, spread);
      if (s + 1 < stations.size()) {
        const double gap = stations[s + 1].front().second - station.back().second;
        require(gap > zGap, label.str() + ": between-station gap below the threshold");
        minGap = std::min(minGap, gap);
      }
    }

    // Every B0 layer needs mapped material on at least one approach surface,
    // and no approach surface may retain prototype material (which proves the
    // map does not resolve for this configuration). Individual surfaces
    // without any material stay allowed: the recorded map predates approach
    // surfaces added by newer ACTS layer builders.
    std::size_t materialSurfaces = 0;
    std::size_t concreteSurfaces = 0;
    if (!material.empty()) {
      std::vector<std::string> materialProblems;
      for (const auto* layer : layers) {
        const auto* approaches = layer->approachDescriptor();
        require(approaches != nullptr, "B0 layer has no approach surfaces");
        bool layerCovered = false;
        for (const auto* surface : approaches->containedSurfaces()) {
          const auto id = surface->geometryId();
          std::ostringstream tag;
          tag << "vol=" << id.volume() << "|lay=" << id.layer() << "|appr=" << id.approach();
          const auto* assigned = surface->surfaceMaterial();
          if (dynamic_cast<const Acts::ProtoSurfaceMaterial*>(assigned) != nullptr) {
            materialProblems.push_back(tag.str() + ": retains prototype material; "
                                                   "map does not match configuration");
          } else if (assigned != nullptr) {
            layerCovered = true;
            ++concreteSurfaces;
          }
          ++materialSurfaces;
        }
        if (!layerCovered) {
          materialProblems.push_back("layer without mapped material on any approach surface");
        }
      }
      require(materialSurfaces > 0, "No B0 material approach surfaces checked");
      if (!materialProblems.empty()) {
        std::ostringstream message;
        message << materialProblems.size() << " B0 material problems:";
        for (const auto& problem : materialProblems) {
          message << "\n  " << problem;
        }
        throw std::runtime_error(message.str());
      }
    }
    std::cout << "Validated " << sensitiveIDs.size() << " B0 sensors in " << layers.size()
              << " layers: center/corner containment and pixel coordinate roundtrips; "
              << stations.size() << " stations (spread " << maxSpread / Acts::UnitConstants::mm
              << " mm < " << zGap / Acts::UnitConstants::mm << " mm < "
              << minGap / Acts::UnitConstants::mm << " mm gap); " << concreteSurfaces << " of "
              << materialSurfaces << " approach surfaces carry mapped material.\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "B0 geometry validation failed: " << error.what() << '\n';
    return 1;
  }
}
