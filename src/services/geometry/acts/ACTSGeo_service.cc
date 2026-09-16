// Copyright 2022, David Lawrence
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//

#include "ACTSGeo_service.h"

#include <Acts/Visualization/ViewConfig.hpp>
#include <JANA/JApplication.h>
#include <JANA/JException.h>
#include <JANA/Services/JServiceLocator.h>
#include <array>
#include <cctype>
#include <exception>
#include <gsl/pointers>
#include <stdexcept>
#include <string>
#include <vector>

#include "ActsGeometryProvider.h"
#include "MaterialMapContract.h"
#include "services/geometry/dd4hep/DD4hep_service.h"
#include "services/log/Log_service.h"

namespace {
std::vector<std::string> splitNames(const std::string& input) {
  std::vector<std::string> result;
  std::size_t start = 0;
  while (start <= input.size()) {
    const auto comma = input.find(',', start);
    auto value = input.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first != std::string::npos) {
      const auto last = value.find_last_not_of(" \t\r\n");
      result.push_back(value.substr(first, last - first + 1));
    }
    if (comma == std::string::npos) {
      break;
    }
    start = comma + 1;
  }
  return result;
}
} // namespace

// Virtual destructor implementation to pin vtable and typeinfo to this
// translation unit
ACTSGeo_service::~ACTSGeo_service() = default;

//----------------------------------------------------------------
// detector
//
/// Return pointer to the dd4hep::Detector object.
/// Call Initialize if needed.
//----------------------------------------------------------------
std::shared_ptr<const ActsGeometryProvider> ACTSGeo_service::actsGeoProvider() {

  try {
    std::call_once(m_init_flag, [this]() {
      // Assemble everything on the first call

      if (m_dd4hepGeo == nullptr) {
        throw JException("ACTSGeo_service m_dd4hepGeo==null which should never be!");
      }

      // Resolve the geometry-declared/default material map separately from the
      // final selected value so that a JANA override cannot hide its provenance.
      std::string geometry_declared_material_map;
      bool geometry_declared_material_map_present = false;
      try {
        geometry_declared_material_map = m_dd4hepGeo->constant<std::string>("material-map");
        geometry_declared_material_map_present = true;
      } catch (const std::runtime_error&) {
        geometry_declared_material_map.clear();
      }

      const std::string fallback_material_map = "calibrations/materials-map.cbor";
      std::string material_map_file = geometry_declared_material_map_present
                                          ? geometry_declared_material_map
                                          : fallback_material_map;
      m_app->SetDefaultParameter("acts:MaterialMap", material_map_file,
                                 "JSON/CBOR material map file path");

      bool require_geometry_material_map = false;
      m_app->SetDefaultParameter(
          "acts:RequireGeometryMaterialMap", require_geometry_material_map,
          "Require selected material-map content to match the map declared by DD4hep geometry");

      std::string expected_material_map_sha256;
      try {
        expected_material_map_sha256 =
            m_dd4hepGeo->constant<std::string>("material-map-sha256");
      } catch (const std::runtime_error&) {
        expected_material_map_sha256.clear();
      }
      m_app->SetDefaultParameter(
          "acts:ExpectedMaterialMapSHA256", expected_material_map_sha256,
          "Optional expected SHA-256 for the selected and geometry-declared material map");

      bool require_material_coverage = false;
      m_app->SetDefaultParameter(
          "acts:RequireMaterialCoverage", require_material_coverage,
          "Fail initialization when required detector approach surfaces lack ACTS material");
      std::string required_material_detector_constants;
      m_app->SetDefaultParameter(
          "acts:RequiredMaterialDetectorConstants", required_material_detector_constants,
          "Comma-separated DD4hep constant names whose detector IDs require complete approach-surface material");

      std::vector<unsigned int> required_material_extra_ids;
      for (const auto& name : splitNames(required_material_detector_constants)) {
        try {
          required_material_extra_ids.push_back(
              static_cast<unsigned int>(m_dd4hepGeo->constant<unsigned long>(name) & 0xffUL));
        } catch (const std::exception& error) {
          throw std::runtime_error("Cannot resolve required material detector constant '" + name +
                                   "': " + error.what());
        }
      }

      const eicrecon::acts_material::MaterialMapSelection selection{
          .geometryDeclaredPath = geometry_declared_material_map,
          .selectedPath = material_map_file,
          .expectedSha256 = expected_material_map_sha256,
          .geometryDeclared = geometry_declared_material_map_present,
          .usedFallback = !geometry_declared_material_map_present,
          .requireGeometryMaterialMap = require_geometry_material_map,
      };
      const auto validation = eicrecon::acts_material::validateMaterialMapSelection(selection);
      const std::string selection_source =
          material_map_file != (geometry_declared_material_map_present
                                    ? geometry_declared_material_map
                                    : fallback_material_map)
              ? "parameter override"
              : (geometry_declared_material_map_present ? "geometry declaration" : "generic fallback");
      m_log->info(
          "ACTS material-map contract: declared='{}' selected='{}' source='{}' fallback={} "
          "requireGeometryMatch={} expectedSHA256='{}' selectedSHA256='{}'",
          geometry_declared_material_map_present ? geometry_declared_material_map : "<none>",
          material_map_file, selection_source, !geometry_declared_material_map_present,
          require_geometry_material_map,
          validation.expectedSha256.empty() ? "<none>" : validation.expectedSha256,
          validation.selectedSha256.empty() ? "<not checked>" : validation.selectedSha256);

      if (require_material_coverage) {
        if (required_material_extra_ids.empty()) {
          throw std::invalid_argument(
              "acts:RequireMaterialCoverage=true requires acts:RequiredMaterialDetectorConstants");
        }
        m_log->info("ACTS strict material coverage enabled for {} detector ID(s)",
                    required_material_extra_ids.size());
      }

      // Reading the geometry may take a long time and if the JANA ticker is enabled, it will keep printing
      // while no other output is coming which makes it look like something is wrong. Disable the ticker
      // while parsing and loading the geometry
      auto tickerEnabled = m_app->IsTickerEnabled();
      m_app->SetTicker(false);

      // Create default m_acts_provider
      m_acts_provider = std::make_shared<ActsGeometryProvider>();

      // Set ActsGeometryProvider parameters
      bool objWriteIt = m_acts_provider->getObjWriteIt();
      bool plyWriteIt = m_acts_provider->getPlyWriteIt();
      m_app->SetDefaultParameter("acts:WriteObj", objWriteIt,
                                 "Write tracking geometry as obj files");
      m_app->SetDefaultParameter("acts:WritePly", plyWriteIt,
                                 "Write tracking geometry as ply files");
      m_acts_provider->setObjWriteIt(objWriteIt);
      m_acts_provider->setPlyWriteIt(plyWriteIt);

      double layerEnvelopeR = m_acts_provider->getLayerEnvelopeR();
      double layerEnvelopeZ = m_acts_provider->getLayerEnvelopeZ();
      m_app->SetDefaultParameter("acts:LayerEnvelopeR", layerEnvelopeR,
                                 "Radial pad added to the bounds of each subdetector tracking "
                                 "volume beyond its outermost layers [mm]");
      m_app->SetDefaultParameter("acts:LayerEnvelopeZ", layerEnvelopeZ,
                                 "Longitudinal pad added to the z bounds of each subdetector "
                                 "tracking volume beyond its outermost layers [mm]; must cover "
                                 "r*tan(tilt) of layers tilted off the beam axis and stay below "
                                 "half the z gap to the neighbouring volume");
      m_acts_provider->setLayerEnvelopeR(layerEnvelopeR);
      m_acts_provider->setLayerEnvelopeZ(layerEnvelopeZ);

      std::string outputTag = m_acts_provider->getOutputTag();
      std::string outputDir = m_acts_provider->getOutputDir();
      m_app->SetDefaultParameter("acts:OutputTag", outputTag, "Obj and ply output file tag");
      m_app->SetDefaultParameter("acts:OutputDir", outputDir, "Obj and ply output file dir");
      m_acts_provider->setOutputTag(outputTag);
      m_acts_provider->setOutputDir(outputDir);

      std::array<int, 3> containerView = m_acts_provider->getContainerView().color.rgb;
      std::array<int, 3> volumeView    = m_acts_provider->getVolumeView().color.rgb;
      std::array<int, 3> sensitiveView = m_acts_provider->getSensitiveView().color.rgb;
      std::array<int, 3> passiveView   = m_acts_provider->getPassiveView().color.rgb;
      std::array<int, 3> gridView      = m_acts_provider->getGridView().color.rgb;
      m_app->SetDefaultParameter("acts:ContainerView", containerView, "RGB for container views");
      m_app->SetDefaultParameter("acts:VolumeView", volumeView, "RGB for volume views");
      m_app->SetDefaultParameter("acts:SensitiveView", sensitiveView, "RGB for sensitive views");
      m_app->SetDefaultParameter("acts:PassiveView", passiveView, "RGB for passive views");
      m_app->SetDefaultParameter("acts:GridView", gridView, "RGB for grid views");
      m_acts_provider->setContainerView(containerView);
      m_acts_provider->setVolumeView(volumeView);
      m_acts_provider->setSensitiveView(sensitiveView);
      m_acts_provider->setPassiveView(passiveView);
      m_acts_provider->setGridView(gridView);

      // Initialize m_acts_provider
      m_acts_provider->initialize(m_dd4hepGeo, material_map_file, m_log, m_log,
                                  require_material_coverage, required_material_extra_ids);

      // Enable ticker back
      m_app->SetTicker(tickerEnabled);
    });
  } catch (std::exception& ex) {
    throw JException(ex.what());
  }

  return m_acts_provider;
}

void ACTSGeo_service::acquire_services(JServiceLocator* srv_locator) {

  auto log_service = srv_locator->get<Log_service>();
  m_log            = log_service->logger("acts");

  // DD4Hep geometry
  auto dd4hep_service = srv_locator->get<DD4hep_service>();
  m_dd4hepGeo         = dd4hep_service->detector();
}