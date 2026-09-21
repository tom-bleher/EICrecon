// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024 Whitney Armstrong, Wouter Deconinck, Dmitry Romanov

#if Acts_VERSION_MAJOR < 45
#include <Acts/Geometry/DetectorElementBase.hpp>
#endif
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Material/IMaterialDecorator.hpp>
#include <fmt/format.h>
#if Acts_VERSION_MAJOR >= 45
#include <Acts/Surfaces/SurfacePlacementBase.hpp>
#endif
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/BinningType.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ObjVisualization3D.hpp>
#include <Acts/Visualization/PlyVisualization3D.hpp>
#include <ActsPlugins/DD4hep/ConvertDD4hepDetector.hpp>
#include <ActsPlugins/DD4hep/DD4hepDetectorElement.hpp>
#include <ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp>
#include <ActsPlugins/DD4hep/DD4hepFieldAdapter.hpp>
#include <ActsPlugins/Json/JsonMaterialDecorator.hpp>
#include <ActsPlugins/Json/MaterialMapJsonConverter.hpp>
#include <DD4hep/DetElement.h>
#include <DD4hep/Alignments.h>
#include <DD4hep/DetectorTools.h>
#include <DD4hep/VolumeManager.h>
#include <TGeoManager.h>
#include <boost/container/detail/std_fwd.hpp>
#include <fmt/ostream.h>
#include <spdlog/common.h>
// Formatter for Eigen matrices
#include <Eigen/Core>
#include <exception>
#include <filesystem>
#include <functional>
#include <initializer_list>
#include <set>
#include <sstream>
#include <type_traits>
#include <utility>

#include "ActsGeometryProvider.h"
#include "B0MaterialValidation.h"
#include "extensions/spdlog/SpdlogToActs.h"

template <typename T>
struct fmt::formatter<T, std::enable_if_t<std::is_base_of_v<Eigen::MatrixBase<T>, T>, char>>
    : fmt::ostream_formatter {};

/// Preserve geometry-declared mapping targets before JSON decoration, including
/// passive representing surfaces. Also report layers without approach material.
class EpicJsonMaterialDecorator : public Acts::IMaterialDecorator {
public:
  /// Key identifying a layer: (volume id, layer id)
  using LayerKey = std::pair<Acts::GeometryIdentifier::Value, Acts::GeometryIdentifier::Value>;

  EpicJsonMaterialDecorator(const Acts::MaterialMapJsonConverter::Config& rConfig,
                            const std::string& jFileName, Acts::Logging::Level level,
                            std::shared_ptr<spdlog::logger> logger)
      : m_inner(rConfig, jFileName, level), m_log(std::move(logger)) {}

  void decorate(Acts::Surface& surface) const override {
    // Preserve the XML mapping targets: missing JSON entries leave prototypes
    // in place, which are not physical material usable by the fitter.
    if (eicrecon::isPrototypeMaterial(surface.surfaceMaterial())) {
      m_expectedMaterial.insert(surface.geometryId());
      // Navigation-layer representing surfaces receive their final identifier
      // after this callback. Preserve identity by pointer until closure ends.
      m_expectedSurfaces.push_back(&surface);
    }
    m_inner.decorate(surface);
    const auto id = surface.geometryId();
    m_log->trace("{} assigned to surface with geometryId=(volume={}, boundary={}, layer={}, "
                 "approach={}, sensitive={}, extra={})",
                 (surface.surfaceMaterial() != nullptr) ? "Material" : "No material", id.volume(),
                 id.boundary(), id.layer(), id.approach(), id.sensitive(), id.extra());
    // Only consider approach surfaces
    if (id.approach() == 0) {
      return;
    }
    LayerKey key{id.volume(), id.layer()};
    // Record that this layer was visited
    m_decoratedLayers.insert(key);
    // If material was assigned, record that
    if (surface.surfaceMaterial() != nullptr) {
      m_layersWithMaterial.insert(key);
    }
  }

  void decorate(Acts::TrackingVolume& volume) const override { m_inner.decorate(volume); }

  const std::set<Acts::GeometryIdentifier>& expectedMaterial() const { return m_expectedMaterial; }

  void validateB0Volumes(const std::set<Acts::GeometryIdentifier::Value>& volumes) const {
    eicrecon::validateB0MaterialTargets(eicrecon::b0MaterialTargetsAfterClosure(m_expectedSurfaces),
                                        volumes);
  }

  void validateExternalPassiveLayers(const dd4hep::DetElement& detector,
                                     const std::set<Acts::GeometryIdentifier::Value>& volumes,
                                     const Acts::GeometryContext& gctx) const {
    std::istringstream references(
        ActsPlugins::getParamOr<std::string>("passive_layer_references", detector, ""));
    std::string path;
    while (references >> path) {
      const auto element = dd4hep::detail::tools::findDaughterElement(detector.world(), path);
      if (!element.isValid() || !ActsPlugins::getParamOr<bool>("passive_disc", element, false)) {
        throw std::runtime_error("Invalid B0Tracker external passive disc: " + path);
      }
      const auto& parameters     = ActsPlugins::getParams(element);
      const auto& world          = element.nominal().worldTransformation();
      const auto* rotation       = world.GetRotationMatrix();
      const auto* translation    = world.GetTranslation();
      Acts::Transform3 transform = Acts::Transform3::Identity();
      for (int i = 0; i < 3; ++i) {
        transform.translation()[i] = translation[i] * Acts::UnitConstants::cm;
        for (int j = 0; j < 3; ++j) {
          transform.linear()(i, j) = rotation[3 * i + j];
        }
      }
      eicrecon::validateB0PassiveDiscRepresentation(
          eicrecon::b0MaterialTargetsAfterClosure(m_expectedSurfaces), volumes, gctx, transform,
          parameters.get<double>("passive_disc_r_min"),
          parameters.get<double>("passive_disc_r_max"), path);
    }
  }

  /// Report every decorated layer that never received material on any approach surface.
  void check() const {
    for (const auto& key : m_decoratedLayers) {
      if (m_layersWithMaterial.find(key) == m_layersWithMaterial.end()) {
        m_log->critical(
            "No material assigned to any approach surface in layer (volume={}, layer={})",
            key.first, key.second);
      }
    }
  }

private:
  mutable std::set<Acts::GeometryIdentifier> m_expectedMaterial;
  mutable std::vector<const Acts::Surface*> m_expectedSurfaces;
  Acts::JsonMaterialDecorator m_inner;
  std::shared_ptr<spdlog::logger> m_log;
  /// All (volume, layer) pairs seen during decoration
  mutable std::set<LayerKey> m_decoratedLayers;
  /// Subset of decorated layers that had at least one surface with material
  mutable std::set<LayerKey> m_layersWithMaterial;
};

void ActsGeometryProvider::initialize(const dd4hep::Detector* dd4hep_geo, std::string material_file,
                                      std::shared_ptr<spdlog::logger> log,
                                      std::shared_ptr<spdlog::logger> init_log) {
  // LOGGING
  m_log      = log;
  m_init_log = init_log;

  m_init_log->debug("ActsGeometryProvider initializing...");

  m_init_log->debug("Set TGeoManager and acts_init_log_level log levels");
  // Turn off TGeo printouts if appropriate for the msg level
  if (m_log->level() >= (int)spdlog::level::info) {
    TGeoManager::SetVerboseLevel(0);
  }

  // Set ACTS logging level
  auto acts_init_log_level = eicrecon::SpdlogToActsLevel(m_init_log->level());

  m_dd4hepDetector = dd4hep_geo;

  // Load ACTS materials maps
  std::shared_ptr<const Acts::IMaterialDecorator> materialDeco{nullptr};
  if (!material_file.empty()) {
    m_init_log->info("loading materials map from file: '{}'", material_file);
    // Set up the converter first
    Acts::MaterialMapJsonConverter::Config jsonGeoConvConfig;
    // Set up the json-based decorator, wrapped to report undecorated layers
    materialDeco = std::make_shared<const EpicJsonMaterialDecorator>(
        jsonGeoConvConfig, material_file, acts_init_log_level, m_init_log);
  }

  // Geometry identifier hook to write detector ID to extra field
  class ConvertDD4hepDetectorGeometryIdentifierHook : public Acts::GeometryIdentifierHook {
    Acts::GeometryIdentifier decorateIdentifier(Acts::GeometryIdentifier identifier,
                                                const Acts::Surface& surface) const override {
#if Acts_VERSION_MAJOR >= 45
      const auto* placement = surface.surfacePlacement();
      const auto* dd4hep_det_element =
          dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(placement);
#else
      const auto* dd4hep_det_element = dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(
          surface.associatedDetectorElement());
#endif
      if (dd4hep_det_element == nullptr) {
        return identifier;
      }
      // set 8-bit extra field to 8-bit DD4hep detector ID
      return identifier.withExtra(0xff & dd4hep_det_element->identifier());
    };
  };
  auto geometryIdHook = std::make_shared<ConvertDD4hepDetectorGeometryIdentifierHook>();

  // Convert DD4hep geometry to ACTS
  m_init_log->info("Converting DD4Hep geometry to ACTS...");
  auto logger                  = eicrecon::getSpdlogLogger("CONV", m_log);
  Acts::BinningType bTypePhi   = Acts::equidistant;
  Acts::BinningType bTypeR     = Acts::equidistant;
  Acts::BinningType bTypeZ     = Acts::equidistant;
  double layerEnvelopeR        = m_layerEnvelopeR * Acts::UnitConstants::mm;
  double layerEnvelopeZ        = m_layerEnvelopeZ * Acts::UnitConstants::mm;
  double defaultLayerThickness = Acts::UnitConstants::fm;

  try {
    m_trackingGeo = ActsPlugins::convertDD4hepDetector(
        m_dd4hepDetector->world(), *logger, bTypePhi, bTypeR, bTypeZ, layerEnvelopeR,
        layerEnvelopeZ, defaultLayerThickness, ActsPlugins::sortDetElementsByID, m_trackingGeoCtx,
        materialDeco, geometryIdHook);
  } catch (const std::exception& ex) {
    m_init_log->error("Error during DD4Hep -> ACTS geometry conversion: {}", ex.what());
    m_init_log->info("Set parameter acts:LogLevel=trace to see conversion info and possibly "
                     "identify failing geometry");
    throw;
  }

  m_init_log->info("DD4Hep geometry converted!");

  // Report layers that were visited but never had material assigned
  if (auto epicDeco = std::dynamic_pointer_cast<const EpicJsonMaterialDecorator>(materialDeco)) {
    epicDeco->check();
  }

  // Visit surfaces
  m_init_log->info("Checking surfaces...");
  if (m_trackingGeo) {
    // Write tracking geometry to collection of obj or ply files
    const Acts::TrackingVolume* world = m_trackingGeo->highestTrackingVolume();
    if (m_objWriteIt) {
      m_init_log->info("Writing obj files to {}...", m_outputDir);
      Acts::ObjVisualization3D objVis;
      Acts::GeometryView3D::drawTrackingVolume(objVis, *world, m_trackingGeoCtx, m_containerView,
                                               m_volumeView, m_passiveView, m_sensitiveView,
                                               m_gridView, m_objWriteIt, m_outputTag, m_outputDir);
    }
    if (m_plyWriteIt) {
      m_init_log->info("Writing ply files to {}...", m_outputDir);
      Acts::PlyVisualization3D plyVis;
      Acts::GeometryView3D::drawTrackingVolume(plyVis, *world, m_trackingGeoCtx, m_containerView,
                                               m_volumeView, m_passiveView, m_sensitiveView,
                                               m_gridView, m_plyWriteIt, m_outputTag, m_outputDir);
    }

    m_init_log->debug("visiting all the surfaces  ");
    const auto epicDeco = std::dynamic_pointer_cast<const EpicJsonMaterialDecorator>(materialDeco);
    std::set<Acts::GeometryIdentifier> checkedB0Layers;
    m_trackingGeo->visitSurfaces([this, &epicDeco, &checkedB0Layers](const Acts::Surface* surface) {
      // for now we just require a valid surface
      if (surface == nullptr) {
        m_init_log->info("no surface??? ");
        return;
      }
#if Acts_VERSION_MAJOR >= 45
      const auto* placement   = surface->surfacePlacement();
      const auto* det_element = dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(placement);
#else
      const auto* det_element = dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(
          surface->associatedDetectorElement());
#endif

      if (det_element == nullptr) {
        m_init_log->error("invalid det_element!!! det_element == nullptr ");
        return;
      }

      // Only validate reconstruction with a supplied map. Geometry construction
      // without a decorator is also used to generate new material maps.
      if (epicDeco) {
        for (auto source = det_element->sourceElement(); source.isValid();
             source      = source.parent()) {
          if (std::string(source.name()) != "B0Tracker") {
            continue;
          }
          const auto* layer = surface->associatedLayer();
          if (layer == nullptr || layer->approachDescriptor() == nullptr) {
            throw std::runtime_error("B0Tracker sensitive surface has no material approach layer");
          }
          if (checkedB0Layers.insert(layer->geometryId()).second) {
            eicrecon::validateB0ApproachMaterial(*layer->approachDescriptor(),
                                                 epicDeco->expectedMaterial());
          }
          break;
        }
      }

      // more verbose output is lower enum value
      m_init_log->debug(" det_element->identifier() = {} ", det_element->identifier());
      auto volman   = m_dd4hepDetector->volumeManager();
      auto* vol_ctx = volman.lookupContext(det_element->identifier());
      auto vol_id   = vol_ctx->identifier;

      if (m_init_log->level() <= spdlog::level::debug) {
        auto de = vol_ctx->element;
        m_init_log->debug("  de.path          = {}", de.path());
        m_init_log->debug("  de.placementPath = {}", de.placementPath());
      }

      this->m_surfaces.insert_or_assign(vol_id, surface);
    });
    if (epicDeco) {
      // Passive layers have no sensitive surface to visit above. Validate every
      // XML-declared mapping target in the discovered B0 tracking volumes too.
      std::set<Acts::GeometryIdentifier::Value> b0Volumes;
      for (const auto& layer : checkedB0Layers) {
        b0Volumes.insert(layer.volume());
      }
      if (!b0Volumes.empty()) {
        epicDeco->validateExternalPassiveLayers(m_dd4hepDetector->detector("B0Tracker"), b0Volumes,
                                                m_trackingGeoCtx);
      }
      epicDeco->validateB0Volumes(b0Volumes);
    }
  } else {
    m_init_log->error("m_trackingGeo==null why am I still alive???");
  }

  // Load ACTS magnetic field
  m_init_log->info("Loading magnetic field...");
  m_magneticField = std::make_shared<ActsPlugins::DD4hepFieldAdapter>(m_dd4hepDetector->field());
  auto bCache     = m_magneticField->makeCache(Acts::MagneticFieldContext{});
  for (int z : {0, 500, 1000, 1500, 2000, 3000, 4000}) {
    auto b = m_magneticField->getField({0.0, 0.0, double(z)}, bCache).value();
    m_init_log->debug("B(z = {:>5} [mm]) = {} T", z, b.transpose() / Acts::UnitConstants::T);
  }

  m_init_log->info("ActsGeometryProvider initialization complete");
}

std::shared_ptr<const Acts::MagneticFieldProvider> ActsGeometryProvider::getFieldProvider() const {
  return m_magneticField;
}
