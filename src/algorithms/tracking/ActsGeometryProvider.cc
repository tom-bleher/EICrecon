// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2024 Whitney Armstrong, Wouter Deconinck, Dmitry Romanov

#if Acts_VERSION_MAJOR < 45
#include <Acts/Geometry/DetectorElementBase.hpp>
#endif
#include <Acts/Geometry/GeometryIdentifier.hpp>
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
#include <ActsPlugins/DD4hep/DD4hepFieldAdapter.hpp>
#include <ActsPlugins/Json/JsonMaterialDecorator.hpp>
#include <ActsPlugins/Json/MaterialMapJsonConverter.hpp>
#include <DD4hep/DetElement.h>
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
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "ActsGeometryProvider.h"
#include "extensions/spdlog/SpdlogToActs.h"

template <typename T>
struct fmt::formatter<T, std::enable_if_t<std::is_base_of_v<Eigen::MatrixBase<T>, T>, char>>
    : fmt::ostream_formatter {};

/// @brief Material decorator wrapper that tracks per-layer material assignment
///
/// Wraps Acts::JsonMaterialDecorator and, for each decorate() call, records
/// whether material was assigned to any approach surface in each
/// (volume, layer) pair. After geometry conversion, call check() to emit
/// critical log messages for every layer that was visited but never had
/// material assigned to any of its approach surfaces.
class EpicJsonMaterialDecorator : public Acts::IMaterialDecorator {
public:
  /// Key identifying a layer: (volume id, layer id)
  using LayerKey = std::pair<Acts::GeometryIdentifier::Value, Acts::GeometryIdentifier::Value>;

  EpicJsonMaterialDecorator(const Acts::MaterialMapJsonConverter::Config& rConfig,
                            const std::string& jFileName, Acts::Logging::Level level,
                            std::shared_ptr<spdlog::logger> logger)
      : m_inner(rConfig, jFileName, level), m_log(std::move(logger)) {}

  void decorate(Acts::Surface& surface) const override {
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
  Acts::JsonMaterialDecorator m_inner;
  std::shared_ptr<spdlog::logger> m_log;
  /// All (volume, layer) pairs seen during decoration
  mutable std::set<LayerKey> m_decoratedLayers;
  /// Subset of decorated layers that had at least one surface with material
  mutable std::set<LayerKey> m_layersWithMaterial;
};

void ActsGeometryProvider::initialize(const dd4hep::Detector* dd4hep_geo, std::string material_file,
                                      std::shared_ptr<spdlog::logger> log,
                                      std::shared_ptr<spdlog::logger> init_log,
                                      bool require_material_coverage,
                                      std::vector<unsigned int> required_material_extra_ids) {
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

  std::set<unsigned int> requiredMaterialExtras;
  for (const auto value : required_material_extra_ids) {
    if (value > 0xffU) {
      throw std::invalid_argument("Required ACTS material detector ID exceeds the 8-bit extra field");
    }
    requiredMaterialExtras.insert(value);
  }
  if (require_material_coverage && requiredMaterialExtras.empty()) {
    throw std::invalid_argument(
        "Strict ACTS material coverage requires at least one explicit detector extra ID");
  }

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
    }
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

  // Keep the historical detector-wide diagnostic, even when strict selective
  // validation is disabled.
  if (auto epicDeco = std::dynamic_pointer_cast<const EpicJsonMaterialDecorator>(materialDeco)) {
    epicDeco->check();
  }

  // Optional production-validation contract. Check final ACTS surface IDs after
  // the DD4hep geometry-ID hook has run. The service supplies explicit detector
  // IDs; no detector name is hard-coded here.
  if (require_material_coverage) {
    if (!m_trackingGeo) {
      throw std::runtime_error("Strict ACTS material coverage requested without a tracking geometry");
    }
    std::map<unsigned int, std::size_t> requiredSurfaces;
    std::map<unsigned int, std::size_t> materialSurfaces;
    std::vector<Acts::GeometryIdentifier> missing;
    m_trackingGeo->visitSurfaces([&](const Acts::Surface* surface) {
      if (surface == nullptr) {
        return;
      }
      const auto id = surface->geometryId();
      const auto extra = static_cast<unsigned int>(id.extra());
      if (id.approach() == 0 || requiredMaterialExtras.find(extra) == requiredMaterialExtras.end()) {
        return;
      }
      ++requiredSurfaces[extra];
      if (surface->surfaceMaterial() != nullptr) {
        ++materialSurfaces[extra];
      } else {
        missing.push_back(id);
      }
    });

    for (const auto extra : requiredMaterialExtras) {
      const auto count = requiredSurfaces[extra];
      if (count == 0) {
        throw std::runtime_error(fmt::format(
            "Strict ACTS material coverage found no approach surfaces for detector extra ID {}",
            extra));
      }
      m_init_log->info("strict material coverage detector extra={} material={}/{} approach surfaces",
                       extra, materialSurfaces[extra], count);
    }
    for (const auto& id : missing) {
      m_init_log->critical(
          "Strict material coverage missing geometryId=(volume={}, boundary={}, layer={}, "
          "approach={}, sensitive={}, extra={})",
          id.volume(), id.boundary(), id.layer(), id.approach(), id.sensitive(), id.extra());
    }
    if (!missing.empty()) {
      throw std::runtime_error(fmt::format(
          "Strict ACTS material coverage failed: {} required approach surfaces have no material",
          missing.size()));
    }
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
    m_trackingGeo->visitSurfaces([this](const Acts::Surface* surface) {
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