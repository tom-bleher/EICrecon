// SPDX-License-Identifier: MPL-2.0
// Geometry orchestration follows ACTS ConvertDD4hepDetector.cpp (MPL-2.0).
// Copyright (C) 2016 CERN for the benefit of the ACTS project.
#include "B0WindowTrackingGeometry.h"
#include <Acts/Geometry/DiscLayer.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/ILayerBuilder.hpp>
#include <Acts/Geometry/TrackingGeometryBuilder.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Material/ProtoSurfaceMaterial.hpp>
#include <Acts/Material/IMaterialDecorator.hpp>
#include <Acts/Material/ISurfaceMaterial.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <Acts/Surfaces/SurfaceArray.hpp>
#include <Acts/Utilities/BinUtility.hpp>
#include <DD4hep/Detector.h>
#include <TGeoBBox.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoTube.h>
#include <TGeoMatrix.h>
#include <algorithm>
#include <array>
#include <span>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace b0window {
bool isWindowSurface(const Acts::Surface& surface, dd4hep::DetElement window,
                     const Acts::GeometryContext& context) {
  const auto* layer  = surface.associatedLayer();
  const auto* bounds = dynamic_cast<const Acts::RadialBounds*>(&surface.bounds());
  const std::span<const double, 3> position(window.nominal().worldTransformation().GetTranslation(),
                                            3);
  const double z = position[2] * Acts::UnitConstants::cm;
  return (layer != nullptr) && layer->layerType() == Acts::passive && (bounds != nullptr) &&
         (surface.center(context) - Acts::Vector3(0., 0., z)).norm() < 1e-6;
}
namespace {
  dd4hep::DetElement findWindow(dd4hep::DetElement element) {
    if (std::string(element.name()) == "B0Window") {
      return element;
    }
    for (const auto& [name, child] : element.children()) {
      auto found = findWindow(child);
      if (found.isValid()) {
        return found;
      }
    }
    return {};
  }
  // Keep the physical apertures exact even when a material-map bin spans an edge.
  class WindowMaterial final : public Acts::ISurfaceMaterial {
  public:
    WindowMaterial(std::shared_ptr<const Acts::ISurfaceMaterial> base, dd4hep::DetElement window)
        : Acts::ISurfaceMaterial(
              base->factor(Acts::Direction::Negative(), Acts::MaterialUpdateStage::PreUpdate),
              base->mappingType())
        , m_base(std::move(base))
        , m_shape(window.volume().solid().ptr())
        , m_matrix(window.nominal().worldTransformation()) {}
    Acts::ISurfaceMaterial& scale(double factor) override {
      if (factor != 1.) {
        throw std::logic_error("Window prototype material scaling is unsupported");
      }
      return *this;
    }
    const Acts::MaterialSlab& materialSlab(const Acts::Vector2& lp) const override {
      // Carrier coordinates are global-axis (r, phi); map bins use the actual
      // offset window frame, so resolve both mask and material through global.
      const std::span<const double, 3> position(m_matrix.GetTranslation(), 3);
      return materialSlab(Acts::Vector3(lp[0] * std::cos(lp[1]), lp[0] * std::sin(lp[1]),
                                        position[2] * Acts::UnitConstants::cm));
    }
    const Acts::MaterialSlab& materialSlab(const Acts::Vector3& gp) const override {
      const std::array<double, 3> global = {gp.x() / Acts::UnitConstants::cm,
                                            gp.y() / Acts::UnitConstants::cm,
                                            gp.z() / Acts::UnitConstants::cm};
      std::array<double, 3> local{};
      m_matrix.MasterToLocal(global.data(), local.data());
      local[2] = 0.;
      return m_shape->Contains(local.data()) ? m_base->materialSlab(gp) : m_vacuum;
    }
    std::ostream& toStream(std::ostream& stream) const override {
      return stream << "B0 window aperture mask: " << *m_base;
    }

  private:
    std::shared_ptr<const Acts::ISurfaceMaterial> m_base;
    TGeoShape* m_shape;
    TGeoHMatrix m_matrix;
    Acts::MaterialSlab m_vacuum = Acts::MaterialSlab::Nothing();
  };
  class WindowDecorator final : public Acts::IMaterialDecorator {
  public:
    WindowDecorator(std::shared_ptr<const Acts::IMaterialDecorator> base, dd4hep::DetElement window)
        : m_base(std::move(base)), m_window(window) {}
    void decorate(Acts::TrackingVolume& volume) const override { m_base->decorate(volume); }
    void decorate(Acts::Surface& surface) const override {
      m_base->decorate(surface);
      if (isWindowSurface(surface, m_window, Acts::GeometryContext())) {
        const auto* material = surface.surfaceMaterial();
        if ((material == nullptr) ||
            (dynamic_cast<const Acts::ProtoSurfaceMaterial*>(material) != nullptr)) {
          throw std::runtime_error(
              "B0 window material is missing; acts:B0WindowMaterial requires a matching "
              "regenerated material map");
        }
        surface.assignSurfaceMaterial(
            std::make_shared<WindowMaterial>(surface.surfaceMaterialSharedPtr(), m_window));
      }
    }

  private:
    std::shared_ptr<const Acts::IMaterialDecorator> m_base;
    dd4hep::DetElement m_window;
  };
  class WindowVolumeHelper final : public Acts::CylinderVolumeHelper {
  public:
    WindowVolumeHelper(const Config& config, std::unique_ptr<const Acts::Logger> logger,
                       dd4hep::DetElement window, double padding,
                       std::shared_ptr<unsigned int> inserted)
        : Acts::CylinderVolumeHelper(config, logger->clone("B0WindowVolumeHelper"))
        , m_window(window)
        , m_padding(padding)
        , m_inserted(std::move(inserted))
        , m_logger(std::move(logger)) {}
    Acts::MutableTrackingVolumePtr
    createGapTrackingVolume(const Acts::GeometryContext& context,
                            Acts::MutableTrackingVolumeVector& volumes,
                            std::shared_ptr<const Acts::IVolumeMaterial> material, double rMin,
                            double rMax, double zMin, double zMax, unsigned int layers,
                            bool cylinder, const std::string& name) const override {
      if (name != "acts_beampipe_central::fGap") {
        return Acts::CylinderVolumeHelper::createGapTrackingVolume(
            context, volumes, std::move(material), rMin, rMax, zMin, zMax, layers, cylinder, name);
      }
      if (cylinder || layers != 1 || !volumes.empty() || material) {
        throw std::runtime_error("Window insertion requires the verified empty axial gap");
      }
      const auto& matrix        = m_window.nominal().worldTransformation();
      Acts::Transform3 physical = Acts::Transform3::Identity();
      const std::span<const double, 9> rotation(matrix.GetRotationMatrix(), 9);
      const std::span<const double, 3> position(matrix.GetTranslation(), 3);
      for (int i = 0; i < 3; ++i) {
        physical.translation()[i] = position[i] * Acts::UnitConstants::cm;
        for (int j = 0; j < 3; ++j) {
          physical.linear()(i, j) = rotation[3 * i + j];
        }
      }
      if (!physical.linear().isApprox(Acts::RotationMatrix3::Identity())) {
        throw std::runtime_error("Window gap prototype requires the actual unrotated window");
      }
      auto* box      = dynamic_cast<TGeoBBox*>(m_window.volume().solid().ptr());
      auto* outer    = dynamic_cast<TGeoCompositeShape*>(m_window.volume().solid().ptr());
      auto* inner    = (outer != nullptr)
                           ? dynamic_cast<TGeoCompositeShape*>(outer->GetBoolNode()->GetLeftShape())
                           : nullptr;
      auto* aperture = (inner != nullptr)
                           ? dynamic_cast<TGeoTube*>(inner->GetBoolNode()->GetRightShape())
                           : nullptr;
      if ((box == nullptr) || (aperture == nullptr)) {
        throw std::runtime_error("Unexpected physical window solid");
      }
      const double halfThickness = box->GetDZ() * Acts::UnitConstants::cm;
      const double radius        = std::max(box->GetDX(), box->GetDY()) * Acts::UnitConstants::cm;
      const double holeRadius    = aperture->GetRmax() * Acts::UnitConstants::cm;
      const double z             = physical.translation().z();
      const double upstream      = z - halfThickness - m_padding;
      const double downstream    = z + halfThickness + m_padding;
      const bool fitsGap         = zMin < upstream && downstream < zMax && rMin < rMax;
      if (!fitsGap) {
        throw std::runtime_error("Physical window slice is outside the native gap");
      }
      // Preserve the original tracking domain. The physical disc outside its
      // outer radius is outside this prototype's acceptance, as in the baseline.
      auto carrier = Acts::Transform3(Acts::Translation3(0., 0., z));
      auto windowLayer =
          Acts::DiscLayer::create(carrier, std::make_shared<Acts::RadialBounds>(rMin, rMax),
                                  nullptr, 2 * halfThickness, nullptr, Acts::passive);
      Acts::BinUtility bins(96, -std::numbers::pi, std::numbers::pi, Acts::closed,
                            Acts::AxisDirection::AxisPhi, physical);
      bins += Acts::BinUtility(64, holeRadius, radius, Acts::open, Acts::AxisDirection::AxisR);
      windowLayer->surfaceRepresentation().assignSurfaceMaterial(
          std::make_shared<Acts::ProtoSurfaceMaterial>(bins));
      Acts::TrackingVolumeVector parts;
      parts.push_back(Acts::CylinderVolumeHelper::createTrackingVolume(
          context, {}, {}, nullptr, rMin, rMax, zMin, upstream, name + "::BeforeWindow"));
      parts.push_back(Acts::CylinderVolumeHelper::createTrackingVolume(
          context, {windowLayer}, {}, nullptr, rMin, rMax, upstream, downstream,
          "B0WindowMaterial"));
      parts.push_back(Acts::CylinderVolumeHelper::createTrackingVolume(
          context, {}, {}, nullptr, rMin, rMax, downstream, zMax, name + "::AfterWindow"));
      for (const auto& part : parts) {
        if (!part) {
          throw std::runtime_error("Window gap subdivision failed");
        }
      }
      auto result = Acts::CylinderVolumeHelper::createContainerTrackingVolume(context, parts);
      if (!result) {
        throw std::runtime_error("Window gap gluing failed");
      }
      ++*m_inserted;
      ACTS_DEBUG("Inserted B0 window material volume at z="
                 << z << " mm in " << name << "; tracking bounds unchanged: r=[" << rMin << ", "
                 << rMax << "], z=[" << zMin << ", " << zMax << "] mm");
      return result;
    }

  private:
    dd4hep::DetElement m_window;
    double m_padding;
    std::shared_ptr<unsigned int> m_inserted;
    std::unique_ptr<const Acts::Logger> m_logger;
    const Acts::Logger& logger() const { return *m_logger; }
  };

} // namespace
std::unique_ptr<const Acts::TrackingGeometry> convertDD4hepDetector(
    dd4hep::DetElement world, const Acts::Logger& logger, Acts::BinningType phi,
    Acts::BinningType r, Acts::BinningType z, double envelopeR, double envelopeZ, double thickness,
    const std::function<void(std::vector<dd4hep::DetElement>&)>& sorter,
    const Acts::GeometryContext& context, std::shared_ptr<const Acts::IMaterialDecorator> material,
    std::shared_ptr<const Acts::GeometryIdentifierHook> hook,
    const ActsPlugins::DD4hepLayerBuilder::ElementFactory& factory) {
  auto window = findWindow(world);
  if (!window.isValid()) {
    throw std::runtime_error("B0 window material requires the DD4hep B0Window detector");
  }
  std::vector<dd4hep::DetElement> detectors;
  ActsPlugins::collectSubDetectors_dd4hep(world, detectors, logger);
  sorter(detectors);
  std::vector<std::shared_ptr<const Acts::CylinderVolumeBuilder>> builders;
  std::shared_ptr<const Acts::CylinderVolumeBuilder> beamPipe;
  bool added    = false;
  auto inserted = std::make_shared<unsigned int>(0);
  for (auto detector : detectors) {
    auto builder = ActsPlugins::volumeBuilder_dd4hep(detector, logger, phi, r, z, envelopeR,
                                                     envelopeZ, thickness, factory);
    if (!builder) {
      continue;
    }
    const bool isBeamPipe = builder->getConfiguration().buildToRadiusZero;
    if (detector.name() == std::string("B0TrackerSubAssembly")) {
      added = true;
    }
    if (isBeamPipe) {
      auto config                 = builder->getConfiguration();
      auto nativeHelper           = ActsPlugins::cylinderVolumeHelper_dd4hep(logger);
      config.trackingVolumeHelper = std::make_shared<WindowVolumeHelper>(
          nativeHelper->getConfiguration(), logger.clone("WindowVolumeHelper"), window, envelopeZ,
          inserted);
      builder = std::make_shared<Acts::CylinderVolumeBuilder>(config,
                                                              logger.clone("WindowVolumeBuilder"));
    }
    if (isBeamPipe) {
      if (beamPipe) {
        throw std::runtime_error("Multiple beam pipes");
      }
      beamPipe = builder;
    } else {
      builders.push_back(builder);
    }
  }
  if (!added) {
    throw std::runtime_error("B0 tracker assembly not found");
  }
  if (beamPipe) {
    builders.push_back(beamPipe);
  }
  Acts::TrackingGeometryBuilder::Config config;
  config.trackingVolumeHelper = ActsPlugins::cylinderVolumeHelper_dd4hep(logger);
  config.materialDecorator =
      material ? std::make_shared<WindowDecorator>(std::move(material), window) : nullptr;
  config.geometryIdentifierHook = std::move(hook);
  for (const auto& builder : builders) {
    config.trackingVolumeBuilders.emplace_back(
        [builder](const Acts::GeometryContext& c, const Acts::TrackingVolumePtr& inner,
                  const std::shared_ptr<const Acts::VolumeBounds>&) {
          return builder->trackingVolume(c, inner);
        });
  }
  auto geometry = Acts::TrackingGeometryBuilder(config).trackingGeometry(context);
  if (*inserted != 1) {
    throw std::runtime_error("Expected exactly one native window gap replacement");
  }
  // ACTS 44 copies a disc's material when creating adjacent navigation layers.
  // These are navigation aids, not additional physical material surfaces.
  geometry->visitVolumes([](const Acts::TrackingVolume* volume) {
    if (!volume->confinedLayers()) {
      return;
    }
    for (const auto& layer : volume->confinedLayers()->arrayObjects()) {
      if (layer->layerType() == Acts::navigation) {
        // Geometry construction owns these generated surfaces; ACTS exposes them as const.
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
        const_cast<Acts::Surface&>(layer->surfaceRepresentation()).assignSurfaceMaterial(nullptr);
      }
    }
  });
  return geometry;
}
} // namespace b0window
