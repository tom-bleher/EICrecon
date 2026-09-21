// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <Acts/Geometry/ApproachDescriptor.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Material/BinnedSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/ProtoSurfaceMaterial.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <fmt/format.h>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace eicrecon {

inline bool isPrototypeMaterial(const Acts::ISurfaceMaterial* material) {
  bool prototype = dynamic_cast<const Acts::ProtoSurfaceMaterial*>(material) != nullptr;
#if Acts_VERSION_MAJOR >= 47
  prototype |= dynamic_cast<const Acts::ProtoGridSurfaceMaterial*>(material) != nullptr;
#endif
  return prototype;
}

/// JSON maps use homogeneous or binned slabs. Empty bins are legitimate, but
/// a complete B0 material approach must contain some finite, physical material.
inline bool hasUsableB0Material(const Acts::ISurfaceMaterial* material) {
  bool hasMaterial = false;
  const auto valid = [&hasMaterial](const Acts::MaterialSlab& slab) {
    if (!std::isfinite(slab.thickness()) || slab.thickness() < 0) {
      return false;
    }
    // Vacuum has infinite radiation/interaction lengths by definition.
    if (slab.material() == Acts::Material::Vacuum()) {
      return true;
    }
    const auto& medium = slab.material();
    if (!medium.parameters().allFinite() || medium.Ar() <= 0 || medium.Z() <= 0 ||
        medium.X0() <= 0 || medium.L0() <= 0 || medium.molarDensity() <= 0 ||
        !std::isfinite(slab.thicknessInX0()) || !std::isfinite(slab.thicknessInL0())) {
      return false;
    }
    hasMaterial |= slab.thickness() > 0;
    return true;
  };
  if (const auto* binned = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(material)) {
    for (const auto& row : binned->fullMaterial()) {
      for (const auto& slab : row) {
        if (!valid(slab)) {
          return false;
        }
      }
    }
    return hasMaterial;
  }
  if (const auto* homogeneous = dynamic_cast<const Acts::HomogeneousSurfaceMaterial*>(material)) {
    return valid(homogeneous->materialSlab(Acts::Vector2::Zero().eval())) && hasMaterial;
  }
  return false;
}

/// Resolve identifiers after closure: some representing surfaces are decorated
/// before their final volume/layer identifier is assigned.
inline std::map<Acts::GeometryIdentifier, const Acts::Surface*>
b0MaterialTargetsAfterClosure(const std::vector<const Acts::Surface*>& surfaces) {
  std::map<Acts::GeometryIdentifier, const Acts::Surface*> targets;
  for (const auto* surface : surfaces) {
    targets.insert_or_assign(surface->geometryId(), surface);
  }
  return targets;
}

/// Validate declared targets, including passive layers without sensors.
/// Legacy geometries without a declaration do not acquire a new requirement.
inline void validateB0MaterialTargets(
    const std::map<Acts::GeometryIdentifier, const Acts::Surface*>& expectedSurfaces,
    const std::set<Acts::GeometryIdentifier::Value>& b0Volumes) {
  for (const auto& [id, surface] : expectedSurfaces) {
    if (b0Volumes.count(id.volume()) == 0) {
      continue;
    }
    const auto* material = surface == nullptr ? nullptr : surface->surfaceMaterial();
    if (!hasUsableB0Material(material)) {
      throw std::runtime_error(
          fmt::format("B0Tracker requires mapped material on geometry-declared surface "
                      "(volume={}, boundary={}, layer={}, approach={}, sensitive={}): {}. "
                      "Check acts:MaterialMap against the geometry.",
                      id.volume(), id.boundary(), id.layer(), id.approach(), id.sensitive(),
                      material == nullptr             ? "missing material"
                      : isPrototypeMaterial(material) ? "unmapped prototype material"
                                                      : "unusable material payload"));
    }
  }
}

/// An older DD4hep converter can silently ignore external passive-layer
/// metadata. Require the declared disc's actual pose and bounds, not just a
/// material count or converter version string.
inline void validateB0PassiveDiscRepresentation(
    const std::map<Acts::GeometryIdentifier, const Acts::Surface*>& expectedSurfaces,
    const std::set<Acts::GeometryIdentifier::Value>& b0Volumes, const Acts::GeometryContext& gctx,
    const Acts::Transform3& transform, double rMin, double rMax, std::string_view path) {
  for (const auto& [id, surface] : expectedSurfaces) {
    if (surface == nullptr || b0Volumes.count(id.volume()) == 0 || id.approach() != 0 ||
        id.sensitive() != 0) {
      continue;
    }
    const auto* bounds = dynamic_cast<const Acts::RadialBounds*>(&surface->bounds());
    if (bounds == nullptr) {
      continue;
    }
#if Acts_VERSION_MAJOR >= 47
    const auto& actual = surface->localToGlobalTransform(gctx);
#else
    const auto& actual = surface->transform(gctx);
#endif
    if ((actual.translation() - transform.translation()).norm() < 1e-6 &&
        (actual.linear() - transform.linear()).norm() < 1e-10 &&
        std::abs(bounds->get(Acts::RadialBounds::eMinR) - rMin) < 1e-6 &&
        std::abs(bounds->get(Acts::RadialBounds::eMaxR) - rMax) < 1e-6) {
      return;
    }
  }
  throw std::runtime_error(fmt::format(
      "B0Tracker external passive layer '{}' has no geometry-declared ACTS material disc "
      "at its DD4hep placement and bounds. Load a DD4hep converter supporting "
      "passive_layer_references and passive_disc.",
      path));
}

/// Expectations come from the geometry's prototype decoration, before the JSON
/// decorator replaces it. Do not assume every geometry maps both layer faces.
inline void validateB0ApproachMaterial(const Acts::ApproachDescriptor& approaches,
                                       const std::set<Acts::GeometryIdentifier>& expectedMaterial) {
  unsigned expected = 0;
  for (const auto* surface : approaches.containedSurfaces()) {
    const auto id = surface->geometryId();
    if (expectedMaterial.count(id) == 0) {
      continue;
    }
    ++expected;
    const auto* material = surface->surfaceMaterial();
    if (!hasUsableB0Material(material)) {
      throw std::runtime_error(fmt::format(
          "B0Tracker requires mapped material on approach surface "
          "(volume={}, layer={}, approach={}): {}. Check acts:MaterialMap against the geometry.",
          id.volume(), id.layer(), id.approach(),
          material == nullptr             ? "missing material"
          : isPrototypeMaterial(material) ? "unmapped prototype material"
                                          : "unusable material payload"));
    }
  }
  if (expected == 0) {
    throw std::runtime_error("B0Tracker layer has no geometry-declared material approaches");
  }
}

} // namespace eicrecon
