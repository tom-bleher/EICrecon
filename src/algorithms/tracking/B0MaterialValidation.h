// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <Acts/Geometry/ApproachDescriptor.hpp>
#include <Acts/Material/BinnedSurfaceMaterial.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/ProtoSurfaceMaterial.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <fmt/format.h>
#include <cmath>
#include <set>
#include <stdexcept>

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
