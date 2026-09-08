// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <Acts/Geometry/ApproachDescriptor.hpp>
#include <Acts/Material/ProtoSurfaceMaterial.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <fmt/format.h>
#include <stdexcept>

namespace eicrecon {

/// B0 maps material onto both disc faces, not the radial cylinder approaches.
/// A leftover prototype is a mapping instruction, not usable fit material.
inline void validateB0ApproachMaterial(const Acts::ApproachDescriptor& approaches) {
  unsigned discs = 0;
  for (const auto* surface : approaches.containedSurfaces()) {
    if (dynamic_cast<const Acts::DiscSurface*>(surface) == nullptr) {
      continue;
    }
    ++discs;
    const auto* material = surface->surfaceMaterial();
    bool prototype       = dynamic_cast<const Acts::ProtoSurfaceMaterial*>(material) != nullptr;
#if Acts_VERSION_MAJOR >= 47
    prototype |= dynamic_cast<const Acts::ProtoGridSurfaceMaterial*>(material) != nullptr;
#endif
    if (material == nullptr || prototype) {
      const auto id = surface->geometryId();
      throw std::runtime_error(fmt::format(
          "B0Tracker requires mapped material on approach surface "
          "(volume={}, layer={}, approach={}): {}. Check acts:MaterialMap against the geometry.",
          id.volume(), id.layer(), id.approach(),
          prototype ? "unmapped prototype material" : "missing material"));
    }
  }
  if (discs != 2) {
    throw std::runtime_error(fmt::format(
        "B0Tracker requires two planar material approaches per layer; found {}", discs));
  }
}

} // namespace eicrecon
