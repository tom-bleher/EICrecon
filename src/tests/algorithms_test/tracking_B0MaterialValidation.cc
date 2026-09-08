// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Geometry/GenericApproachDescriptor.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <memory>
#include <vector>

#include "algorithms/tracking/B0MaterialValidation.h"

TEST_CASE("B0 requires mapped material on both planar approaches", "[tracking][B0Material]") {
  const auto bounds = std::make_shared<Acts::RadialBounds>(1., 20.);
  auto inner = Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), bounds);
  auto outer = Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), bounds);
  inner->assignGeometryId(Acts::GeometryIdentifier{}.withVolume(5).withLayer(2).withApproach(1));
  outer->assignGeometryId(Acts::GeometryIdentifier{}.withVolume(5).withLayer(2).withApproach(2));
  // The radial boundary is deliberately not decorated by B0 layer_material XML.
  auto radial =
      Acts::Surface::makeShared<Acts::CylinderSurface>(Acts::Transform3::Identity(), 20., 1.);
  Acts::GenericApproachDescriptor approaches({inner, outer, radial});
  const auto mapped = std::make_shared<Acts::HomogeneousSurfaceMaterial>();

  SECTION("missing map entries fail") {
    CHECK_THROWS_WITH(eicrecon::validateB0ApproachMaterial(approaches),
                      Catch::Matchers::ContainsSubstring("missing material"));
  }
  SECTION("a partial map does not pass on one decorated face") {
    inner->assignSurfaceMaterial(mapped);
    CHECK_THROWS_WITH(eicrecon::validateB0ApproachMaterial(approaches),
                      Catch::Matchers::ContainsSubstring("approach=2"));
  }
  SECTION("leftover mapping prototypes fail despite non-null pointers") {
    inner->assignSurfaceMaterial(mapped);
    outer->assignSurfaceMaterial(std::make_shared<Acts::ProtoSurfaceMaterial>());
    CHECK_THROWS_WITH(eicrecon::validateB0ApproachMaterial(approaches),
                      Catch::Matchers::ContainsSubstring("unmapped prototype material"));
  }
  SECTION("mapped faces pass with an intentionally unmapped radial boundary") {
    inner->assignSurfaceMaterial(mapped);
    outer->assignSurfaceMaterial(mapped);
    CHECK_NOTHROW(eicrecon::validateB0ApproachMaterial(approaches));
  }
}
