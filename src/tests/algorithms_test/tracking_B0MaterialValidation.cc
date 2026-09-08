// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Geometry/GenericApproachDescriptor.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>

#include "algorithms/tracking/B0MaterialValidation.h"

TEST_CASE("B0 validates every geometry-declared material approach", "[tracking][b0material]") {
  const auto slab = Acts::MaterialSlab(
      Acts::Material::fromMolarDensity(93.7, 465.2, 28.0855, 14., 0.000083), 0.3);
  const auto mapped = std::make_shared<Acts::HomogeneousSurfaceMaterial>(slab);
  const auto first =
      Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), 0., 100.);
  const auto second =
      Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), 0., 100.);
  first->assignGeometryId(Acts::GeometryIdentifier{}.withVolume(2).withLayer(2).withApproach(1));
  second->assignGeometryId(Acts::GeometryIdentifier{}.withVolume(2).withLayer(2).withApproach(2));
  const Acts::GenericApproachDescriptor approaches({first, second});
  const std::set<Acts::GeometryIdentifier> expected{first->geometryId(), second->geometryId()};
  first->assignSurfaceMaterial(mapped);
  second->assignSurfaceMaterial(mapped);
  CHECK_NOTHROW(eicrecon::validateB0ApproachMaterial(approaches, expected));

  SECTION("one missing face cannot be hidden by the other mapped face") {
    second->assignSurfaceMaterial(nullptr);
    CHECK_THROWS_AS(eicrecon::validateB0ApproachMaterial(approaches, expected), std::runtime_error);
  }
  SECTION("leftover mapping prototypes are rejected") {
    second->assignSurfaceMaterial(std::make_shared<Acts::ProtoSurfaceMaterial>());
    CHECK_THROWS_AS(eicrecon::validateB0ApproachMaterial(approaches, expected), std::runtime_error);
  }
  SECTION("a geometry mapping only one face is valid") {
    second->assignSurfaceMaterial(nullptr);
    CHECK_NOTHROW(eicrecon::validateB0ApproachMaterial(approaches, {first->geometryId()}));
  }
  SECTION("a B0 layer needs declared material targets") {
    CHECK_THROWS_AS(eicrecon::validateB0ApproachMaterial(approaches, {}), std::runtime_error);
  }
  SECTION("vacuum bins are allowed but wholly empty maps are rejected") {
    const Acts::BinUtility bins(2, 0., 100., Acts::open, Acts::AxisDirection::AxisR);
    second->assignSurfaceMaterial(std::make_shared<Acts::BinnedSurfaceMaterial>(
        bins, Acts::MaterialSlabVector{Acts::MaterialSlab{}, slab}));
    CHECK_NOTHROW(eicrecon::validateB0ApproachMaterial(approaches, expected));
    second->assignSurfaceMaterial(std::make_shared<Acts::BinnedSurfaceMaterial>(
        bins, Acts::MaterialSlabVector{Acts::MaterialSlab{}, Acts::MaterialSlab{}}));
    CHECK_THROWS_AS(eicrecon::validateB0ApproachMaterial(approaches, expected), std::runtime_error);
    second->assignSurfaceMaterial(
        std::make_shared<Acts::HomogeneousSurfaceMaterial>(Acts::MaterialSlab{}));
    CHECK_THROWS_AS(eicrecon::validateB0ApproachMaterial(approaches, expected), std::runtime_error);
  }
  SECTION("one physical bin cannot hide a corrupt bin") {
    const Acts::BinUtility bins(2, 0., 100., Acts::open, Acts::AxisDirection::AxisR);
    second->assignSurfaceMaterial(std::make_shared<Acts::BinnedSurfaceMaterial>(
        bins,
        Acts::MaterialSlabVector{
            slab, Acts::MaterialSlab(slab.material(), std::numeric_limits<float>::quiet_NaN())}));
    CHECK_THROWS_AS(eicrecon::validateB0ApproachMaterial(approaches, expected), std::runtime_error);
  }
  SECTION("nonfinite payloads are rejected") {
    second->assignSurfaceMaterial(std::make_shared<Acts::HomogeneousSurfaceMaterial>(
        Acts::MaterialSlab(slab.material(), std::numeric_limits<float>::quiet_NaN())));
    CHECK_THROWS_AS(eicrecon::validateB0ApproachMaterial(approaches, expected), std::runtime_error);
  }
}
