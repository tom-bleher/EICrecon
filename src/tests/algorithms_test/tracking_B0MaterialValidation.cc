// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Geometry/GenericApproachDescriptor.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>

#include "algorithms/tracking/B0MaterialValidation.h"

TEST_CASE("B0 material targets use identifiers assigned after decoration",
          "[tracking][b0material]") {
  const auto surface =
      Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), 0., 100.);
  // A NavigationLayer owns a separate representing surface whose ID is still
  // zero when the material decorator sees its copied prototype.
  surface->assignSurfaceMaterial(std::make_shared<Acts::ProtoSurfaceMaterial>());
  const std::vector<const Acts::Surface*> declared{surface.get()};
  const auto finalId = Acts::GeometryIdentifier{}.withVolume(7).withLayer(1);
  surface->assignGeometryId(finalId);
  const auto targets = eicrecon::b0MaterialTargetsAfterClosure(declared);
  REQUIRE(targets.count(finalId) == 1);
  CHECK_THROWS_AS(eicrecon::validateB0MaterialTargets(targets, {7}), std::runtime_error);
}

TEST_CASE("B0 external passive declarations require matching ACTS geometry",
          "[tracking][b0material]") {
  const auto surface =
      Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), 5., 100.);
  const auto id = Acts::GeometryIdentifier{}.withVolume(7).withLayer(12);
  surface->assignGeometryId(id);
  const std::map<Acts::GeometryIdentifier, const Acts::Surface*> expected{{id, surface.get()}};
#if Acts_VERSION_MAJOR >= 45
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
#else
  const auto gctx = Acts::GeometryContext{};
#endif
  auto transform = Acts::Transform3::Identity();
  CHECK_NOTHROW(eicrecon::validateB0PassiveDiscRepresentation(expected, {7}, gctx, transform, 5.,
                                                              100., "window"));
  CHECK_THROWS_AS(eicrecon::validateB0PassiveDiscRepresentation({}, {7}, gctx, transform, 5., 100.,
                                                                "ignored reference"),
                  std::runtime_error);
  CHECK_THROWS_AS(eicrecon::validateB0PassiveDiscRepresentation(expected, {7}, gctx, transform, 5.,
                                                                101., "wrong bounds"),
                  std::runtime_error);
  transform.translation().z() = 1.;
  CHECK_THROWS_AS(eicrecon::validateB0PassiveDiscRepresentation(expected, {7}, gctx, transform, 5.,
                                                                100., "shifted window"),
                  std::runtime_error);
  transform.translation().z() = 0.;
  transform.linear()(0, 0)    = -1.;
  transform.linear()(2, 2)    = -1.;
  CHECK_THROWS_AS(eicrecon::validateB0PassiveDiscRepresentation(expected, {7}, gctx, transform, 5.,
                                                                100., "reversed normal"),
                  std::runtime_error);
}

TEST_CASE("B0 validates declared passive material without sensors", "[tracking][b0material]") {
  const auto window =
      Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), 0., 100.);
  const auto id = Acts::GeometryIdentifier{}.withVolume(7).withLayer(12);
  window->assignGeometryId(id);
  const std::map<Acts::GeometryIdentifier, const Acts::Surface*> expected{{id, window.get()}};

  SECTION("missing or prototype passive material is rejected") {
    CHECK_THROWS_AS(eicrecon::validateB0MaterialTargets(expected, {7}), std::runtime_error);
    window->assignSurfaceMaterial(std::make_shared<Acts::ProtoSurfaceMaterial>());
    CHECK_THROWS_AS(eicrecon::validateB0MaterialTargets(expected, {7}), std::runtime_error);
  }
  SECTION("a mapped passive target is accepted") {
    window->assignSurfaceMaterial(
        std::make_shared<Acts::HomogeneousSurfaceMaterial>(Acts::MaterialSlab(
            Acts::Material::fromMolarDensity(17.6, 168., 55.85, 26., 0.000141), 8.)));
    CHECK_NOTHROW(eicrecon::validateB0MaterialTargets(expected, {7}));
  }
  SECTION("mapped sensor approaches do not satisfy a missing passive target") {
    const auto approach =
        Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), 0., 100.);
    const auto approachId = id.withLayer(14).withApproach(1);
    approach->assignGeometryId(approachId);
    approach->assignSurfaceMaterial(
        std::make_shared<Acts::HomogeneousSurfaceMaterial>(Acts::MaterialSlab(
            Acts::Material::fromMolarDensity(93.7, 465.2, 28.0855, 14., 0.000083), 0.3)));
    auto both = expected;
    both.emplace(approachId, approach.get());
    CHECK_THROWS_AS(eicrecon::validateB0MaterialTargets(both, {7}), std::runtime_error);
  }
  SECTION("unrelated detector volumes do not become B0 requirements") {
    CHECK_NOTHROW(eicrecon::validateB0MaterialTargets(expected, {4}));
  }
  SECTION("legacy geometries have no undeclared passive requirement") {
    CHECK_NOTHROW(eicrecon::validateB0MaterialTargets({}, {7}));
  }
}

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
