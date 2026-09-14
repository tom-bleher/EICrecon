// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2016 CERN for the benefit of the ACTS project.

#pragma once
#include <ActsPlugins/DD4hep/ConvertDD4hepDetector.hpp>
namespace b0window {
bool isWindowSurface(const Acts::Surface& surface, dd4hep::DetElement window,
                     const Acts::GeometryContext& context);
std::unique_ptr<const Acts::TrackingGeometry> convertDD4hepDetector(
    dd4hep::DetElement world, const Acts::Logger& logger, Acts::BinningType phi = Acts::equidistant,
    Acts::BinningType r = Acts::equidistant, Acts::BinningType z = Acts::equidistant,
    double envelopeR = Acts::UnitConstants::mm, double envelopeZ = 5 * Acts::UnitConstants::mm,
    double thickness = Acts::UnitConstants::fm,
    const std::function<void(std::vector<dd4hep::DetElement>&)>& sorter =
        ActsPlugins::sortDetElementsByID,
    const Acts::GeometryContext& context                     = Acts::GeometryContext(),
    std::shared_ptr<const Acts::IMaterialDecorator> material = nullptr,
    std::shared_ptr<const Acts::GeometryIdentifierHook> hook =
        std::make_shared<Acts::GeometryIdentifierHook>(),
    const ActsPlugins::DD4hepLayerBuilder::ElementFactory& factory =
        ActsPlugins::DD4hepLayerBuilder::defaultDetectorElementFactory);
} // namespace b0window
