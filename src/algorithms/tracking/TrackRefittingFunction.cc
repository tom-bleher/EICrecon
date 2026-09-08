// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include "CKFTracking.h"

#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#include <utility>

namespace eicrecon {
namespace {
  using Propagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
  using Fitter     = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;

  struct TrackFitterFunctionImpl : CKFTracking::TrackFitterFunction {
    Fitter fitter;
    explicit TrackFitterFunctionImpl(Fitter&& value) : fitter(std::move(value)) {}
    CKFTracking::TrackFitterResult operator()(const std::vector<Acts::SourceLink>& measurements,
                                              const ActsExamples::TrackParameters& initial,
                                              const CKFTracking::TrackFitterOptions& options,
                                              ActsExamples::TrackContainer& tracks) const override {
      return fitter.fit(measurements.begin(), measurements.end(), initial, options, tracks);
    }
  };
} // namespace

std::shared_ptr<CKFTracking::TrackFitterFunction>
CKFTracking::makeTrackFitterFunction(std::shared_ptr<const Acts::TrackingGeometry> geometry,
                                     std::shared_ptr<const Acts::MagneticFieldProvider> field,
                                     const Acts::Logger& logger) {
  Acts::Navigator::Config config{.trackingGeometry = std::move(geometry)};
  config.resolvePassive   = false;
  config.resolveMaterial  = true;
  config.resolveSensitive = true;
  Propagator propagator(Acts::EigenStepper<>(std::move(field)),
                        Acts::Navigator(config, logger.cloneWithSuffix("RefitNavigator")));
  return std::make_shared<TrackFitterFunctionImpl>(
      Fitter(std::move(propagator), logger.cloneWithSuffix("Refit")));
}
} // namespace eicrecon
