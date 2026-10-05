// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck

#pragma once

#include <Acts/Definitions/Direction.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Result.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>
#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/TrackParameters.hpp>
#endif
#include <ActsExamples/EventData/Track.hpp>
#include <algorithms/algorithm.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <Eigen/Core>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "CKFTrackingConfig.h"
#include "algorithms/interfaces/ActsSvc.h"
#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"

namespace eicrecon {

using CKFTrackingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackSeedCollection, edm4eic::Measurement2DCollection>,
    algorithms::Output<Acts::ConstVectorMultiTrajectory*, Acts::ConstVectorTrackContainer*>>;

/** Fitting algorithm implementation .
 *
 * \ingroup tracking
 */

class CKFTracking : public CKFTrackingAlgorithm, public WithPodConfig<eicrecon::CKFTrackingConfig> {
public:
  /// Track finder function that takes input measurements, initial trackstate
  /// and track finder options and returns some track-finder-specific result.
  using TrackFinderOptions = Acts::CombinatorialKalmanFilterOptions<ActsExamples::TrackContainer>;
  using TrackFinderResult  = Acts::Result<std::vector<ActsExamples::TrackContainer::TrackProxy>>;

  /// Find function that takes the above parameters
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class CKFTrackingFunction {
  public:
    virtual ~CKFTrackingFunction() = default;

    virtual TrackFinderResult operator()(const ActsExamples::TrackParameters&,
                                         const TrackFinderOptions&,
                                         ActsExamples::TrackContainer&) const = 0;
  };

  /// Create the track finder function implementation.
  /// The magnetic field is intentionally given by-value since the variantresults
  /// contains shared_ptr anyways.
  static std::shared_ptr<CKFTrackingFunction>
  makeCKFTrackingFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                          std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
                          const Acts::Logger& logger);

  CKFTracking(std::string_view name)
      : CKFTrackingAlgorithm{name,
                             {"inputTrackParameters", "inputMeasurements"},
                             {"outputActsTrackStates", "outputActsTracks"},
                             "Combinatorial Kalman Filter track finding"} {}

  /// Count longitudinal groups using only accepted measurement surfaces.
  static std::size_t countMeasurementGroups(const ActsExamples::TrackContainer::TrackProxy& track,
                                            const Acts::GeometryContext& gctx, double zGap);

  /// Extrapolate a smoothed track upstream to a reference surface.
  ///
  /// Starts from the first measurement state and propagates explicitly
  /// backward. Unlike Acts' firstOrLast strategy, the direction does not come
  /// from a straight-line intersection, which misleads for dipole-bent
  /// forward tracks. Uses ordinary SurfaceReached: ForcedSurfaceReached
  /// accepts arbitrarily negative intersections and can reverse the stepping
  /// direction towards a downstream perigee instead. Omit EndOfWorldReached
  /// from the options when the reference can lie outside the tracking
  /// geometry; material beyond the modeled geometry is then absent by
  /// construction.
  template <typename propagator_t, typename propagator_options_t>
  static Acts::Result<void> extrapolateBackwardToReference(
      ActsExamples::TrackContainer::TrackProxy& track, const Acts::Surface& referenceSurface,
      const propagator_t& propagator, propagator_options_t options, const Acts::Logger& logger) {
    auto firstMeasurement = Acts::findFirstMeasurementState(track);
    if (!firstMeasurement.ok()) {
      return firstMeasurement.error();
    }
    const auto parameters = track.createParametersFromState(*firstMeasurement);
    options.direction     = Acts::Direction::Backward();
#if Acts_VERSION_MAJOR >= 46
    auto result = propagator.template propagate<propagator_options_t, Acts::SurfaceReached>(
        parameters, referenceSurface, options);
#else
    auto result =
        propagator.template propagate<Acts::BoundTrackParameters, propagator_options_t,
                                      Acts::SurfaceReached>(parameters, referenceSurface, options);
#endif
    if (!result.ok()) {
      return result.error();
    }
    track.setReferenceSurface(referenceSurface.getSharedPtr());
    track.parameters() = result->endParameters.value().parameters();
    track.covariance() = result->endParameters.value().covariance().value();
    return Acts::Result<void>::success();
  }

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  std::shared_ptr<const Acts::Logger> m_acts_logger{nullptr};
  std::shared_ptr<CKFTrackingFunction> m_trackFinderFunc;
  std::shared_ptr<const ActsGeometryProvider> m_geoSvc{
      algorithms::ActsSvc::instance().acts_geometry_provider()};
  std::shared_ptr<const Acts::MagneticFieldProvider> m_BField{m_geoSvc->getFieldProvider()};

  Acts::MeasurementSelector::Config m_sourcelinkSelectorCfg;

  /// Private access to the logging instance
  const Acts::Logger& acts_logger() const { return *m_acts_logger; }
};

} // namespace eicrecon
