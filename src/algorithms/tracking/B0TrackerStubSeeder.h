// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <memory>
#include <string>
#include <string_view>

#include "B0TrackerStubSeederConfig.h"
#include "algorithms/interfaces/WithPodConfig.h"

namespace dd4hep {
namespace DDSegmentation {
  class BitFieldCoder;
}
namespace rec {
  class CellIDPositionConverter;
}
} // namespace dd4hep

class ActsGeometryProvider;

namespace eicrecon {

using B0TrackerStubSeederAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::Measurement2DCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

/// Dedicated seeder for the B0 tracker (dipole spectrometer at z ~ 6 m).
///
/// Fits one hit per station with a straight line (non-bend plane) and a
/// parabola (bend plane) in the ion-rotated frame, extracts direction,
/// momentum (from the dipole sagitta) and charge, back-extrapolates
/// analytically to the origin, and emits seed parameters on the origin
/// perigee surface that CKFTracking expects.
class B0TrackerStubSeeder : public B0TrackerStubSeederAlgorithm,
                            public WithPodConfig<B0TrackerStubSeederConfig> {
public:
  B0TrackerStubSeeder(std::string_view name)
      : B0TrackerStubSeederAlgorithm{name,
                                     {"inputMeasurements"},
                                     {"outputSeeds", "outputTrackParameters"},
                                     "stub-based track seeds for the B0 dipole tracker"} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  const dd4hep::rec::CellIDPositionConverter* m_converter{nullptr};
  const dd4hep::DDSegmentation::BitFieldCoder* m_decoder{nullptr};
  std::shared_ptr<const ActsGeometryProvider> m_acts_context;
};

} // namespace eicrecon
