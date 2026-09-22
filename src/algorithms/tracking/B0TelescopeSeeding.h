// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#pragma once

#include "B0TelescopeActs.h"
#if EICRECON_HAS_B0_TELESCOPE
#include <algorithms/algorithm.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <memory>
#include <string_view>
#include "ActsGeometryProvider.h"
#include "B0TelescopeSeedingConfig.h"
#include "algorithms/interfaces/WithPodConfig.h"

namespace dd4hep::DDSegmentation {
class BitFieldCoder;
}

namespace eicrecon {
using B0TelescopeSeedingAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackerHitCollection>,
    algorithms::Output<edm4eic::TrackSeedCollection, edm4eic::TrackParametersCollection>>;

class B0TelescopeSeeding : public B0TelescopeSeedingAlgorithm,
                          public WithPodConfig<B0TelescopeSeedingConfig> {
public:
  explicit B0TelescopeSeeding(std::string_view name)
      : B0TelescopeSeedingAlgorithm{name, {"inputHits"}, {"outputSeeds", "outputParameters"},
                                    "B0 telescope seeding with ACTS station groups"} {}
  void init() final;
  void process(const Input&, const Output&) const final;

private:
  std::shared_ptr<const ActsGeometryProvider> m_geo;
  const dd4hep::DDSegmentation::BitFieldCoder* m_decoder{};
  long m_system{};
  double m_rotation{};
};
} // namespace eicrecon
#endif
