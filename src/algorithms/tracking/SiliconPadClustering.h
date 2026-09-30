// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <DDRec/CellIDPositionConverter.h>
#include <DDSegmentation/CartesianGridXY.h>
#include <algorithms/algorithm.h>
#include <edm4eic/TrackerHitCollection.h>
#include <string>
#include <string_view>

#include "SiliconPadClusteringConfig.h"
#include "algorithms/interfaces/WithPodConfig.h"

namespace eicrecon {

using SiliconPadClusteringAlgorithm =
    algorithms::Algorithm<algorithms::Input<edm4eic::TrackerHitCollection>,
                          algorithms::Output<edm4eic::TrackerHitCollection>>;

/**
 * Merges touching pad hits into one cluster hit per particle crossing. Along each axis the
 * position is taken from the two leading columns of the cluster by inverting the trapezoid pad
 * response, so it is placed between the pads rather than at the charge centroid.
 */
class SiliconPadClustering : public SiliconPadClusteringAlgorithm,
                             public WithPodConfig<SiliconPadClusteringConfig> {

public:
  SiliconPadClustering(std::string_view name)
      : SiliconPadClusteringAlgorithm{name,
                                      {"inputPadHits"},
                                      {"outputClusterHits"},
                                      "Merge touching pad hits into cluster hits with a sub-pad "
                                      "position from charge sharing."} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  const dd4hep::rec::CellIDPositionConverter* m_converter = nullptr;
  const dd4hep::DDSegmentation::CartesianGridXY* m_grid   = nullptr;
};

} // namespace eicrecon
