// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "ActsToTracks.h"
#include "ActsGeometryProvider.h"
#include <memory>
#include <string_view>
namespace eicrecon {
/// Separate exporter for local-plane parameters. It never silently describes
/// a B0-plane position/momentum as a vertex/perigee measurement.
class B0TelescopeTrackExport : public ActsToTracksAlgorithm, public WithPodConfig<NoConfig> {
public:
  explicit B0TelescopeTrackExport(std::string_view name)
      : ActsToTracksAlgorithm{name, {"measurements", "seeds", "states", "tracks", "rawAssociations"},
                                {"trajectories", "parameters", "outputTracks", "links", "associations"},
                                "Export local B0 tracks with their actual reference surface"} {}
  void init() final;
  void process(const Input&, const Output&) const final;
private:
  std::shared_ptr<const ActsGeometryProvider> m_provider;
};
} // namespace eicrecon
