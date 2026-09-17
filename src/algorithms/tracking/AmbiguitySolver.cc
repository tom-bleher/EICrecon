// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Minjung Kim, Barak Schmookler
#include "AmbiguitySolver.h"
#include "B0ReconstructionCounters.h"

#include <Acts/AmbiguityResolution/GreedyAmbiguityResolution.hpp>
#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/EventData/ProxyAccessor.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackStatePropMask.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <ActsExamples/EventData/IndexSourceLink.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <boost/container/flat_set.hpp>
#include <spdlog/common.h>
#include <Eigen/Core>
#include <Eigen/LU> // IWYU pragma: keep
#include <any>
#include <cstddef>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "Acts/Utilities/Logger.hpp"
#include "AmbiguitySolverConfig.h"
#include "extensions/spdlog/SpdlogFormatters.h" // IWYU pragma: keep
#include "extensions/spdlog/SpdlogToActs.h"

namespace eicrecon {

Acts::GreedyAmbiguityResolution::Config
transformConfig(const eicrecon::AmbiguitySolverConfig& cfg) {
  Acts::GreedyAmbiguityResolution::Config result;
  result.maximumSharedHits = cfg.maximum_shared_hits;
  result.maximumIterations = cfg.maximum_iterations;
  result.nMeasurementsMin  = cfg.n_measurements_min;
  return result;
}

static std::size_t sourceLinkHash(const Acts::SourceLink& a) {
  return static_cast<std::size_t>(a.get<ActsExamples::IndexSourceLink>().index());
}

static bool sourceLinkEquality(const Acts::SourceLink& a, const Acts::SourceLink& b) {
  return a.get<ActsExamples::IndexSourceLink>().index() ==
         b.get<ActsExamples::IndexSourceLink>().index();
}

void AmbiguitySolver::init() {
  // Convert algorithm log level to Acts log level
  const auto spdlog_level = static_cast<spdlog::level::level_enum>(this->level());
  const auto acts_level   = eicrecon::SpdlogToActsLevel(spdlog_level);

  // Create Acts logger with appropriate level
  m_acts_logger = Acts::getDefaultLogger("AmbiguitySolver", acts_level);
  m_acts_cfg    = transformConfig(m_cfg);
  m_core = std::make_unique<Acts::GreedyAmbiguityResolution>(m_acts_cfg, acts_logger().clone());
}

void AmbiguitySolver::process(const Input& input, const Output& output) const {
  const auto [input_track_states, input_tracks] = input;
  auto [output_track_states, output_tracks]     = output;

  // Construct ConstTrackContainer from underlying containers
  auto trackStateContainer =
      std::make_shared<Acts::ConstVectorMultiTrajectory>(*input_track_states);
  auto trackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(*input_tracks);
  ActsExamples::ConstTrackContainer input_trks(trackContainer, trackStateContainer);

  // B0 CKF failure diagnostics (see b0counters::ckfdiag): the unfiltered input
  // may carry CKF-rejected candidates and stateless find-failure markers.
  // They must not enter ambiguity resolution; only accepted candidates are
  // resolved, so the output matches a CKF that drops rejects. Containers
  // without the column (central tracking, older files, unit tests) resolve
  // every row exactly as before.
  Acts::ConstProxyAccessor<unsigned int> ckfStatus(b0counters::ckfdiag::kStatusColumn);
  // hasColumn also covers empty containers, where probing a first track is impossible.
  const bool haveCkfStatus = input_trks.hasColumn(b0counters::ckfdiag::kStatusColumn);

  Acts::GreedyAmbiguityResolution::State state;
  ActsExamples::TrackContainer acceptedTracks{std::make_shared<Acts::VectorTrackContainer>(),
                                              std::make_shared<Acts::VectorMultiTrajectory>()};
  if (haveCkfStatus) {
    acceptedTracks.ensureDynamicColumns(input_trks);
    for (auto track : input_trks) {
      unsigned status = b0counters::ckfdiag::kAccepted;
      try {
        status = ckfStatus(track);
      } catch (...) {
      }
      if (!b0counters::ckfdiag::resolveRow(true, status)) {
        continue;
      }
      auto destProxy = acceptedTracks.makeTrack();
      destProxy.copyFrom(track);
    }
  }

  ActsExamples::TrackContainer solvedTracks{std::make_shared<Acts::VectorTrackContainer>(),
                                            std::make_shared<Acts::VectorMultiTrajectory>()};
  std::size_t nResolvedIn = 0;
  if (!haveCkfStatus) {
    m_core->computeInitialState(input_trks, state, &sourceLinkHash, &sourceLinkEquality);
    m_core->resolve(state);
    solvedTracks.ensureDynamicColumns(input_trks);
    for (auto iTrack : state.selectedTracks) {
      auto destProxy = solvedTracks.makeTrack();
      auto srcProxy  = input_trks.getTrack(state.trackTips.at(iTrack));
      destProxy.copyFrom(srcProxy);
    }
    nResolvedIn = input_trks.size();
  } else {
    m_core->computeInitialState(acceptedTracks, state, &sourceLinkHash, &sourceLinkEquality);
    m_core->resolve(state);
    solvedTracks.ensureDynamicColumns(acceptedTracks);
    for (auto iTrack : state.selectedTracks) {
      auto destProxy = solvedTracks.makeTrack();
      auto srcProxy  = acceptedTracks.getTrack(state.trackTips.at(iTrack));
      destProxy.copyFrom(srcProxy);
    }
    nResolvedIn = acceptedTracks.size();
  }

  if (b0counters::isB0AlgorithmName(this->name())) {
    const auto chain = b0counters::chainFromAlgorithmName(this->name());
    const auto nIn   = nResolvedIn;
    const auto nOut  = solvedTracks.size();
    b0counters::incrementAmbiguity(chain, b0counters::AmbiguityStat::tracksIn, nIn);
    b0counters::incrementAmbiguity(chain, b0counters::AmbiguityStat::tracksOut, nOut);
    if (nIn > nOut) {
      b0counters::incrementAmbiguity(chain, b0counters::AmbiguityStat::ambiguityRejected,
                                     nIn - nOut);
    }
  }

  // Allocate new const containers and assign pointers to outputs
  *output_track_states =
      new Acts::ConstVectorMultiTrajectory(std::move(solvedTracks.trackStateContainer()));
  *output_tracks = new Acts::ConstVectorTrackContainer(std::move(solvedTracks.container()));
}

} // namespace eicrecon
