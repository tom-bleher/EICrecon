// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/MCRecoTrackParticleAssociationCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrajectoryCollection.h>
#include <edm4eic/VertexCollection.h>
#include <cstdint>
#include <memory>
#include <string>

#include "algorithms/tracking/TrackVertexConstraint.h"
#include "algorithms/tracking/TrackVertexConstraintConfig.h"
#include "extensions/jana/JOmniFactory.h"

namespace eicrecon {

class TrackVertexConstraint_factory
    : public JOmniFactory<TrackVertexConstraint_factory, TrackVertexConstraintConfig> {

private:
  using AlgoT = eicrecon::TrackVertexConstraint;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::Track> m_tracks_input{this};
  PodioInput<edm4eic::MCRecoTrackParticleAssociation> m_track_assocs_input{this};
  PodioInput<edm4eic::Vertex> m_vertices_input{this};
  PodioInput<edm4eic::Track> m_vertex_tracks_input{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};
  PodioOutput<edm4eic::Trajectory> m_trajectories_output{this};
  PodioOutput<edm4eic::Track> m_tracks_output{this};
  PodioOutput<edm4eic::MCRecoTrackParticleAssociation> m_track_assocs_output{this};

  ParameterRef<std::string> m_beamSetting{
      this, "beamSetting", config().beamSetting,
      "beam setting of the luminous region with the beam energies of the geometry: hiDiv or hiAcc"};
  ParameterRef<double> m_beamSpotX{this, "beamSpotX", config().beamSpotX, "beam spot x [mm]"};
  ParameterRef<double> m_beamSpotY{this, "beamSpotY", config().beamSpotY, "beam spot y [mm]"};
  ParameterRef<double> m_beamSpotZ{this, "beamSpotZ", config().beamSpotZ, "beam spot z [mm]"};
  ParameterRef<double> m_beamSpotSigmaX{
      this, "beamSpotSigmaX", config().beamSpotSigmaX,
      "beam spot spread in x [mm], 0: that of the beam configuration"};
  ParameterRef<double> m_beamSpotSigmaY{
      this, "beamSpotSigmaY", config().beamSpotSigmaY,
      "beam spot spread in y [mm], 0: that of the beam configuration"};
  ParameterRef<double> m_beamSpotSigmaZ{
      this, "beamSpotSigmaZ", config().beamSpotSigmaZ,
      "beam spot spread in z [mm], 0: that of the beam configuration"};
  ParameterRef<double> m_maxChi2{this, "maxChi2", config().maxChi2,
                                 "max chi2 of a track's (d0, z0) with the vertex"};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void Process(int32_t /* run_number */, uint64_t /* event_number */) override {
    m_algo->process(
        {m_tracks_input(), m_track_assocs_input(), m_vertices_input(), m_vertex_tracks_input()},
        {m_trackparams_output().get(), m_trajectories_output().get(), m_tracks_output().get(),
         m_track_assocs_output().get()});
  }
};

} // namespace eicrecon
