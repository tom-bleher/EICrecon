// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/MCRecoTrackParticleAssociationCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrajectoryCollection.h>
#include <edm4eic/VertexCollection.h>
#include <Eigen/Core>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/TrackVertexConstraintConfig.h"

namespace eicrecon {

using TrackVertexConstraintAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::TrackCollection, edm4eic::MCRecoTrackParticleAssociationCollection,
                      edm4eic::VertexCollection, edm4eic::TrackCollection>,
    algorithms::Output<edm4eic::TrackParametersCollection, edm4eic::TrajectoryCollection,
                       edm4eic::TrackCollection,
                       edm4eic::MCRecoTrackParticleAssociationCollection>>;

/// Constrains tracks with the event vertex.
///
/// The vertex starts as the beam spot, the luminous region of the beam configuration, and is
/// updated with the first reconstructed vertex or, without one, with the (d0, z0) of the
/// compatible vertex tracks (e.g. central tracks), best measured first. Each track is then updated
/// with the vertex projected onto its (d0, z0) at the perigee of the origin. Tracks measured far
/// from the vertex, like those of the B0 tracker, have (d0, z0) uncertain by millimetres, and the
/// vertex fixes their direction at the interaction point.
///
/// A track incompatible with the vertex, or not expressed at the perigee of the origin, is copied
/// unchanged. Tracks from decays that point back to the vertex (e.g. the proton of a forward
/// Lambda) pass the compatibility test, so decay analyses should use the unconstrained tracks.
///
/// The outputs have one entry per input track, in order: the parameters, trajectories that carry
/// them (otherwise copies of the input trajectories, as the parameters of a track are read from its
/// trajectory downstream), the tracks (momentum and position from the constrained parameters; chi2
/// and ndf include the vertex; the input track as their segment), and the MC associations of the
/// input tracks pointing to the new ones.
class TrackVertexConstraint : public TrackVertexConstraintAlgorithm,
                              public WithPodConfig<TrackVertexConstraintConfig> {
public:
  TrackVertexConstraint(std::string_view name)
      : TrackVertexConstraintAlgorithm{
            name,
            {"inputTracks", "inputTrackAssociations", "inputVertices", "inputVertexTracks"},
            {"outputTrackParameters", "outputTrajectories", "outputTracks",
             "outputTrackAssociations"},
            "constrain tracks with the beam spot and the event vertex"} {}

  /// Resolves the beam spot: the configured spreads, or those of the beam configuration of the
  /// geometry
  void init() final;
  void process(const Input& input, const Output& output) const final;

  /// Vertex estimate (x, y, z) and its covariance
  struct Vertex {
    Eigen::Vector3d position;
    Eigen::Matrix3d covariance;
  };

  /// The beam spot, updated with the first vertex or the compatible vertex tracks
  Vertex estimateVertex(const edm4eic::VertexCollection& vertices,
                        const edm4eic::TrackCollection& vertexTracks) const;

  /// The track's parameters at the perigee of the origin updated with the vertex, and the chi2 of
  /// the update, or nothing if the track is not at that perigee or not compatible with the vertex
  std::optional<std::pair<edm4eic::MutableTrackParameters, double>>
  constrain(const edm4eic::Track& track, const Vertex& vertex) const;

  /// Luminous region spread of the afterburner beam configurations at IP6, from 5000 events per
  /// configuration; nothing for an unknown configuration
  struct BeamSpot {
    Eigen::Vector3d sigma;
  };
  static std::optional<BeamSpot> beamSpot(double electronEnergy, double ionEnergy,
                                          const std::string& setting);

private:
  Eigen::Vector3d m_spot{Eigen::Vector3d::Zero()};
  Eigen::Vector3d m_spotSigma{Eigen::Vector3d::Zero()};
};

} // namespace eicrecon
