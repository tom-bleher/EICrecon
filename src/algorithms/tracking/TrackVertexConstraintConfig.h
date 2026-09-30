// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <string>

namespace eicrecon {

struct TrackVertexConstraintConfig {
  /// Beam settings of the luminous region: "hiDiv" (high divergence) or "hiAcc" (high acceptance).
  /// With the beam energies of the geometry (constants ElectronBeamEnergy, IonBeamEnergy) they
  /// select the luminous region of the afterburner configuration.
  std::string beamSetting = "hiDiv";
  /// Luminous region centre [mm] (the centres of the beam configurations are at the origin) and
  /// spread [mm]; a spread of 0 takes that of the beam configuration
  double beamSpotX      = 0.;
  double beamSpotY      = 0.;
  double beamSpotZ      = 0.;
  double beamSpotSigmaX = 0.;
  double beamSpotSigmaY = 0.;
  double beamSpotSigmaZ = 0.;
  /// Maximum chi2 (2 degrees of freedom) of a track's (d0, z0) with the vertex, for a central
  /// track to enter the vertex and for a track to be constrained
  double maxChi2 = 9.;
};

} // namespace eicrecon
