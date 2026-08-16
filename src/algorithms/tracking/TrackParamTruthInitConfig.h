// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <Evaluator/DD4hepUnits.h>

struct TrackParamTruthInitConfig {

  double maxVertexX     = 80 * dd4hep::mm;
  double maxVertexY     = 80 * dd4hep::mm;
  double maxVertexZ     = 200 * dd4hep::mm;
  double minMomentum    = 100 * dd4hep::MeV;
  double maxEtaForward  = 6.0;
  double maxEtaBackward = 4.1;
  double momentumSmear  = 0.1;

  // Seed covariance (diagonal; edm4eic units: mm^2, rad^2, (e/GeV)^2, ns^2).
  // The loc entries pin the seed to the true vertex at the mm level, which is
  // information no hit-based seeder has; widen them to the beam-spot /
  // bunch-length scale for a like-for-like reference.
  double locaError   = 1.0;
  double locbError   = 1.0;
  double phiError    = 0.05;
  double thetaError  = 0.01;
  double qOverPError = 0.1;
  double timeError   = 10e9;
};
