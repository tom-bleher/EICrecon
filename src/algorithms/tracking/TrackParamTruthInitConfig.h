// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <Evaluator/DD4hepUnits.h>
#include <cmath>

struct TrackParamTruthInitConfig {

  double maxVertexX     = 80 * dd4hep::mm;
  double maxVertexY     = 80 * dd4hep::mm;
  double maxVertexZ     = 200 * dd4hep::mm;
  double minMomentum    = 100 * dd4hep::MeV;
  double maxEtaForward  = 6.0;
  double maxEtaBackward = 4.1;
  double momentumSmear  = 0.1;
  /** When true the seed perigee is placed at the MCParticle vertex instead
   *  of being extrapolated back to the beam axis.  This is essential for
   *  particles produced at large |z| (e.g. Lambda decay daughters) where
   *  the standard back-extrapolation puts the reference point far behind
   *  the IP, forcing ACTS to propagate many metres before reaching the
   *  detector, which degrades or eliminates CKF efficiency.
   */
  bool useVertexAsPerigee = false;
  /** When true the seed starts at the particle's earliest hit in the input sim hit collections
   *  instead of at its vertex: anchorDistance before that hit along its momentum there, with
   *  that momentum. Particles without such a hit get no seed. This keeps the CKF from having to
   *  propagate a loose prior over metres of field and material before the first measurement
   *  (e.g. ~6 m from the IP to the B0 tracker). Requires at least one sim hit collection input.
   */
  bool anchorAtFirstHit = false;
  double anchorDistance = 10 * dd4hep::mm;
  // Prior uncertainties of the seed (Gaussian sigmas)
  double locError    = 1 * dd4hep::mm;
  double phiError    = std::sqrt(0.05) * dd4hep::rad;
  double thetaError  = std::sqrt(0.01) * dd4hep::rad;
  double qOverPError = std::sqrt(0.1) / dd4hep::GeV;
};
