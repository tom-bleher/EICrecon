// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <edm4eic/unit_system.h>
#include <cstddef>

namespace eicrecon {

struct B0TripletSeedingConfig {
  /// A z gap between consecutive hits larger than this starts a new station. It
  /// must exceed the z spread within a station (up to 12 mm) and stay below the
  /// station spacing (at least 270 mm).
  double stationGap = 50. * Acts::UnitConstants::mm;
  /// Maximum distance of the middle hit from the chord of the outer two,
  /// measured along the local magnetic field, after subtracting the part that
  /// the field integral predicts (see fieldSamples)
  double maxResidual = 1. * Acts::UnitConstants::mm;
  /// Minimum momentum of the three-hit estimate
  double minMomentum = 1. * Acts::UnitConstants::GeV;
  /// Distance upstream of the first hit at which the seed is expressed. The seed is moved there
  /// along a straight line, and its position error grows by the sagitta of the track over it.
  double anchorDistance = 10. * Acts::UnitConstants::mm;
  /// Maximum number of seeds per event. Candidates are ranked by residual alone and the rest are
  /// dropped, so in a busy event genuine seeds can be lost. A warning is logged the 1st, 2nd, 4th,
  /// 8th, ... time it happens.
  std::size_t maxSeeds = 1000;

  /// Hit-time window. A particle is never faster than light, so the later of two hits may only
  /// lag by up to the flight time of the slowest seeded particle, distance * (1 / beta_min - 1),
  /// and may not be early. Both bounds are widened by timeWindowSigmas times the combined hit time
  /// error (TrackerHit::getTimeError(), at least timeErrorFloor each) and by timeWindowTolerance.
  /// beta_min is that of a particle of mass slowParticleMass at minMomentum. With hit time errors
  /// of several ns (before the B0 timing was set to 30 ps) the window rejects nothing.
  /// 0 turns the window off.
  double timeWindowSigmas = 5.;
  /// in EDM4eic units (ns)
  double timeWindowTolerance = 0.;
  /// in EDM4eic units (ns)
  double timeErrorFloor = 0.01;
  /// Mass per unit charge of the slowest particles to seed. minMomentum is a cut on rigidity, so
  /// for nuclear fragments set the mass per unit charge: about 1.876 GeV for deuterons, alpha
  /// particles and carbon nuclei.
  double slowParticleMass = 0.938272 * Acts::UnitConstants::GeV;

  /// Number of points at which the field is sampled between the outer hits of a triplet. The
  /// field changes along the track (the B0 magnet has a quadrupole gradient and is not aligned
  /// with the ion beam), so the field at the middle hit biases the momentum by up to 5 % for
  /// triplets that skip a station, and the chord residual picks up a term from the change of the
  /// field direction. The samples, on the parabola through the hits and weighted by their effect
  /// on the middle hit, give the effective field for the momentum and the predicted residual,
  /// which is subtracted. They are two Gauss-Legendre points per interval on each side of the
  /// middle hit, so the number is rounded up to a multiple of 4; 4 agrees with 400 to 2 um in the
  /// B0 field. 0 uses the field at the middle hit only.
  std::size_t fieldSamples = 4;
  /// With fieldSamples, triplets are first kept if their raw residual is below maxResidual plus
  /// this, to admit tracks whose raw residual is mostly the predicted field term
  double fieldTermAllowance = 2. * Acts::UnitConstants::mm;

  /// Prior uncertainties of the seed parameters. The prior is a deliberately loose diagonal
  /// envelope, not an error estimate: the track finder refits from it, so it only has to cover
  /// the true error. The relative q/p error is a floor that dominates for most seeds. Widening
  /// the angle error to 10 or 30 mrad did not change the efficiency for 2-15 GeV pions. Seeds
  /// built from other particles or from stray hits can have the wrong charge.
  double positionError       = 0.1 * Acts::UnitConstants::mm;
  double angleError          = 2. * Acts::UnitConstants::mrad;
  double qOverPRelativeError = 0.25;
  /// in EDM4eic units (ns), unlike the other values: Acts measures time as a length
  double timeError = 10. * edm4eic::unit::ns;

  /// Particle hypothesis of the seeds, written to TrackParameters::pdg for the track finding,
  /// which uses its mass for the multiple scattering and energy loss and, with hit times, for the
  /// expected time of flight. The times of the three hits measure 1/beta to about 1.3 % with
  /// 30 ps hits, which separates protons from pions up to about 3 GeV. Most B0 tracks are
  /// protons, so a seed is a proton unless its 1/beta lies pionSigmas errors or more below that
  /// of a proton of its momentum, which makes it a pion. Above a few GeV, where they cannot be
  /// told apart, the proton mass still gives the right multiple scattering for protons and changes
  /// it by a few % at most for pions. No hypothesis is written (pdg 0, for which the track finding
  /// uses pions) if the error on 1/beta exceeds maxInverseBetaError, as with hit time errors of
  /// several ns, or if massFromTimeOfFlight is off. Like Track::pdg, which is the pion for every
  /// other track, this is the mass hypothesis of the track fit, not a particle identification.
  bool massFromTimeOfFlight  = true;
  double pionSigmas          = 3.;
  double maxInverseBetaError = 0.1;
};

} // namespace eicrecon
