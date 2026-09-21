// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <limits>

#include "algorithms/tracking/CKFTracking.h"
#include "algorithms/tracking/CKFTrackingConfig.h"

TEST_CASE("CKF fit hypotheses preserve seed charge and momentum", "[tracking][ckfhypothesis]") {
  CHECK(eicrecon::CKFTrackingConfig{}.particleHypothesisPdg == 211);
  const auto surface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  for (const int pdg : {11, 13, 211, 321, 2212}) {
    const auto hypothesis = eicrecon::CKFTracking::makeParticleHypothesis(pdg);
    CHECK(static_cast<int>(hypothesis.absolutePdg()) == pdg);
    CHECK(hypothesis.mass() > 0.0);
    CHECK(hypothesis.absoluteCharge() == Acts::UnitConstants::e);
    for (const double charge : {-1., 1.}) {
      auto parameters                = Acts::BoundVector::Zero().eval();
      parameters[Acts::eBoundTheta]  = 0.03;
      parameters[Acts::eBoundQOverP] = charge / (15. * Acts::UnitConstants::GeV);
      const ActsExamples::TrackParameters track(surface, parameters, std::nullopt, hypothesis);
      CHECK(track.charge() == charge * Acts::UnitConstants::e);
      CHECK(track.absoluteMomentum() == Catch::Approx(15. * Acts::UnitConstants::GeV));
      CHECK(track.particleHypothesis().absolutePdg() == hypothesis.absolutePdg());
    }
  }
  CHECK(eicrecon::CKFTracking::makeParticleHypothesis(211).mass() ==
        Acts::ParticleHypothesis::pion().mass());
  CHECK(eicrecon::CKFTracking::makeParticleHypothesis(2212).mass() ==
        Acts::ParticleHypothesis::proton().mass());
  CHECK(eicrecon::CKFTracking::makeParticleHypothesis(2212).mass() >
        eicrecon::CKFTracking::makeParticleHypothesis(211).mass());
}

TEST_CASE("CKF rejects unsupported signed and neutral hypotheses", "[tracking][ckfhypothesis]") {
  for (const int pdg : {0, -11, -211, -2212, 22, 111, 2112, 999999}) {
    CHECK_THROWS_AS(eicrecon::CKFTracking::makeParticleHypothesis(pdg), std::invalid_argument);
  }
}

TEST_CASE("CKF refit covariance scale rejects invalid priors", "[tracking][ckfhypothesis]") {
  CHECK_NOTHROW(eicrecon::CKFTracking::validateRefitCovarianceScale(0.0));
  CHECK_NOTHROW(eicrecon::CKFTracking::validateRefitCovarianceScale(1.0));
  CHECK_NOTHROW(eicrecon::CKFTracking::validateRefitCovarianceScale(100.0));
  CHECK_THROWS_AS(eicrecon::CKFTracking::validateRefitCovarianceScale(-1.0), std::invalid_argument);
  CHECK_THROWS_AS(eicrecon::CKFTracking::validateRefitCovarianceScale(0.5), std::invalid_argument);
  CHECK_THROWS_AS(
      eicrecon::CKFTracking::validateRefitCovarianceScale(std::numeric_limits<double>::infinity()),
      std::invalid_argument);
  CHECK_THROWS_AS(
      eicrecon::CKFTracking::validateRefitCovarianceScale(std::numeric_limits<double>::quiet_NaN()),
      std::invalid_argument);
}
