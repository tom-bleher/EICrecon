// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/Cov6f.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4hep/Vector3f.h>
#include <cmath>
#include <initializer_list>
#include <optional>

#include "algorithms/tracking/B0TrackSeeding.h"
#include "algorithms/tracking/B0TrackSeedingConfig.h"

using eicrecon::B0TrackSeeding;
using eicrecon::B0TrackSeedingConfig;

namespace {

const Acts::Vector3 kField{0., 1.3 * Acts::UnitConstants::T, 0.};

/// Position at z of a track of charge q and momentum p starting at (0, 0, 5800) along +z in kField;
/// the force q v x B points along -x for q > 0
Acts::Vector3 helix(double z, double q, double p) {
  const double radius = p / kField.norm();
  const double dz     = z - 5800.;
  return {-q * (radius - std::sqrt(radius * radius - dz * dz)), 0., z};
}

} // namespace

TEST_CASE("B0TrackSeeding recovers charge and momentum on a helix", "[B0TrackSeeding]") {
  const B0TrackSeedingConfig cfg;
  for (const double q : {+1., -1.}) {
    for (const double p : {5., 20., 100.}) {
      edm4eic::TrackSeedCollection seeds;
      edm4eic::TrackParametersCollection parameters;
      const auto seed = B0TrackSeeding::addSeed(cfg, helix(5900., q, p), helix(6270., q, p),
                                                helix(6880., q, p), kField, 3.F, seeds, parameters);
      REQUIRE(seed.has_value());
      CHECK(parameters[0].getQOverP() == Catch::Approx(q / p).epsilon(0.01));
      CHECK(parameters[0].getTime() == Catch::Approx(3.));
      CHECK(seeds[0].getPerigee().z == Catch::Approx(5900. - cfg.anchorDistance).margin(0.1));
      CHECK(std::sqrt(parameters[0].getCovariance()(4, 4)) ==
            Catch::Approx(cfg.qOverPRelativeError / p).epsilon(0.01));
    }
  }
}
