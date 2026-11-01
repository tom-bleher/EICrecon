// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <limits>

#include "algorithms/tracking/B0TrackSeeding.h"
#include "algorithms/tracking/B0TrackSeedingConfig.h"

using eicrecon::B0TrackSeeding;
using eicrecon::B0TrackSeedingConfig;

namespace {

const Acts::Vector3 kField{0., 1.3 * Acts::UnitConstants::T, 0.};

/// A hit at a position, in time with a particle of the given 1/beta from z = 5850 mm, with 30 ps
/// time errors (Acts measures time as a length)
B0TrackSeeding::Hit hit(const Acts::Vector3& position, double inverseBeta = 1.) {
  return {.position  = position,
          .field     = kField,
          .time      = (position.z() - 5850.) * inverseBeta,
          .timeError = 0.03 * Acts::UnitConstants::ns,
          .variance  = 4e-4,
          .index     = 0};
}

/// Position at z of a track of charge q and momentum p starting at (0, 0, 5800) along +z in kField;
/// the force q v x B points along -x for q > 0
Acts::Vector3 helix(double z, double q, double p) {
  const double radius = p / kField.norm();
  const double dz     = z - 5800.;
  return {-q * (radius - std::sqrt(radius * radius - dz * dz)), 0., z};
}

} // namespace

TEST_CASE("B0TrackSeeding rotates the hits into a frame along the field and the beam",
          "[B0TrackSeeding]") {
  // The ion beam at the crossing angle, and the B0 field off its midplane, where the quadrupole
  // term turns it about the beam
  const double crossing = -0.025;
  const Acts::Vector3 beam{std::sin(crossing), 0., std::cos(crossing)};
  const Acts::Vector3 field{-8.1224 * Acts::UnitConstants::T / Acts::UnitConstants::m * 80.,
                            1.18 * Acts::UnitConstants::T, 0.};
  const auto rotation = B0TrackSeeding::frameRotation(field, beam);
  CHECK((rotation * rotation.transpose()).isApprox(Acts::RotationMatrix3::Identity(), 1e-12));
  CHECK(rotation.determinant() == Catch::Approx(1.));
  // z' along the field, x' along the beam as far as it is across the field
  CHECK((rotation * field.normalized()).isApprox(Acts::Vector3::UnitZ(), 1e-12));
  const Acts::Vector3 beamPrimed = rotation * beam;
  CHECK(beamPrimed.y() == Catch::Approx(0.).margin(1e-12));
  CHECK(beamPrimed.x() > 0.99);
}

TEST_CASE("B0TrackSeeding keeps hits within the time window", "[B0TrackSeeding]") {
  B0TrackSeedingConfig cfg;
  const auto first = hit({0., 0., 5850.});
  // Over 1030 mm: light passes, a 1.2 GeV proton (lag 0.92 ns) passes, the window of a 1 GeV
  // proton ends at 1.28 ns plus 5 sigma of 0.042 ns
  CHECK(B0TrackSeeding::inTimeWindow(cfg, first, hit({0., 0., 6880.})));
  auto lagged = [&](double ns) {
    auto h = hit({0., 0., 6880.});
    h.time += ns * Acts::UnitConstants::ns;
    return h;
  };
  CHECK(B0TrackSeeding::inTimeWindow(cfg, first, lagged(0.92)));
  CHECK(B0TrackSeeding::inTimeWindow(cfg, first, lagged(1.45)));
  CHECK_FALSE(B0TrackSeeding::inTimeWindow(cfg, first, lagged(1.55)));
  // Nothing is faster than light; hits from another bunch crossing are rejected
  CHECK_FALSE(B0TrackSeeding::inTimeWindow(cfg, first, lagged(-0.3)));
  CHECK_FALSE(B0TrackSeeding::inTimeWindow(cfg, first, lagged(10.)));
  // With hit time errors of several ns the window rejects nothing
  auto coarse      = lagged(-10.);
  coarse.timeError = 8. * Acts::UnitConstants::ns;
  CHECK(B0TrackSeeding::inTimeWindow(cfg, first, coarse));
  // A hit without a valid time error is never rejected, and the window can be turned off
  auto invalid      = lagged(30.);
  invalid.timeError = std::numeric_limits<double>::quiet_NaN();
  CHECK(B0TrackSeeding::inTimeWindow(cfg, first, invalid));
  auto off             = cfg;
  off.timeWindowSigmas = 0.;
  CHECK(B0TrackSeeding::inTimeWindow(off, first, lagged(-30.)));
}

TEST_CASE("B0TrackSeeding recovers charge and momentum on a helix", "[B0TrackSeeding]") {
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  B0TrackSeedingConfig cfg;
  for (const double q : {+1., -1.}) {
    for (const double p : {5., 20., 100.}) {
      const auto seed = B0TrackSeeding::makeSeed(cfg, gctx, hit(helix(5900., q, p)),
                                                 hit(helix(6270., q, p)), hit(helix(6880., q, p)));
      REQUIRE(seed.has_value());
      CHECK(seed->parameters[Acts::eBoundQOverP] * Acts::UnitConstants::GeV ==
            Catch::Approx(q / p).epsilon(0.01));
      // On a perigee anchorDistance upstream of the first hit, with a positive prior
      CHECK(seed->anchor.z() == Catch::Approx(5900. - cfg.anchorDistance).margin(0.1));
      CHECK((seed->errors.array() > 0.).all());
    }
  }
  // Below minMomentum there is no seed
  CHECK_FALSE(B0TrackSeeding::makeSeed(cfg, gctx, hit(helix(5900., 1., 0.5)),
                                       hit(helix(6270., 1., 0.5)), hit(helix(6880., 1., 0.5)))
                  .has_value());
}
