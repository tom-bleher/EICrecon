// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Eigen/Core>
#include <cmath>

#include "algorithms/tracking/B0Seeding2.h"

using eicrecon::B0Seeding2;

TEST_CASE("B0Seeding2 rotates the hits into a frame along the field and the beam", "[B0Seeding2]") {
  // The ion beam at the crossing angle, and the B0 field off its midplane, where the quadrupole
  // term turns it about the beam
  const double crossing = -0.025;
  const Acts::Vector3 beam{std::sin(crossing), 0., std::cos(crossing)};
  const double gradient = -8.1224 * Acts::UnitConstants::T / Acts::UnitConstants::m;
  const double y        = 80. * Acts::UnitConstants::mm;
  const Acts::Vector3 field{gradient * y, 1.18 * Acts::UnitConstants::T, 0.};

  const auto rotation = B0Seeding2::frameRotation(field, beam);
  // A rotation
  CHECK((rotation * rotation.transpose()).isApprox(Acts::RotationMatrix3::Identity(), 1e-12));
  CHECK(rotation.determinant() == Catch::Approx(1.));
  // z' along the field, x' along the beam as far as it is across the field
  CHECK((rotation * field.normalized()).isApprox(Acts::Vector3::UnitZ(), 1e-12));
  const Acts::Vector3 beamPrimed = rotation * beam;
  CHECK(beamPrimed.y() == Catch::Approx(0.).margin(1e-12));
  CHECK(beamPrimed.x() > 0.99);
  // so a hit 6 m along the beam is at r' = 6 m, the stations become layers in r'
  const Acts::Vector3 hit = rotation * (6. * Acts::UnitConstants::m * beam);
  CHECK(std::hypot(hit.x(), hit.y()) ==
        Catch::Approx(6. * Acts::UnitConstants::m * std::abs(beamPrimed.x())).epsilon(1e-12));
}
