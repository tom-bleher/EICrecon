// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <sstream>

#include "algorithms/digi/PadResponseTable.h"
#include "algorithms/digi/TrapezoidPadResponse.h"

using Catch::Matchers::WithinAbs;
using eicrecon::trapezoidPadShare;
using eicrecon::trapezoidSharedOffset;

TEST_CASE("trapezoid pad response", "[TrapezoidPadResponse]") {
  const double pitch     = 0.5;
  const double electrode = GENERATE(0.0, 0.1, 0.15, 0.3);

  SECTION("metal collects everything, mid-gap splits evenly, neighbour metal gets nothing") {
    CHECK_THAT(trapezoidPadShare(0.5 * electrode, pitch, electrode), WithinAbs(1.0, 1e-12));
    CHECK_THAT(trapezoidPadShare(0.5 * pitch, pitch, electrode), WithinAbs(0.5, 1e-12));
    CHECK_THAT(trapezoidPadShare(pitch - 0.5 * electrode, pitch, electrode), WithinAbs(0.0, 1e-12));
    CHECK(trapezoidPadShare(1.5 * pitch, pitch, electrode) == 0.0);
  }

  SECTION("shares of a pad row sum to one and invert back to the deposit position") {
    const double offset = GENERATE(0.0, 0.03, 0.07, 0.12, 0.2, 0.249);
    double sum          = 0;
    for (int k = -2; k <= 2; ++k) {
      sum += trapezoidPadShare(offset - k * pitch, pitch, electrode);
    }
    CHECK_THAT(sum, WithinAbs(1.0, 1e-12));

    const double lead      = trapezoidPadShare(offset, pitch, electrode);
    const double neighbour = trapezoidPadShare(offset - pitch, pitch, electrode);
    if (neighbour > 0) {
      const double f = neighbour / (lead + neighbour);
      CHECK_THAT(trapezoidSharedOffset(f, pitch, electrode), WithinAbs(offset, 1e-12));
    }
  }
}

TEST_CASE("electrode as wide as the pitch means no sharing", "[TrapezoidPadResponse]") {
  CHECK(trapezoidPadShare(0.2, 0.5, 0.5) == 1.0);
  CHECK(trapezoidPadShare(0.3, 0.5, 0.5) == 0.0);
}

TEST_CASE("signal loss grows with the distance from the metal", "[TrapezoidPadResponse]") {
  using eicrecon::trapezoidSignalFraction;
  const double loss = 0.5;
  CHECK_THAT(trapezoidSignalFraction(1.0, 1.0, loss), WithinAbs(1.0, 1e-12)); // on the metal
  CHECK_THAT(trapezoidSignalFraction(0.5, 1.0, loss), WithinAbs(0.5, 1e-12)); // mid-gap
  CHECK_THAT(trapezoidSignalFraction(0.5, 0.5, loss), WithinAbs(0.5, 1e-12)); // corner, saturated
  CHECK_THAT(trapezoidSignalFraction(0.25, 1.0, loss), WithinAbs(0.75, 1e-12));
  CHECK_THAT(trapezoidSignalFraction(0.875, 1.0, loss), WithinAbs(0.875, 1e-12));
}

TEST_CASE("a pad response table reproduces the function it was filled from", "[PadResponseTable]") {
  const double pitch = 0.5, electrode = 0.1;
  const int n = 50, radius = 1;
  std::stringstream text;
  text << "# trapezoid, 100 um electrode\n"
       << pitch << " " << pitch << " " << n << " " << n << " " << radius << "\n";
  for (int ix = 0; ix < n; ++ix) {
    for (int iy = 0; iy < n; ++iy) {
      const double ux = (ix + 0.5) / n * pitch - pitch / 2, uy = (iy + 0.5) / n * pitch - pitch / 2;
      for (int di = -radius; di <= radius; ++di) {
        for (int dj = -radius; dj <= radius; ++dj) {
          text << trapezoidPadShare(ux - di * pitch, pitch, electrode) *
                      trapezoidPadShare(uy - dj * pitch, pitch, electrode)
               << " ";
        }
      }
      text << "\n";
    }
  }
  const eicrecon::PadResponseTable table(text);
  // Deposit offsets from the pad centre, away from the kinks at the metal edge
  for (const double dx : {0.0, 0.12, 0.2, 0.3, -0.35, 0.6}) {
    for (const double dy : {0.0, -0.15, 0.4}) {
      const double expected =
          trapezoidPadShare(dx, pitch, electrode) * trapezoidPadShare(dy, pitch, electrode);
      CHECK_THAT(table.fraction(dx, dy), WithinAbs(expected, 1e-3));
    }
  }
  CHECK(table.fraction(1.2, 0.0) == 0.0); // beyond the tabulated neighbourhood
}
