// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <sstream>
#include <stdexcept>

#include "algorithms/tracking/SharedOffsetCurve.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("shared offset curve interpolates and clamps", "[SharedOffsetCurve]") {
  std::stringstream text("# f offset_mm\n0.0 0.05\n0.25 0.15\n0.5 0.25\n");
  const eicrecon::SharedOffsetCurve curve(text);
  CHECK_THAT(curve.offset(0.0), WithinAbs(0.05, 1e-12));
  CHECK_THAT(curve.offset(0.125), WithinAbs(0.10, 1e-12));
  CHECK_THAT(curve.offset(0.4), WithinAbs(0.21, 1e-12));
  CHECK_THAT(curve.offset(0.7), WithinAbs(0.25, 1e-12));
}

TEST_CASE("shared offset curve rejects unordered rows", "[SharedOffsetCurve]") {
  std::stringstream text("0.2 0.1\n0.1 0.2\n");
  CHECK_THROWS_AS(eicrecon::SharedOffsetCurve(text), std::runtime_error);
}
