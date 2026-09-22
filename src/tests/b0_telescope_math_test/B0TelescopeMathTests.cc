// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#include "B0TelescopeMath.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
std::size_t checks = 0;
void check(bool condition, const char* message) {
  ++checks;
  if (!condition) {
    throw std::runtime_error(message);
  }
}
bool close(double a, double b, double tolerance = 1.e-10) {
  return std::abs(a - b) <= tolerance;
}
}

int main() {
  using namespace eicrecon::b0telescope;
  try {
    for (int layer = 1; layer <= 8; ++layer) {
      check(stationFromLayer(layer, 2) == static_cast<unsigned>((layer - 1) / 2),
             "realistic front/back station mapping");
    }
    for (int layer = 1; layer <= 4; ++layer) {
      check(stationFromLayer(layer, 1) == static_cast<unsigned>(layer - 1), "ideal station mapping");
    }
    for (int layer : {-1, 0, 9, 100}) {
      check(!stationFromLayer(layer, 2), "invalid layer must not wrap into a station");
    }
    check(!stationFromLayer(1, 0), "zero layers per station");
    check(!stationFromLayer(1, 3), "unsupported layer layout");
    check(stationTriples.size() == 4, "all four leave-one-station-out triples");
    for (const auto& t : stationTriples) {
      check(t[0] < t[1] && t[1] < t[2], "three distinct ordered physical stations");
    }
    const Point lab{-160., 30., 6100.};
    const auto local = toTelescopeFrame(lab, -0.025);
    const auto restored = toTelescopeFrame(local, 0.025);
    check(close(lab.x, restored.x) && close(lab.z, restored.z), "rotation round trip");
    check(close(lab.y, local.y), "rotation preserves y");
    check(close(std::hypot(lab.x, lab.z), std::hypot(local.x, local.z)), "rotation preserves length");
    check(compatibleDoublet({-150., 10., 6000.}, {-155., 11., 6250.}, 50., 0.1),
           "accept a track with increasing radius");
    check(compatibleDoublet({-150., 10., 6000.}, {-145., 11., 6250.}, 50., 0.1),
           "accept a track with decreasing radius");
    check(!compatibleDoublet({-150., 10., 6000.}, {-149., 11., 6007.}, 50., 0.1),
           "front/back hits alone cannot make a telescope doublet");
    check(!compatibleDoublet({0., 0., 6000.}, {0., 0., 6000.}, 50., 0.1), "zero lever arm");
    check(!compatibleDoublet({0., 0., 6250.}, {0., 0., 6000.}, 50., 0.1), "wrong ordering");
    check(!compatibleDoublet({0., 0., 6000.}, {100., 0., 6250.}, 50., 0.1), "slope rejection");
    const Point a{1000., 2000., 6000.}, b{1003., 2006., 6300.}, c{1008., 2016., 6800.};
    const auto displaced = tripletScore(a, b, c, 20., 5.);
    check(displaced && close(*displaced, 0.), "do not impose an interaction-point cut");
    check(!tripletScore(a, {1003., 2012., 6300.}, c, 20., 5.), "non-bending mismatch");
    check(!tripletScore(a, {1030., 2006., 6300.}, c, 20., 5.), "bending mismatch");
    check(!tripletScore(a, c, b, 20., 5.), "triplet station ordering");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    check(!tripletScore(a, {nan, 0., 6300.}, c, 20., 5.), "NaN position");
    check(!tripletScore(a, {0., inf, 6300.}, c, 20., 5.), "infinite position");
    check(!tripletScore(a, b, c, 0., 5.), "zero residual window");
    for (int i = 1; i <= 1000; ++i) {
      const double curvature = i * 1.e-8;
      const double dz = 750.;
      const double f = 276.086 / dz;
      const Point p0{-150., 7., 6099.75};
      const Point pm{p0.x + 0.01 * f * dz + curvature * f * f * dz * dz,
                      p0.y + 0.003 * f * dz, p0.z + f * dz};
      const Point p2{p0.x + 0.01 * dz + curvature * dz * dz,
                      p0.y + 0.003 * dz, p0.z + dz};
      const auto score = tripletScore(p0, pm, p2, 20., 5.);
      const double expected = std::pow(curvature * f * (f - 1.) * dz * dz / 20., 2);
      check(score && close(*score, expected), "quadratic sagitta on unequally spaced stations");
      const auto mirrored = tripletScore({-p0.x, p0.y, p0.z}, {-pm.x, pm.y, pm.z},
                                          {-p2.x, p2.y, p2.z}, 20., 5.);
      check(mirrored && close(*score, *mirrored), "sagitta selection is charge-sign symmetric");
    }
    check(withinCombinationBudget({0, 0, 0, 0}, 0), "empty event");
    check(withinCombinationBudget({1, 1, 1, 1}, 4), "exact four-triplet budget");
    check(!withinCombinationBudget({1, 1, 1, 1}, 3), "do not silently return partial station coverage");
    check(withinCombinationBudget({2, 2, 2, 2}, 32), "eight face choices per station triple");
    check(!withinCombinationBudget({2, 2, 2, 2}, 31), "combinatorial cap");
    check(withinCombinationBudget({1, 0, 1, 1}, 1), "missing station still has one valid triple");
    check(withinCombinationBudget({1, 1, 0, 0}, 0), "two stations cannot seed");
    const auto huge = std::numeric_limits<std::size_t>::max();
    check(!withinCombinationBudget({huge, huge, huge, huge}, 100000), "overflow-safe budget");
    check(withinCombinationBudget({huge, huge, 0, 0}, 0), "empty product does not overflow");
    std::cout << checks << " B0 telescope geometry checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << "B0 telescope geometry test failed: " << error.what() << '\n';
    return 1;
  }
}
