// SPDX-License-Identifier: LGPL-3.0-or-later
#include "B0WeakPrior.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

void check(bool condition) { if (!condition) throw std::runtime_error("weak-prior test failed"); }
int main() {
  std::array<double, 6> variances{1e6, 1e6, 1, 1, 1, 1e6};
  const auto covariance = eicrecon::b0::weakPrior(variances);
  check(covariance(0, 0) == 1e6 && covariance(4, 4) == 1 && covariance(0, 4) == 0);
  check(eicrecon::b0::weakPrior(variances, 100)(2, 2) == 100);
  for (const double bad : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity()}) {
    bool caught = false;
    try { (void)eicrecon::b0::weakPrior(variances, bad); }
    catch (const std::invalid_argument&) { caught = true; }
    check(caught);
  }
  // Linear scalar toy: information accounting, not a B0 physics benchmark.
  const double information = 8 / .0004;
  const double measurementVariance = 1 / information;
  for (double prior : {1e3, 1e5, 1e7}) {
    const double fittedVariance = 1 / (1 / prior + information);
    check(std::abs(fittedVariance / measurementVariance - 1) < 1e-7);
  }
  const double reusedPosterior = 1 / (information + information);
  check(std::abs(reusedPosterior / measurementVariance - .5) < 1e-12);
  auto invalid = variances; invalid[3] = -1;
  bool caught = false;
  try { (void)eicrecon::b0::weakPrior(invalid); } catch (const std::invalid_argument&) { caught = true; }
  check(caught);
  std::cout << "weak initialization and information-accounting tests passed\n";
}
