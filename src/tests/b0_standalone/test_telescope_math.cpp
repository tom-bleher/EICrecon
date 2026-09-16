// SPDX-License-Identifier: LGPL-3.0-or-later
#include "B0TelescopeMath.h"
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace eicrecon::b0;
void check(bool pass, const char* why) { if (!pass) throw std::runtime_error(why); }
std::vector<TelescopeHit> fixture(double q = .05) {
  std::vector<TelescopeHit> hits;
  for (std::size_t i = 0; i < 4; ++i) {
    const double z = 250 * i;
    TelescopeHit h;
    h.index = i; h.station = i; h.surface = 100 + i;
    h.position = {40 + .02 * z - .5 * bendConstant * q * z * z, 20 + .01 * z, 6000 + z};
    h.covariance.diagonal() << .0004, .0004, 0;
    h.time = h.position.norm() / speedOfLight; h.timeVariance = .0009;
    hits.push_back(h);
  }
  return hits;
}
int main() {
  Field dipole = [](const Eigen::Vector3d&) { return Eigen::Vector3d(0, 1, 0); };
  for (const bool integral : {false, true}) {
    auto fit = fitTelescope(fixture(), dipole, {integral, 4, 16});
    check(fit.valid, "dipole fit");
    check(std::abs(fit.state(4) - .05) < .0001, "curvature units/sign");
    check(std::abs(fit.state(2) - .02) < .0001, "direction");
    check(fit.covariance.llt().info() == Eigen::Success, "positive covariance");
  }
  auto zero = fitTelescope(fixture(0), dipole);
  zero.state(4) = 0;
  const auto both = chargeProposals(zero);
  check(both.size() == 2 && both[0].state(4) < 0 && both[1].state(4) > 0, "finite zero-curvature proposals");
  check((both[0].state.head<4>() - both[1].state.head<4>()).norm() > 0, "charge changes correlated state");
  check(chargeProposals(zero, 3, 10000, -1).size() == 1, "forced sign");
  auto noField = fitTelescope(fixture(), [](const Eigen::Vector3d&) { return Eigen::Vector3d::Zero(); });
  check(!noField.valid, "unmeasured q/p is not a precision estimate");
  auto bad = fixture(); bad[0].covariance(0, 0) = -1;
  check(!fitTelescope(bad, dipole).valid, "invalid covariance");
  auto displaced = fixture();
  for (auto& hit : displaced) hit.position.x() += 80;
  check(fitTelescope(displaced, dipole).valid, "no hidden origin cut");
  auto hits = fixture();
  const double ds = (hits[1].position - hits[0].position).norm();
  hits[0].time = 43; hits[1].time = 43 + ds / speedOfLight;
  check(timeCompatible(hits[0], hits[1], 0, .938, 5, .01), "unknown production time cancels");
  hits[1].time += 2;
  check(!timeCompatible(hits[0], hits[1], 0, .938, 5, .01), "out-of-time rejection");
  hits[1].time = std::numeric_limits<double>::quiet_NaN();
  check(timeCompatible(hits[0], hits[1], 0, .938, 5, .01), "missing timing is explicit permissive policy");
  auto result = findTelescopeCandidates(fixture(), dipole);
  check(!result.candidates.empty() && result.candidates.front().stations == 4, "four-station extension");
  auto reversed = fixture(); std::reverse(reversed.begin(), reversed.end());
  auto reordered = findTelescopeCandidates(reversed, dipole);
  check(result.candidates.front().hits == reordered.candidates.front().hits, "order-independent association");
  SearchOptions limited; limited.maxTrials = 2;
  auto truncated = findTelescopeCandidates(fixture(), dipole, limited);
  check(truncated.trials == 2 && truncated.truncated, "hard endpoint/candidate budget");
  SearchOptions minimum; minimum.minStations = 4;
  auto three = fixture(); three.pop_back();
  check(findTelescopeCandidates(three, dipole, minimum).candidates.empty(), "distinct-station minimum");
  auto doubled = fixture();
  auto copy = doubled[0]; copy.index = 90; copy.surface = 190; copy.position.z() += 1;
  copy.position.x() += .02; copy.position.y() += .01; doubled.push_back(copy);
  auto extended = findTelescopeCandidates(doubled, dipole);
  check(extended.candidates.front().hits.size() == 5, "front/back are independent measurements");
  std::cout << "B0 telescope mathematics and search tests passed\n";
}
