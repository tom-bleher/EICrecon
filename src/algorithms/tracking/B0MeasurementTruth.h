// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace eicrecon::b0 {
using TruthId = std::pair<int, int>;
using TruthVotes = std::map<TruthId, double>;
using RawTruthVotes = std::map<TruthId, TruthVotes>;
template <typename Object> TruthId truthId(const Object& object) {
  const auto id = object.getObjectID(); return {id.collectionID, id.index};
}
// Use measurement channel weights, retain unassociated noise in the purity
// denominator, and normalize legacy unnormalized associations only downward.
// Never promote a tiny AC response tail to a unit-weight truth match.
template <typename Measurement>
TruthVotes measurementTruth(const Measurement& measurement, const RawTruthVotes& rawTruth) {
  TruthVotes result;
  if (measurement.hits_size() == 0) return result;
  const auto weights = measurement.getWeights();
  if (!weights.empty() && weights.size() != measurement.hits_size())
    throw std::invalid_argument("B0 measurement hits/weights size mismatch");
  double denominator = weights.empty() ? double(measurement.hits_size()) : 0;
  if (!weights.empty()) for (double w : weights) {
    if (!std::isfinite(w) || w < 0) throw std::invalid_argument("Invalid B0 measurement weight");
    denominator += w;
  }
  if (!(denominator > 0) || !std::isfinite(denominator)) throw std::invalid_argument("Invalid B0 measurement weight sum");
  std::set<TruthId> seen;
  std::size_t i = 0;
  for (const auto& hit : measurement.getHits()) {
    if (!seen.insert(truthId(hit)).second) throw std::invalid_argument("Duplicate B0 measurement hit");
    const double fraction = (weights.empty() ? 1 : weights[i]) / denominator; ++i;
    const auto found = rawTruth.find(truthId(hit.getRawHit()));
    if (found == rawTruth.end()) continue;
    double sum = 0;
    for (const auto& [id, w] : found->second) {
      (void)id;
      if (!std::isfinite(w) || w < 0) throw std::invalid_argument("Invalid B0 raw truth weight");
      sum += w;
    }
    if (!std::isfinite(sum)) throw std::invalid_argument("Invalid B0 raw truth sum");
    for (const auto& [id, w] : found->second) result[id] += fraction * w / std::max(1.0, sum);
  }
  return result;
}
} // namespace eicrecon::b0
