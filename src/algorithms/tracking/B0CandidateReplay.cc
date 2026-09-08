// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include "algorithms/tracking/B0CandidateReplay.h"

#include <fstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

namespace eicrecon {

B0ReplayFile loadB0ReplayFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::invalid_argument("B0 candidate replay: cannot open " + path);
  }
  nlohmann::json data;
  try {
    in >> data;
  } catch (const std::exception& exc) {
    throw std::invalid_argument(std::string("B0 candidate replay: invalid JSON: ") + exc.what());
  }
  if (!data.is_object() || data.value("format_version", 0) != 1) {
    throw std::invalid_argument("B0 candidate replay: expected format_version 1 object");
  }
  if (!data.contains("events") || !data["events"].is_array()) {
    throw std::invalid_argument("B0 candidate replay: missing events array");
  }
  B0ReplayFile out;
  out.format_version = 1;
  for (const auto& event : data["events"]) {
    if (!event.is_object() || !event.contains("run") || !event.contains("event") ||
        !event.contains("candidates")) {
      throw std::invalid_argument("B0 candidate replay: malformed event");
    }
    B0ReplayEvent rec;
    rec.run        = event["run"].get<std::int64_t>();
    rec.event      = event["event"].get<std::int64_t>();
    const auto key = std::make_pair(rec.run, rec.event);
    if (out.events.contains(key)) {
      throw std::invalid_argument("B0 candidate replay: duplicate run/event (" +
                                  std::to_string(rec.run) + ", " + std::to_string(rec.event) + ")");
    }
    if (!event["candidates"].is_array()) {
      throw std::invalid_argument("B0 candidate replay: candidates must be an array");
    }
    for (const auto& candidate : event["candidates"]) {
      if (!candidate.is_object() || !candidate.contains("hits") || !candidate["hits"].is_array()) {
        throw std::invalid_argument("B0 candidate replay: malformed candidate");
      }
      std::vector<B0ReplayHit> hits;
      for (const auto& hit : candidate["hits"]) {
        if (!hit.is_object() || !hit.contains("index")) {
          throw std::invalid_argument("B0 candidate replay: malformed hit");
        }
        B0ReplayHit spec;
        spec.index = hit["index"].get<int>();
        if (hit.contains("cellID") && !hit["cellID"].is_null()) {
          spec.cellID = hit["cellID"].get<std::uint64_t>();
        }
        hits.push_back(spec);
      }
      rec.candidates.push_back(std::move(hits));
    }
    out.events.emplace(key, std::move(rec));
  }
  return out;
}

std::vector<std::vector<podio::ObjectID>>
bindReplayCandidates(const edm4eic::TrackerHitCollection& hits,
                     const std::vector<std::vector<B0ReplayHit>>& candidates) {
  std::vector<std::vector<podio::ObjectID>> out;
  out.reserve(candidates.size());
  for (const auto& candidate : candidates) {
    std::vector<podio::ObjectID> ids;
    ids.reserve(candidate.size());
    for (const auto& spec : candidate) {
      if (spec.index < 0 || static_cast<std::size_t>(spec.index) >= hits.size()) {
        throw std::invalid_argument("B0 candidate replay: hit index " + std::to_string(spec.index) +
                                    " is outside RecHits size " + std::to_string(hits.size()));
      }
      const auto hit = hits[static_cast<std::size_t>(spec.index)];
      const auto oid = hit.getObjectID();
      if (oid.index != spec.index) {
        throw std::invalid_argument(
            "B0 candidate replay: ObjectID.index " + std::to_string(oid.index) +
            " disagrees with collection order " + std::to_string(spec.index));
      }
      if (spec.cellID.has_value() && static_cast<std::uint64_t>(hit.getCellID()) != *spec.cellID) {
        throw std::invalid_argument("B0 candidate replay: cellID mismatch at index " +
                                    std::to_string(spec.index));
      }
      ids.push_back(oid);
    }
    out.push_back(std::move(ids));
  }
  return out;
}

} // namespace eicrecon
