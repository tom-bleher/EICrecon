// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "B0ACLGADResponse.h"
#include <charconv>
#include <istream>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>

namespace eicrecon::b0::aclgad {
struct Calibration { Config config; std::string status, reference; };
inline std::string trim(std::string s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  return first == std::string::npos ? "" : s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
}
inline Calibration readCalibration(std::istream& input, bool allowExperimental) {
  std::map<std::string, std::string> values;
  std::string line;
  std::size_t bytes = 0;
  while (std::getline(input, line)) {
    bytes += line.size();
    if (bytes > 65536) throw std::invalid_argument("B0 calibration is too large");
    line = trim(line.substr(0, line.find('#')));
    if (line.empty()) continue;
    const auto separator = line.find('=');
    if (separator == std::string::npos || !values.emplace(trim(line.substr(0, separator)), trim(line.substr(separator + 1))).second)
      throw std::invalid_argument("Malformed or duplicate B0 calibration key");
  }
  if (input.bad()) throw std::runtime_error("Cannot read B0 calibration");
  auto take = [&](const char* name) {
    const auto found = values.find(name);
    if (found == values.end() || found->second.empty()) throw std::invalid_argument(std::string("Missing B0 calibration key: ") + name);
    auto value = found->second; values.erase(found); return value;
  };
  if (take("schemaVersion") != "1") throw std::invalid_argument("Unsupported B0 response schema");
  Calibration result; result.status = take("status"); result.reference = take("reference");
  if (result.status != "calibrated" && result.status != "experimental") throw std::invalid_argument("Invalid B0 calibration status");
  if (result.status == "experimental" && !allowExperimental)
    throw std::invalid_argument("Experimental B0 response requires explicit allowExperimental=true");
  auto number = [&](const char* name, auto& output) {
    const std::string text = take(name);
    using T = std::decay_t<decltype(output)>;
    if constexpr (std::is_integral_v<T>) {
      auto parsed = std::from_chars(text.data(), text.data() + text.size(), output);
      if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        throw std::invalid_argument(std::string("Invalid B0 integer: ") + name);
    } else {
      std::istringstream stream(text); stream.imbue(std::locale::classic());
      stream >> output;
      if (!stream || stream.peek() != std::char_traits<char>::eof() || !std::isfinite(output))
        throw std::invalid_argument(std::string("Invalid B0 number: ") + name);
    }
  };
  auto boolean = [&](const char* name, bool& output) {
    const auto value = take(name);
    if (value != "true" && value != "false") throw std::invalid_argument(std::string("Invalid B0 boolean: ") + name);
    output = value == "true";
  };
  auto& c = result.config;
#define B0_NUMBER(member) number(#member, c.member)
  B0_NUMBER(pitchX); B0_NUMBER(pitchY); B0_NUMBER(offsetX); B0_NUMBER(offsetY);
  B0_NUMBER(sigmaSharingX); B0_NUMBER(sigmaSharingY); B0_NUMBER(pairEnergy); B0_NUMBER(gain);
  B0_NUMBER(channelGainSigma); B0_NUMBER(noiseSigma); B0_NUMBER(threshold); B0_NUMBER(adcStep);
  B0_NUMBER(maxADC); B0_NUMBER(gateStart); B0_NUMBER(gateEnd);
  B0_NUMBER(sensorTimeSigma); B0_NUMBER(channelTimeSigma); B0_NUMBER(clockTimeSigma);
  B0_NUMBER(timeWalk); B0_NUMBER(timeStep); B0_NUMBER(clusterTimeWindow);
  B0_NUMBER(modelSigmaX); B0_NUMBER(modelSigmaY); B0_NUMBER(modelCorrelation);
  B0_NUMBER(singleChannelSigmaX); B0_NUMBER(singleChannelSigmaY); B0_NUMBER(edgeErrorScale);
  B0_NUMBER(calibrationSeed); B0_NUMBER(maxDeposits); B0_NUMBER(maxChannels); B0_NUMBER(maxWork);
#undef B0_NUMBER
  boolean("noiseOnlySensors", c.noiseOnlySensors);
  boolean("rejectSaturatedClusters", c.rejectSaturatedClusters);
  if (!values.empty()) throw std::invalid_argument("Unknown B0 calibration key: " + values.begin()->first);
  validate(c);
  return result;
}
} // namespace eicrecon::b0::aclgad
