// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace eicrecon::b0::aclgad {

// All lengths are mm, times ns, deposited energies eV, and amplitudes equivalent
// output electrons. This is a phenomenological pad-response model, not a
// microscopic charge-transport or waveform simulation.
struct Config {
  double pitchX = 0.5, pitchY = 0.5;
  double offsetX = 0.25, offsetY = 0.25;
  double sigmaSharingX = 0.2, sigmaSharingY = 0.2;
  double pairEnergy = 3.6, gain = 20.0;
  double channelGainSigma = 0.02, noiseSigma = 500.0;
  double threshold = 2500.0, adcStep = 100.0;
  std::int32_t maxADC = 32767;
  double gateStart = 0.0, gateEnd = 100.0;
  double sensorTimeSigma = 0.025, channelTimeSigma = 0.015, clockTimeSigma = 0.01;
  double timeWalk = 0.0; // ns*electron; reconstruction corrects with measured amplitude
  double timeStep = 0.001, clusterTimeWindow = 1.0;
  double modelSigmaX = 0.05, modelSigmaY = 0.05, modelCorrelation = 0.0;
  double singleChannelSigmaX = 0.15, singleChannelSigmaY = 0.15;
  double edgeErrorScale = 2.0;
  bool noiseOnlySensors = true;
  bool rejectSaturatedClusters = true;
  std::uint64_t calibrationSeed = 1;
  std::size_t maxDeposits = 100000, maxChannels = 1000000, maxWork = 20000000;
};
struct Sensor { std::uint64_t id; double halfX, halfY; };
struct Deposit {
  std::uint64_t id, sensor; // unique simulation-hit identity and unique sensor volume ID
  double x, y, energy, time;
};
struct Channel {
  std::uint64_t sensor;
  int ix, iy;
  double x, y, charge, signal, time, independentTimeVariance, chargeVariance;
  std::int32_t adc;
  bool saturated;
  std::vector<std::pair<std::uint64_t, double>> truth; // pre-noise signal fractions
};
struct Cluster {
  std::uint64_t sensor;
  std::vector<std::size_t> channels;
  double x = 0, y = 0, xx = 0, yy = 0, xy = 0;
  double time = 0, timeVariance = 0, xt = 0, yt = 0, charge = 0;
  bool edge = false, saturated = false;
};
struct Result {
  std::vector<Channel> channels;
  std::vector<Cluster> clusters;
  std::size_t work = 0, outsideGate = 0, zeroEnergy = 0, belowThreshold = 0;
  std::size_t saturatedChannels = 0, rejectedSaturatedClusters = 0;
  double inputSignal = 0, sharedSignal = 0; // before channel gain/noise; difference is edge loss
};

inline std::uint64_t mix(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
inline double uniform(std::uint64_t key) {
  return (static_cast<double>(mix(key) >> 12) + 0.5) * 0x1.0p-52;
}
inline double normal(std::uint64_t key) {
  // Open-interval uniform avoids log(0); no mutable RNG or thread-local state.
  return std::sqrt(-2.0 * std::log(uniform(key))) *
         std::cos(6.2831853071795864769 * uniform(key ^ 0xa0761d6478bd642fULL));
}
inline void validate(const Config& c) {
  for (double v : {c.pitchX, c.pitchY, c.sigmaSharingX, c.sigmaSharingY, c.pairEnergy,
                   c.gain, c.adcStep, c.timeStep, c.clusterTimeWindow,
                   c.modelSigmaX, c.modelSigmaY, c.singleChannelSigmaX, c.singleChannelSigmaY})
    if (!std::isfinite(v) || !(v > 0)) throw std::invalid_argument("B0 AC-LGAD: positive finite parameter required");
  for (double v : {c.channelGainSigma, c.noiseSigma, c.threshold, c.sensorTimeSigma,
                   c.channelTimeSigma, c.clockTimeSigma, c.timeWalk})
    if (!std::isfinite(v) || v < 0) throw std::invalid_argument("B0 AC-LGAD: negative/non-finite noise or threshold");
  if (!std::isfinite(c.offsetX) || !std::isfinite(c.offsetY) ||
      !std::isfinite(c.gateStart) || !std::isfinite(c.gateEnd) || c.gateEnd <= c.gateStart ||
      !std::isfinite(c.modelCorrelation) || std::abs(c.modelCorrelation) >= 1 ||
      !std::isfinite(c.edgeErrorScale) || c.edgeErrorScale < 1 || c.channelGainSigma > 1 ||
      c.maxADC <= 0 || c.maxDeposits == 0 || c.maxChannels == 0 || c.maxWork == 0 ||
      !std::isfinite(c.maxADC * c.adcStep))
    throw std::invalid_argument("B0 AC-LGAD: invalid configuration");
}
inline std::pair<int, int> channelRange(double half, double pitch, double offset) {
  if (!std::isfinite(half) || half <= 0 || !std::isfinite(pitch) || pitch <= 0 || !std::isfinite(offset)) throw std::invalid_argument("Invalid B0 sensor extent");
  // Require complete cells: the profile offsets align physical channel edges to
  // sensor edges. A half cell is not silently added beyond the sensor.
  const double low = (-half + pitch / 2 - offset) / pitch;
  const double high = (half - pitch / 2 - offset) / pitch;
  if (low < -32768 || high > 32767 || high < low ||
      std::abs(low - std::round(low)) > 1e-7 || std::abs(high - std::round(high)) > 1e-7)
    throw std::invalid_argument("B0 AC-LGAD: sensor edges do not align with physical cells");
  return {static_cast<int>(std::round(low)), static_cast<int>(std::round(high))};
}
inline double fraction(double coordinate, double center, double pitch, double sigma) {
  const double s = std::sqrt(2.0) * sigma;
  return std::max(0.0, 0.5 * (std::erf((center + pitch / 2 - coordinate) / s) -
                              std::erf((center - pitch / 2 - coordinate) / s)));
}

inline Result simulate(const std::vector<Sensor>& sensors, std::vector<Deposit> deposits,
                       const Config& c, std::uint64_t eventSeed) {
  validate(c);
  if (deposits.size() > c.maxDeposits) throw std::runtime_error("B0 AC-LGAD deposit budget exceeded");
  Result out;
  auto spend = [&] {
    if (out.work == c.maxWork) throw std::runtime_error("B0 AC-LGAD work budget exceeded; event rejected, not truncated");
    ++out.work;
  };
  std::map<std::uint64_t, Sensor> sensorMap;
  for (const auto& s : sensors) {
    channelRange(s.halfX, c.pitchX, c.offsetX); channelRange(s.halfY, c.pitchY, c.offsetY);
    if (!sensorMap.emplace(s.id, s).second) throw std::invalid_argument("Duplicate B0 sensor identity");
  }
  std::sort(deposits.begin(), deposits.end(), [](const Deposit& a, const Deposit& b) { return a.id < b.id; });
  std::map<std::uint64_t, std::vector<Deposit>> grouped;
  for (std::size_t i = 0; i < deposits.size(); ++i) {
    const auto& d = deposits[i];
    if (i && d.id == deposits[i - 1].id) throw std::invalid_argument("Duplicate B0 simulation-hit identity");
    const auto found = sensorMap.find(d.sensor);
    if (found == sensorMap.end() || !std::isfinite(d.x) || !std::isfinite(d.y) ||
        !std::isfinite(d.energy) || d.energy < 0 || !std::isfinite(d.time))
      throw std::invalid_argument("Invalid B0 deposit or missing sensor");
    if (d.x < -found->second.halfX || d.x > found->second.halfX ||
        d.y < -found->second.halfY || d.y > found->second.halfY)
      throw std::invalid_argument("B0 deposit outside active sensor");
    if (d.time < c.gateStart || d.time >= c.gateEnd) { ++out.outsideGate; continue; }
    if (d.energy == 0) { ++out.zeroEnergy; continue; }
    grouped[d.sensor].push_back(d);
  }
  const double clockShift = c.clockTimeSigma * normal(eventSeed ^ 0xe7037ed1a0b428dbULL);
  using Grid = std::map<std::pair<int, int>, std::size_t>;
  std::map<std::uint64_t, Grid> grids;
  std::map<std::uint64_t, std::tuple<int, int, int, int>> ranges;
  for (const auto& [sensorId, sensor] : sensorMap) {
    const auto found = grouped.find(sensorId);
    if (found == grouped.end() && !c.noiseOnlySensors) continue;
    const auto xr = channelRange(sensor.halfX, c.pitchX, c.offsetX);
    const auto yr = channelRange(sensor.halfY, c.pitchY, c.offsetY);
    ranges.emplace(sensorId, std::make_tuple(xr.first, xr.second, yr.first, yr.second));
    struct Signal { double q = 0, qt = 0; std::vector<std::pair<std::uint64_t, double>> truth; };
    std::map<std::pair<int, int>, Signal> signal;
    if (found != grouped.end()) for (const auto& deposit : found->second) {
      const double q = deposit.energy / c.pairEnergy * c.gain;
      if (!std::isfinite(q)) throw std::overflow_error("B0 signal overflow");
      out.inputSignal += q;
      const double time = deposit.time + c.sensorTimeSigma * normal(eventSeed ^ mix(deposit.id) ^ 0x8ebc6af09c88c6e3ULL);
      for (int ix = xr.first; ix <= xr.second; ++ix) {
        const double wx = fraction(deposit.x, ix * c.pitchX + c.offsetX, c.pitchX, c.sigmaSharingX);
        for (int iy = yr.first; iy <= yr.second; ++iy) {
          spend();
          const double part = q * wx * fraction(deposit.y, iy * c.pitchY + c.offsetY, c.pitchY, c.sigmaSharingY);
          if (part <= 0) continue;
          auto& s = signal[{ix, iy}]; s.q += part; s.qt += part * time;
          s.truth.emplace_back(deposit.id, part); out.sharedSignal += part;
        }
      }
    }
    for (int ix = xr.first; ix <= xr.second; ++ix) for (int iy = yr.first; iy <= yr.second; ++iy) {
      spend();
      const auto key = mix(sensorId) ^ mix(static_cast<std::uint32_t>(ix)) ^
                       mix(static_cast<std::uint32_t>(iy) + 0xd1b54a32d192ed03ULL);
      const auto entry = signal.find({ix, iy});
      const double q = entry == signal.end() ? 0 : entry->second.q;
      const double gain = std::exp(c.channelGainSigma * normal(c.calibrationSeed ^ key) -
                                   0.5 * c.channelGainSigma * c.channelGainSigma);
      const double measured = q * gain + c.noiseSigma * normal(eventSeed ^ key ^ 0x589965cc75374cc3ULL);
      if (!std::isfinite(measured)) throw std::overflow_error("B0 channel amplitude overflow");
      const bool saturated = measured >= c.maxADC * c.adcStep;
      const auto adc = static_cast<std::int32_t>(std::llround(std::clamp(measured / c.adcStep, 0.0, double(c.maxADC))));
      const double charge = adc * c.adcStep;
      // Threshold is applied ONCE, to the accumulated, noisy digitized channel.
      if (charge <= 0 || charge < c.threshold) { ++out.belowThreshold; continue; }
      if (out.channels.size() == c.maxChannels) throw std::runtime_error("B0 AC-LGAD output channel budget exceeded");
      const double baseTime = q > 0 ? entry->second.qt / q :
          c.gateStart + (c.gateEnd - c.gateStart) * uniform(eventSeed ^ key ^ 0x1d8e4e27c47d124fULL);
      double time = baseTime + clockShift + c.channelTimeSigma * normal(eventSeed ^ key ^ 0xeb44accab455d165ULL);
      if (q > 0) time += c.timeWalk / std::max(q * gain, c.adcStep);
      time = std::round(time / c.timeStep) * c.timeStep;
      if (!std::isfinite(time)) throw std::overflow_error("B0 channel time overflow");
      Channel channel{sensorId, ix, iy, ix * c.pitchX + c.offsetX, iy * c.pitchY + c.offsetY,
          charge, q * gain, time, c.channelTimeSigma * c.channelTimeSigma + c.timeStep * c.timeStep / 12,
          c.noiseSigma * c.noiseSigma + c.adcStep * c.adcStep / 12, adc, saturated, {}};
      // Positive electronics noise is an unassociated fraction, never credited
      // to a vanishing response tail. Negative noise cannot create weight > 1.
      if (q > 0) for (const auto& [id, part] : entry->second.truth)
        channel.truth.emplace_back(id, part * gain / std::max(q * gain, charge));
      grids[sensorId][{ix, iy}] = out.channels.size();
      out.channels.push_back(std::move(channel));
      if (saturated) ++out.saturatedChannels;
    }
  }
  std::set<std::size_t> visited;
  for (const auto& [sensor, grid] : grids) {
    const auto [xmin, xmax, ymin, ymax] = ranges.at(sensor);
    for (const auto& [coordinate, index] : grid) {
      (void)coordinate;
      if (!visited.insert(index).second) continue;
      Cluster cluster; cluster.sensor = sensor; cluster.channels.push_back(index);
      double tmin = out.channels[index].time, tmax = tmin;
      for (std::size_t cursor = 0; cursor < cluster.channels.size(); ++cursor) {
        const auto& ch = out.channels[cluster.channels[cursor]];
        for (const auto& [dx, dy] : std::array<std::pair<int, int>, 4>{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}}) {
          spend();
          auto next = grid.find({ch.ix + dx, ch.iy + dy});
          if (next == grid.end() || visited.count(next->second)) continue;
          const double time = out.channels[next->second].time;
          if (std::max(tmax, time) - std::min(tmin, time) > c.clusterTimeWindow) continue;
          visited.insert(next->second); cluster.channels.push_back(next->second);
          tmin = std::min(tmin, time); tmax = std::max(tmax, time);
        }
      }
      std::sort(cluster.channels.begin(), cluster.channels.end());
      for (auto i : cluster.channels) {
        const auto& ch = out.channels[i];
        cluster.charge += ch.charge; cluster.x += ch.charge * ch.x; cluster.y += ch.charge * ch.y;
        cluster.time += ch.charge * (ch.time - c.timeWalk / ch.charge);
        cluster.edge |= ch.ix == xmin || ch.ix == xmax || ch.iy == ymin || ch.iy == ymax;
        cluster.saturated |= ch.saturated;
      }
      if (cluster.saturated && c.rejectSaturatedClusters) { ++out.rejectedSaturatedClusters; continue; }
      cluster.x /= cluster.charge; cluster.y /= cluster.charge; cluster.time /= cluster.charge;
      const bool single = cluster.channels.size() == 1;
      const double scale = cluster.edge ? c.edgeErrorScale : 1;
      const double sx = scale * (single ? c.singleChannelSigmaX : c.modelSigmaX);
      const double sy = scale * (single ? c.singleChannelSigmaY : c.modelSigmaY);
      cluster.xx = sx * sx; cluster.yy = sy * sy; cluster.xy = c.modelCorrelation * sx * sy;
      // Sensor jitter is shared by channels from a crossing. Do not falsely
      // improve it (or the event clock) by sqrt(number of pads).
      cluster.timeVariance = c.sensorTimeSigma * c.sensorTimeSigma + c.clockTimeSigma * c.clockTimeSigma;
      for (auto i : cluster.channels) {
        const auto& ch = out.channels[i]; const double w = ch.charge / cluster.charge;
        const double dx = (ch.x - cluster.x) / cluster.charge;
        const double dy = (ch.y - cluster.y) / cluster.charge;
        // For sum(q_i*t_i - timeWalk)/sum(q_i), d t / d q_i = (t_i-t)/sum(q_i).
        const double dt = (ch.time - cluster.time) / cluster.charge;
        cluster.xx += dx * dx * ch.chargeVariance; cluster.yy += dy * dy * ch.chargeVariance;
        cluster.xy += dx * dy * ch.chargeVariance;
        cluster.xt += dx * dt * ch.chargeVariance; cluster.yt += dy * dt * ch.chargeVariance;
        cluster.timeVariance += w * w * ch.independentTimeVariance + dt * dt * ch.chargeVariance;
      }
      const double determinant2 = cluster.xx * cluster.yy - cluster.xy * cluster.xy;
      const double determinant3 = determinant2 * cluster.timeVariance -
          cluster.yy * cluster.xt * cluster.xt - cluster.xx * cluster.yt * cluster.yt +
          2 * cluster.xy * cluster.xt * cluster.yt;
      if (!std::isfinite(cluster.x + cluster.y + cluster.time + cluster.charge) ||
          !std::isfinite(cluster.xx + cluster.yy + cluster.xy + cluster.timeVariance + cluster.xt + cluster.yt) ||
          !(determinant2 > 0) || !(determinant3 > 0))
        throw std::runtime_error("Invalid B0 reconstructed covariance");
      out.clusters.push_back(std::move(cluster));
    }
  }
  return out;
}
} // namespace eicrecon::b0::aclgad
