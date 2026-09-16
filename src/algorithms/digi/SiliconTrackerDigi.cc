// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Wouter Deconinck, Sylvester Joosten, Dmitry Romanov

#include <Evaluator/DD4hepUnits.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/Vector3d.h>
#include <edm4hep/Vector3f.h>
#include <podio/detail/Link.h>
#include <podio/detail/LinkCollectionImpl.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "SiliconTrackerDigi.h"
#include "algorithms/digi/SiliconTrackerDigiConfig.h"

namespace eicrecon {

namespace {
struct CellAccumulator {
  double energyDeposit{0.0};
  double earliestTrueTime{0.0};
  bool hasTime{false};
};

std::int32_t channelTimeStamp(std::uint64_t event_seed, std::uint64_t cell_id,
                              double true_time, double time_resolution) {
  // Give every channel its own deterministic random stream. This keeps the
  // reported time independent of SimTrackerHit insertion order and of how many
  // other cells happen to be present in the event.
  std::seed_seq seed_sequence{static_cast<std::uint32_t>(event_seed),
                              static_cast<std::uint32_t>(event_seed >> 32),
                              static_cast<std::uint32_t>(cell_id),
                              static_cast<std::uint32_t>(cell_id >> 32)};
  std::default_random_engine generator(seed_sequence);
  std::normal_distribution<double> gaussian;
  const double measured_time = true_time + gaussian(generator) * time_resolution;
  const double timestamp_ps  = measured_time * 1.0e3;
  if (!std::isfinite(timestamp_ps) ||
      timestamp_ps < static_cast<double>(std::numeric_limits<std::int32_t>::min()) ||
      timestamp_ps > static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
    throw std::overflow_error("SiliconTrackerDigi channel timestamp is outside int32 ps range");
  }
  return static_cast<std::int32_t>(std::llround(timestamp_ps));
}
} // namespace

void SiliconTrackerDigi::init() {
  if (!std::isfinite(m_cfg.threshold) || m_cfg.threshold < 0.0) {
    throw std::invalid_argument("SiliconTrackerDigi threshold must be finite and non-negative");
  }
  if (!std::isfinite(m_cfg.timeResolution) || m_cfg.timeResolution < 0.0) {
    throw std::invalid_argument("SiliconTrackerDigi timeResolution must be finite and non-negative");
  }
}

void SiliconTrackerDigi::process(const SiliconTrackerDigi::Input& input,
                                 const SiliconTrackerDigi::Output& output) const {

  const auto [headers, sim_hits]       = input;
  auto [raw_hits, links, associations] = output;

  const auto event_seed = static_cast<std::uint64_t>(m_uid.getUniqueID(*headers, name()));

  // Accumulate the complete signal and the deterministic physical first-arrival
  // time in each readout cell before applying electronics effects. The generic
  // timing model is deliberately simple: earliest unsmeared contribution plus
  // one Gaussian channel-level measurement error. Detector-specific time walk,
  // pulse shaping, threshold crossing and common-clock effects belong in a
  // detector response model (for example the opt-in B0 AC-LGAD chain).
  std::unordered_map<std::uint64_t, CellAccumulator> cell_hit_map;

  for (const auto& sim_hit : *sim_hits) {
    debug("--------------------");
    debug("Hit cellID   = {}", sim_hit.getCellID());
    debug("   position  = ({:.2f}, {:.2f}, {:.2f})", sim_hit.getPosition().x,
          sim_hit.getPosition().y, sim_hit.getPosition().z);
    debug("   xy_radius = {:.2f}", std::hypot(sim_hit.getPosition().x, sim_hit.getPosition().y));
    debug("   momentum  = ({:.2f}, {:.2f}, {:.2f})", sim_hit.getMomentum().x,
          sim_hit.getMomentum().y, sim_hit.getMomentum().z);
    debug("   edep = {:.2f}", sim_hit.getEDep());
    debug("   time = {:.4f}[ns]", sim_hit.getTime());
    debug("   particle time = {}[ns]", sim_hit.getParticle().getTime());

    if (!std::isfinite(sim_hit.getTime())) {
      throw std::invalid_argument("SiliconTrackerDigi received a non-finite SimTrackerHit time");
    }

    auto& cell = cell_hit_map[sim_hit.getCellID()];
    cell.energyDeposit += sim_hit.getEDep();
    if (!cell.hasTime) {
      cell.earliestTrueTime = sim_hit.getTime();
      cell.hasTime          = true;
    } else {
      cell.earliestTrueTime = std::min(cell.earliestTrueTime, static_cast<double>(sim_hit.getTime()));
    }
  }

  for (const auto& [cell_id, cell] : cell_hit_map) {
    if (cell.energyDeposit < m_cfg.threshold) {
      debug("Cell {} summed edep {:.2f} is below threshold of {:.2f} [keV]", cell_id,
            cell.energyDeposit, m_cfg.threshold / dd4hep::keV);
      continue;
    }
    if (!cell.hasTime) {
      throw std::logic_error("SiliconTrackerDigi accumulated a signal without a channel time");
    }

    const auto hit_time_stamp =
        channelTimeStamp(event_seed, cell_id, cell.earliestTrueTime, m_cfg.timeResolution);
    debug("Cell {} earliest true time {:.4f} ns -> channel timestamp {} ps", cell_id,
          cell.earliestTrueTime, hit_time_stamp);

    edm4eic::MutableRawTrackerHit raw_hit_value{
        cell_id, static_cast<std::int32_t>(std::llround(cell.energyDeposit * 1e6)), hit_time_stamp};
    raw_hits->push_back(raw_hit_value);
    auto raw_hit = raw_hits->at(raw_hits->size() - 1);

    for (const auto& sim_hit : *sim_hits) {
      if (cell_id != sim_hit.getCellID()) {
        continue;
      }

      // Energy-fraction weights make truth relations describe the signal that
      // actually contributed to the emitted channel hit. A zero-energy channel
      // can only pass when the threshold is zero; retain its relations with
      // zero weight rather than divide by zero.
      const double weight =
          cell.energyDeposit > 0.0 ? sim_hit.getEDep() / cell.energyDeposit : 0.0;

      auto link = links->create();
      link.setFrom(raw_hit);
      link.setTo(sim_hit);
      link.setWeight(weight);

      auto hitassoc = associations->create();
      hitassoc.setWeight(weight);
      hitassoc.setRawHit(raw_hit);
      hitassoc.setSimHit(sim_hit);
    }
  }
}

} // namespace eicrecon
