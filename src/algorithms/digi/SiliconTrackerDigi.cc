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
#include <memory>
#include <random>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "SiliconTrackerDigi.h"
#include "algorithms/digi/SiliconTrackerDigiConfig.h"

namespace eicrecon {

namespace {
struct CellAccumulator {
  double energyDeposit{0.0};
  std::int32_t earliestTimeStamp{0};
  bool hasTimeStamp{false};
};
} // namespace

void SiliconTrackerDigi::init() {}

void SiliconTrackerDigi::process(const SiliconTrackerDigi::Input& input,
                                 const SiliconTrackerDigi::Output& output) const {

  const auto [headers, sim_hits]       = input;
  auto [raw_hits, links, associations] = output;

  // local random generator
  auto seed = m_uid.getUniqueID(*headers, name());
  std::default_random_engine generator(seed);
  std::normal_distribution<double> gaussian;

  // Accumulate the complete signal in each readout cell before applying the
  // electronics threshold. Thresholding individual simulation contributions
  // loses a channel whose summed signal is above threshold and can make truth
  // associations inconsistent with the emitted raw charge.
  std::unordered_map<std::uint64_t, CellAccumulator> cell_hit_map;

  for (const auto& sim_hit : *sim_hits) {

    // time smearing
    double time_smearing = gaussian(generator) * m_cfg.timeResolution;
    double result_time   = sim_hit.getTime() + time_smearing;
    auto hit_time_stamp  = (std::int32_t)(result_time * 1e3);

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
    debug("   time smearing: {:.4f}, resulting time = {:.4f} [ns]", time_smearing, result_time);
    debug("   hit_time_stamp: {} [~ps]", hit_time_stamp);

    auto& cell = cell_hit_map[sim_hit.getCellID()];
    cell.energyDeposit += sim_hit.getEDep();
    if (!cell.hasTimeStamp) {
      cell.earliestTimeStamp = hit_time_stamp;
      cell.hasTimeStamp      = true;
    } else {
      cell.earliestTimeStamp = std::min(cell.earliestTimeStamp, hit_time_stamp);
    }
  }

  for (const auto& [cell_id, cell] : cell_hit_map) {
    if (cell.energyDeposit < m_cfg.threshold) {
      debug("Cell {} summed edep {:.2f} is below threshold of {:.2f} [keV]", cell_id,
            cell.energyDeposit, m_cfg.threshold / dd4hep::keV);
      continue;
    }

    edm4eic::MutableRawTrackerHit raw_hit_value{
        cell_id, (std::int32_t)std::llround(cell.energyDeposit * 1e6), cell.earliestTimeStamp};
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
