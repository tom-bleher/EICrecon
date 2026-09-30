// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "SiliconPadClustering.h"

#include <DD4hep/Detector.h>
#include <DD4hep/Objects.h>
#include <DD4hep/Readout.h>
#include <DD4hep/Segmentations.h>
#include <DD4hep/Shapes.h>
#include <DD4hep/VolumeManager.h>
#include <Evaluator/DD4hepUnits.h>
#include <algorithms/geo.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4hep/Vector3f.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/digi/TrapezoidPadResponse.h"

namespace eicrecon {

void SiliconPadClustering::init() {
  m_converter    = algorithms::GeoSvc::instance().cellIDPositionConverter();
  const auto seg = algorithms::GeoSvc::instance().detector()->readout(m_cfg.readout).segmentation();
  m_grid         = dynamic_cast<const dd4hep::DDSegmentation::CartesianGridXY*>(seg.segmentation());
  if (m_grid == nullptr) {
    throw std::runtime_error("SiliconPadClustering needs a CartesianGridXY readout, not " +
                             seg.type() + " for " + m_cfg.readout);
  }
  if (m_cfg.electrode_x >= m_grid->gridSizeX() || m_cfg.electrode_y >= m_grid->gridSizeY()) {
    warning("Electrodes are at least as wide as the {} pads: positions are pad centres",
            m_cfg.readout);
  }
}

void SiliconPadClustering::process(const Input& input, const Output& output) const {
  using dd4hep::mm;
  const auto [pads] = input;
  auto [clusters]   = output;

  const auto* decoder = m_grid->decoder();
  const auto& fx      = m_grid->fieldNameX();
  const auto& fy      = m_grid->fieldNameY();

  // The pad and its eight neighbours
  auto around = [&](dd4hep::CellID cellID) {
    const auto ix = decoder->get(cellID, fx);
    const auto iy = decoder->get(cellID, fy);
    std::vector<dd4hep::CellID> cells;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        dd4hep::CellID neighbour = cellID;
        decoder->set(neighbour, fx, ix + dx);
        decoder->set(neighbour, fy, iy + dy);
        cells.push_back(neighbour);
      }
    }
    return cells;
  };

  // With a seed (discriminator) threshold, a pad below it is only read out next to a fired pad
  std::unordered_set<dd4hep::CellID> fired;
  for (const auto& pad : *pads) {
    if (pad.getEdep() >= m_cfg.seed_threshold) {
      fired.insert(pad.getCellID());
    }
  }
  std::vector<bool> read(pads->size(), true);
  if (m_cfg.seed_threshold > 0) {
    for (std::size_t i = 0; i < pads->size(); ++i) {
      const auto cells = around((*pads)[i].getCellID());
      read[i]          = std::ranges::any_of(cells, [&fired](auto c) { return fired.contains(c); });
    }
  }

  // Union-find over read pads that touch and are compatible in time
  std::vector<std::size_t> parent(pads->size());
  std::iota(parent.begin(), parent.end(), std::size_t{0});
  auto root = [&parent](std::size_t i) {
    while (parent[i] != i) {
      i = parent[i] = parent[parent[i]];
    }
    return i;
  };
  std::unordered_map<dd4hep::CellID, std::vector<std::size_t>> pads_in_cell;
  for (std::size_t i = 0; i < pads->size(); ++i) {
    if (read[i]) {
      pads_in_cell[(*pads)[i].getCellID()].push_back(i);
    }
  }
  for (std::size_t i = 0; i < pads->size(); ++i) {
    if (!read[i]) {
      continue;
    }
    const auto& pad = (*pads)[i];
    for (const auto neighbour : around(pad.getCellID())) {
      const auto it = pads_in_cell.find(neighbour);
      if (it == pads_in_cell.end()) {
        continue;
      }
      for (const auto j : it->second) {
        if (std::abs(pad.getTime() - (*pads)[j].getTime()) < m_cfg.deltaT) {
          parent[root(i)] = root(j);
        }
      }
    }
  }
  std::vector<std::vector<std::size_t>> groups;
  std::unordered_map<std::size_t, std::size_t> group_of_root;
  for (std::size_t i = 0; i < pads->size(); ++i) {
    if (!read[i]) {
      continue;
    }
    const auto [it, inserted] = group_of_root.try_emplace(root(i), groups.size());
    if (inserted) {
      groups.emplace_back();
    }
    groups[it->second].push_back(i);
  }

  // Position and resolution along one axis from the energy summed per column of pads: the
  // leading column and its larger adjacent column share the deposit as in trapezoidPadShare.
  // A cluster wider than two columns (delta ray, overlapping particle) gets pitch/sqrt(12), as
  // does a single pad at the sensor edge (half_width), whose outer neighbour does not exist
  auto axis = [](const std::map<long, double>& columns, long seed_index, double seed_centre,
                 double pitch, double electrode, double single_resolution, double shared_resolution,
                 double half_width) {
    constexpr double sqrt_12 = 3.4641016151;
    const auto lead =
        std::ranges::max_element(columns, {}, &std::map<long, double>::value_type::second);
    const double centre   = seed_centre + static_cast<double>(lead->first - seed_index) * pitch;
    double neighbour_edep = 0;
    long direction        = 0;
    for (const long d : {-1L, 1L}) {
      const auto it = columns.find(lead->first + d);
      if (it != columns.end() && it->second > neighbour_edep) {
        neighbour_edep = it->second;
        direction      = d;
      }
    }
    if (direction == 0) {
      const bool edge = std::abs(centre) + pitch > half_width;
      return std::pair{centre, single_resolution > 0 && columns.size() == 1 && !edge
                                   ? single_resolution
                                   : pitch / sqrt_12};
    }
    const double f = neighbour_edep / (lead->second + neighbour_edep);
    return std::pair{
        centre + static_cast<double>(direction) * trapezoidSharedOffset(f, pitch, electrode),
        shared_resolution > 0 && columns.size() == 2 ? shared_resolution : pitch / sqrt_12};
  };

  for (const auto& group : groups) {
    const std::size_t seed = *std::ranges::max_element(
        group, {}, [&pads](std::size_t i) { return (*pads)[i].getEdep(); });
    const auto& seed_pad = (*pads)[seed];
    if (seed_pad.getEdep() < m_cfg.seed_threshold) {
      continue;
    }
    const auto seed_id  = seed_pad.getCellID();
    const auto* context = m_converter->findContext(seed_id);
    if (context == nullptr) {
      error("No volume context for cell ID {:x}", seed_id);
      continue;
    }

    std::map<long, double> columns;
    std::map<long, double> rows;
    float edep = 0;
    for (const auto i : group) {
      const auto& pad = (*pads)[i];
      columns[decoder->get(pad.getCellID(), fx)] += pad.getEdep();
      rows[decoder->get(pad.getCellID(), fy)] += pad.getEdep();
      edep += pad.getEdep();
    }

    const auto seed_local = m_grid->position(seed_id);
    // Sensor half-widths (no edge treatment when the sensor is not a box)
    double half_x = std::numeric_limits<double>::max();
    double half_y = std::numeric_limits<double>::max();
    try {
      const dd4hep::Box box = context->element.solid();
      half_x                = box.x();
      half_y                = box.y();
    } catch (const std::exception&) {
    }
    const auto [u, sigma_u] =
        axis(columns, decoder->get(seed_id, fx), seed_local.x(), m_grid->gridSizeX(),
             m_cfg.electrode_x, m_cfg.single_pad_resolution, m_cfg.shared_pad_resolution, half_x);
    const auto [v, sigma_v] =
        axis(rows, decoder->get(seed_id, fy), seed_local.y(), m_grid->gridSizeY(),
             m_cfg.electrode_y, m_cfg.single_pad_resolution, m_cfg.shared_pad_resolution, half_y);
    const auto global = context->localToWorld(dd4hep::Position(u, v, seed_local.z()));

    // Covariance in local segmentation coordinates, as in TrackerHitReconstruction
    auto cluster = clusters->create(
        seed_id,
        edm4hep::Vector3f{static_cast<float>(global.x() / mm), static_cast<float>(global.y() / mm),
                          static_cast<float>(global.z() / mm)},
        edm4eic::CovDiag3f{(sigma_u / mm) * (sigma_u / mm), (sigma_v / mm) * (sigma_v / mm), 0.},
        seed_pad.getTime(), seed_pad.getTimeError(), edep, 0.F);
    cluster.setRawHit(seed_pad.getRawHit());
    trace("cluster of {} pads at local ({:.4f}, {:.4f}) mm, sigma ({:.4f}, {:.4f}) mm",
          group.size(), u / mm, v / mm, sigma_u / mm, sigma_v / mm);
  }
}

} // namespace eicrecon
