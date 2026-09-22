// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher
#include "B0TelescopeSeedFinder.h"

#if EICRECON_HAS_B0_TELESCOPE
#include "B0TelescopeMath.h"
#include <Acts/Definitions/Units.hpp>
#include <algorithm>
#include <limits>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace eicrecon::b0telescope {
namespace {

Point position(const SpacePoint& sp, double angle) {
  constexpr double mm = Acts::UnitConstants::mm;
  return toTelescopeFrame({sp.xy()[0] / mm, sp.xy()[1] / mm, sp.zr()[0] / mm}, angle);
}

class TelescopeDoubletFinder final : public Acts::DoubletSeedFinder {
public:
  TelescopeDoubletFinder(bool bottom, double angle, const B0TelescopeSeedingConfig& cfg)
      : m_bottom(bottom), m_angle(angle), m_cfg(cfg),
        m_contract(Config{.spacePointsSortedByRadius = false}, Acts::UnitConstants::T) {}

  const DerivedConfig& config() const override { return m_contract; }

  void createDoublets(const SpacePoint& middle, const Acts::MiddleSpInfo&,
                       DoubletView<SpacePoints::ConstSubset> candidates,
                       Acts::DoubletsForMiddleSp& output) const override {
    build(middle, candidates, output);
  }
  void createDoublets(const SpacePoint& middle, const Acts::MiddleSpInfo&,
                       DoubletView<SpacePoints::ConstRange> candidates,
                       Acts::DoubletsForMiddleSp& output) const override {
    build(middle, candidates, output);
  }

private:
  template <typename Range>
  void build(const SpacePoint& middle, const Range& candidates,
              Acts::DoubletsForMiddleSp& output) const {
    const auto m = position(middle, m_angle);
    for (const auto other : candidates) {
      const auto p = position(other, m_angle);
      if (compatibleDoublet(m_bottom ? p : m, m_bottom ? m : p,
                             m_cfg.minDeltaZ / edm4eic::unit::mm, m_cfg.maxSlope)) {
        // Only indices are consumed by our triplet finder/filter. The solenoidal
        // conformal-map auxiliaries have no telescope meaning and are not used.
        output.emplace_back(other.index(), 0.F, 0.F, 0.F, 0.F, 0.F, 0.F, 0.F);
      }
    }
  }
  bool m_bottom;
  double m_angle;
  const B0TelescopeSeedingConfig& m_cfg;
  DerivedConfig m_contract;
};

class TelescopeTripletFinder final : public Acts::TripletSeedFinder {
public:
  TelescopeTripletFinder(double angle, const B0TelescopeSeedingConfig& cfg)
      : m_angle(angle), m_cfg(cfg),
        m_contract(Config{.sortedByCotTheta = false}, Acts::UnitConstants::T) {}

  const DerivedConfig& config() const override { return m_contract; }

  TripletReturn<Acts::DoubletsForMiddleSp::Range>
  createTripletTopCandidates(const SpacePoints& points, const SpacePoint& middle,
                             const Acts::DoubletsForMiddleSp::Proxy& bottom,
                             TripletView<Acts::DoubletsForMiddleSp::Range> tops,
                             Acts::TripletTopCandidates& output) const override {
    return build<Acts::DoubletsForMiddleSp::Range>(points, middle, bottom, tops, output);
  }
  TripletReturn<Acts::DoubletsForMiddleSp::Subset>
  createTripletTopCandidates(const SpacePoints& points, const SpacePoint& middle,
                             const Acts::DoubletsForMiddleSp::Proxy& bottom,
                             TripletView<Acts::DoubletsForMiddleSp::Subset> tops,
                             Acts::TripletTopCandidates& output) const override {
    return build<Acts::DoubletsForMiddleSp::Subset>(points, middle, bottom, tops, output);
  }
  TripletReturn<Acts::DoubletsForMiddleSp::Subset2>
  createTripletTopCandidates(const SpacePoints& points, const SpacePoint& middle,
                             const Acts::DoubletsForMiddleSp::Proxy& bottom,
                             TripletView<Acts::DoubletsForMiddleSp::Subset2> tops,
                             Acts::TripletTopCandidates& output) const override {
    return build<Acts::DoubletsForMiddleSp::Subset2>(points, middle, bottom, tops, output);
  }

private:
  template <typename Range>
  TripletReturn<Range> build(const SpacePoints& points, const SpacePoint& middle,
                              const Acts::DoubletsForMiddleSp::Proxy& bottom, TripletView<Range> tops,
                              Acts::TripletTopCandidates& output) const {
    const auto a = position(points[bottom.spacePointIndex()], m_angle);
    const auto b = position(middle, m_angle);
    for (const auto top : tops) {
      if (tripletScore(a, b, position(points[top.spacePointIndex()], m_angle),
                        m_cfg.maxResidualX / edm4eic::unit::mm,
                        m_cfg.maxResidualY / edm4eic::unit::mm)) {
        // Do not interpret these as collider curvature/impact parameters.
        // Our filter evaluates the Cartesian score directly from the points.
        output.emplace_back(top.spacePointIndex(), 0.F, 0.F);
      }
    }
    if constexpr (!mutableTripletViews<Acts::TripletSeedFinder>) {
      return tops; // No radial/cot(theta) cursor pruning in the telescope policy.
    }
  }
  double m_angle;
  const B0TelescopeSeedingConfig& m_cfg;
  DerivedConfig m_contract;
};

class TelescopeSeedFilter final : public Acts::ITripletSeedFilter {
public:
  TelescopeSeedFilter(double angle, const B0TelescopeSeedingConfig& cfg)
      : m_angle(angle), m_cfg(cfg) {}

  bool sufficientTopDoublets(const SpacePoints&, const SpacePoint&,
                             const Acts::DoubletsForMiddleSp& tops) const override {
    m_candidates.clear(); // Also reset when ACTS subsequently finds no bottom.
    return !tops.empty();
  }

  void filterTripletTopCandidates(const SpacePoints& points, const SpacePoint& middle,
                                  const Acts::DoubletsForMiddleSp::Proxy& bottom,
                                  const Acts::TripletTopCandidates& tops) const override {
    for (const auto top : tops) {
      const std::array<std::uint32_t, 3> ids{
          bottom.spacePointIndex(), middle.index(), top.spacePoint()};
      const auto score = tripletScore(position(points[ids[0]], m_angle),
                                      position(points[ids[1]], m_angle),
                                      position(points[ids[2]], m_angle),
                                      m_cfg.maxResidualX / edm4eic::unit::mm,
                                      m_cfg.maxResidualY / edm4eic::unit::mm);
      if (score) {
        m_candidates.push_back({*score, ids});
      }
    }
  }

  void filterTripletsMiddleFixed(const SpacePoints&, Seeds& output) const override {
    std::sort(m_candidates.begin(), m_candidates.end());
    m_candidates.erase(std::unique(m_candidates.begin(), m_candidates.end()), m_candidates.end());
    const auto count = std::min(m_candidates.size(), m_cfg.maxSeedsPerMiddle);
    for (std::size_t i = 0; i < count; ++i) {
      auto seed = output.createSeed();
      seed.assignSpacePointIndices(m_candidates[i].second);
      seed.quality() = static_cast<float>(-m_candidates[i].first);
      // This API field is not a B0 vertex measurement. No vertex cut uses it,
      // and the parameter estimator below ignores it completely.
      seed.vertexZ() = std::numeric_limits<float>::quiet_NaN();
    }
    m_candidates.clear();
  }

private:
  double m_angle;
  const B0TelescopeSeedingConfig& m_cfg;
  // ACTS specifies a stateful const filter interface. This object is event-local.
  mutable std::vector<std::pair<double, std::array<std::uint32_t, 3>>> m_candidates;
};

} // namespace

void findTelescopeSeeds(const SpacePoints& points, const std::array<std::uint32_t, 5>& offsets,
                        double detectorRotation, const B0TelescopeSeedingConfig& cfg,
                        Seeds& output) {
  if (offsets.front() != 0 || offsets.back() != points.size() ||
      !std::is_sorted(offsets.begin(), offsets.end())) {
    throw std::invalid_argument("B0 station offsets do not partition the space points");
  }
  TelescopeDoubletFinder bottom(true, detectorRotation, cfg), top(false, detectorRotation, cfg);
  TelescopeTripletFinder triplet(detectorRotation, cfg);
  TelescopeSeedFilter filter(detectorRotation, cfg);
  Acts::TripletSeeder seeder;
  Acts::TripletSeeder::Cache cache;
  for (const auto& stations : stationTriples) {
    std::array<SpacePoints::ConstRange, 1> bottoms{
        points.range({offsets[stations[0]], offsets[stations[0] + 1]})};
    const auto middle = points.range({offsets[stations[1]], offsets[stations[1] + 1]});
    std::array<SpacePoints::ConstRange, 1> tops{
        points.range({offsets[stations[2]], offsets[stations[2] + 1]})};
    // Explicit station groups, not a cylindrical grid. Radius/cot(theta) sorted
    // fast paths are disabled in both policies; radiusRangeForMiddle is unused.
    seeder.createSeedsFromGroups(cache, bottom, top, triplet, filter, points,
                                  std::span<SpacePoints::ConstRange>(bottoms), middle,
                                  std::span<SpacePoints::ConstRange>(tops), {0.F, 0.F}, output);
  }
}

} // namespace eicrecon::b0telescope
#endif
