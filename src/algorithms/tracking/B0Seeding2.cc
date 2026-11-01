// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "B0Seeding2.h"

#include <Acts/ActsVersion.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>
#include <fmt/format.h>
#include <spdlog/logger.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "B0TripletSeedingHelpers.h"
#include "extensions/spdlog/SpdlogToActs.h"

#if B0SEEDING2_HAS_SEEDING2
#include <Acts/Definitions/Direction.hpp>
#include <Acts/EventData/SeedContainer2.hpp>
#include <Acts/EventData/SeedProxy2.hpp>
#include <Acts/EventData/SpacePointColumns.hpp>
#include <Acts/EventData/SpacePointContainer2.hpp>
#include <Acts/EventData/SpacePointProxy2.hpp>
#include <Acts/EventData/Types.hpp>
#include <Acts/Seeding2/BroadTripletSeedFilter.hpp>
#include <Acts/Seeding2/DoubletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeedFinder.hpp>
#include <Acts/Seeding2/TripletSeeder.hpp>
#endif

namespace eicrecon {

namespace {
  enum class FieldModel { Scalar, Effective, PerMiddle };
} // namespace

#if B0SEEDING2_HAS_SEEDING2
namespace {

  struct Finders {
    std::unique_ptr<Acts::DoubletSeedFinder> bottom;
    std::unique_ptr<Acts::DoubletSeedFinder> top;
    std::unique_ptr<Acts::TripletSeedFinder> triplet;
  };

  /// Third pair of the hit-time window: the bottom and the top hit. The two pairs with the middle
  /// hit are applied when the doublet candidates are selected.
  class TimeWindowCuts final : public Acts::ITripletSeedCuts {
  public:
    TimeWindowCuts(const B0TripletSeedingConfig& cfg,
                   const std::vector<b0_triplet_seeding::Point>& points)
        : m_cfg(cfg), m_points(points) {}

    float seedWeight(const Acts::ConstSpacePointProxy2& /*bottom*/,
                     const Acts::ConstSpacePointProxy2& /*middle*/,
                     const Acts::ConstSpacePointProxy2& /*top*/) const override {
      return 0.F;
    }
    bool singleSeedCut(float /*weight*/, const Acts::ConstSpacePointProxy2& bottom,
                       const Acts::ConstSpacePointProxy2& /*middle*/,
                       const Acts::ConstSpacePointProxy2& top) const override {
      return b0_triplet_seeding::inTimeWindow(m_cfg, m_points[bottom.copyFromIndex()],
                                              m_points[top.copyFromIndex()]);
    }
    void cutPerMiddleSp(std::span<Acts::TripletCandidate2>& /*candidates*/) const override {}

  private:
    const B0TripletSeedingConfig& m_cfg;
    const std::vector<b0_triplet_seeding::Point>& m_points;
  };

} // namespace
#endif

struct B0Seeding2::Impl {
  FieldModel model = FieldModel::Effective;
  Acts::Vector3 beamAxis{0., 0., 1.};
  Acts::Vector3 fixedFieldAxis{0., 1., 0.};
  Acts::RotationMatrix3 fixedRotation = Acts::RotationMatrix3::Identity();
#if B0SEEDING2_HAS_SEEDING2
  Finders fixed;
  std::vector<Finders> binned;
  Acts::BroadTripletSeedFilter::Config filterConfig;
  std::unique_ptr<const Acts::Logger> actsLogger;
  Acts::TripletSeeder seeder;
#endif
};

Acts::RotationMatrix3 B0Seeding2::frameRotation(const Acts::Vector3& fieldAxis,
                                                const Acts::Vector3& beamAxis) {
  const Acts::Vector3 zAxis = fieldAxis.normalized();
  const Acts::Vector3 xAxis = (beamAxis - beamAxis.dot(zAxis) * zAxis).normalized();
  const Acts::Vector3 yAxis = zAxis.cross(xAxis);
  Acts::RotationMatrix3 rotation;
  rotation.row(0) = xAxis.transpose();
  rotation.row(1) = yAxis.transpose();
  rotation.row(2) = zAxis.transpose();
  return rotation;
}

void B0Seeding2::init() {
#if !B0SEEDING2_HAS_SEEDING2
  throw std::runtime_error("B0Seeding2 needs Acts Seeding2 (Acts 45.3 to 46.x)");
#else
  if (!m_geoSvc || !m_geoSvc->getFieldProvider()) {
    throw std::runtime_error("B0Seeding2: no magnetic field provider");
  }
  auto impl = std::make_shared<Impl>();
  if (m_cfg.fieldModel == "scalar") {
    impl->model = FieldModel::Scalar;
  } else if (m_cfg.fieldModel == "effective") {
    impl->model = FieldModel::Effective;
  } else if (m_cfg.fieldModel == "perMiddle") {
    impl->model = FieldModel::PerMiddle;
  } else {
    throw std::runtime_error(fmt::format(
        "B0Seeding2: fieldModel must be scalar, effective or perMiddle, got {}", m_cfg.fieldModel));
  }
  if (!(m_cfg.bFieldInZ > 0.) || !(m_cfg.bFieldBinWidth > 0.) || !(m_cfg.minPt > 0.) ||
      m_cfg.triplet.maxSeeds == 0) {
    throw std::runtime_error("B0Seeding2: bFieldInZ, bFieldBinWidth, minPt and maxSeeds must be "
                             "positive");
  }

  // The frame: x' along the ion beam, z' along the field (global y for the dipole)
  const auto* detector = m_geoSvc->dd4hepDetector();
  double angle         = m_cfg.frameAngle;
  if (!std::isfinite(angle)) {
    angle = detector->constant<double>("CrossingAngle") / dd4hep::rad;
  }
  impl->beamAxis = Acts::Vector3{std::sin(angle), 0., std::cos(angle)};
  // Orient z' with the field at the centre of the B0 magnet, so that the scalar field model has
  // the right sign. The cuts do not depend on the orientation.
  try {
    const Acts::Vector3 centre{
        detector->constant<double>("B0PF_XPosition") / dd4hep::mm * Acts::UnitConstants::mm, 0.,
        detector->constant<double>("B0PF_CenterPosition") / dd4hep::mm * Acts::UnitConstants::mm};
    auto cache = m_geoSvc->getFieldProvider()->makeCache(m_geoSvc->getActsMagneticFieldContext());
    const auto field = m_geoSvc->getFieldProvider()->getField(centre, cache);
    if (field.ok() && field->y() < 0.) {
      impl->fixedFieldAxis = Acts::Vector3{0., -1., 0.};
    }
    info("B0Seeding2: crossing angle {} rad, field at the B0 magnet centre ({}, 0, {}) mm = ({}, "
         "{}, {}) T",
         angle, centre.x(), centre.z(), field.ok() ? field->x() / Acts::UnitConstants::T : 0.,
         field.ok() ? field->y() / Acts::UnitConstants::T : 0.,
         field.ok() ? field->z() / Acts::UnitConstants::T : 0.);
  } catch (const std::exception& e) {
    warning("B0Seeding2: no B0 magnet constants ({}), frame z' along +y", e.what());
  }
  impl->fixedRotation = frameRotation(impl->fixedFieldAxis, impl->beamAxis);

  auto makeFinders = [&](double bFieldTesla) {
    const auto bField = static_cast<float>(bFieldTesla * Acts::UnitConstants::T);
    Acts::DoubletSeedFinder::Config bottomConfig;
    bottomConfig.spacePointsSortedByRadius    = false;
    bottomConfig.candidateDirection           = Acts::Direction::Backward();
    bottomConfig.deltaRMin                    = m_cfg.deltaRMinBottom;
    bottomConfig.deltaRMax                    = m_cfg.deltaRMaxBottom;
    bottomConfig.deltaZMin                    = m_cfg.deltaZMin;
    bottomConfig.deltaZMax                    = m_cfg.deltaZMax;
    bottomConfig.impactMax                    = m_cfg.impactMax;
    bottomConfig.interactionPointCut          = false;
    bottomConfig.collisionRegionMin           = m_cfg.collisionRegionMin;
    bottomConfig.collisionRegionMax           = m_cfg.collisionRegionMax;
    bottomConfig.cotThetaMax                  = m_cfg.cotThetaMax;
    bottomConfig.minPt                        = m_cfg.minPt;
    bottomConfig.helixCutTolerance            = m_cfg.helixCutTolerance;
    Acts::DoubletSeedFinder::Config topConfig = bottomConfig;
    topConfig.candidateDirection              = Acts::Direction::Forward();
    topConfig.deltaRMin                       = m_cfg.deltaRMinTop;
    topConfig.deltaRMax                       = m_cfg.deltaRMaxTop;
    Acts::TripletSeedFinder::Config tripletConfig;
    tripletConfig.useStripInfo      = false;
    tripletConfig.sortedByCotTheta  = m_cfg.sortedByCotTheta;
    tripletConfig.minPt             = m_cfg.minPt;
    tripletConfig.sigmaScattering   = m_cfg.sigmaScattering;
    tripletConfig.radLengthPerSeed  = m_cfg.radLengthPerSeed;
    tripletConfig.impactMax         = m_cfg.impactMax;
    tripletConfig.helixCutTolerance = m_cfg.helixCutTolerance;
    return Finders{.bottom = Acts::DoubletSeedFinder::create(
                       Acts::DoubletSeedFinder::DerivedConfig(bottomConfig, bField)),
                   .top = Acts::DoubletSeedFinder::create(
                       Acts::DoubletSeedFinder::DerivedConfig(topConfig, bField)),
                   .triplet = Acts::TripletSeedFinder::create(
                       Acts::TripletSeedFinder::DerivedConfig(tripletConfig, bField))};
  };
  impl->fixed = makeFinders(m_cfg.bFieldInZ);
  if (impl->model == FieldModel::PerMiddle) {
    // |B| up to 4 T, finders at the bin centres
    const auto bins = static_cast<std::size_t>(std::ceil(4. / m_cfg.bFieldBinWidth));
    for (std::size_t i = 0; i < bins; ++i) {
      impl->binned.push_back(makeFinders((static_cast<double>(i) + 0.5) * m_cfg.bFieldBinWidth));
    }
  }

  auto& filter                     = impl->filterConfig;
  filter.deltaInvHelixDiameter     = m_cfg.deltaInvHelixDiameter;
  filter.deltaRMin                 = m_cfg.filterDeltaRMin;
  filter.compatSeedWeight          = m_cfg.compatSeedWeight;
  filter.compatSeedLimit           = m_cfg.compatSeedLimit;
  filter.impactWeightFactor        = m_cfg.impactWeightFactor;
  filter.zOriginWeightFactor       = m_cfg.zOriginWeightFactor;
  filter.maxSeedsPerSpM            = m_cfg.maxSeedsPerSpM;
  filter.maxSeedsPerSpMConf        = m_cfg.maxSeedsPerSpMConf;
  filter.maxQualitySeedsPerSpMConf = m_cfg.maxQualitySeedsPerSpMConf;
  filter.seedConfirmation          = false;

  auto spdlogLogger = spdlog::default_logger()->clone(std::string(this->name()));
  spdlogLogger->set_level(static_cast<spdlog::level::level_enum>(this->level()));
  impl->actsLogger = eicrecon::getSpdlogLogger(std::string(this->name()), spdlogLogger);
  m_impl           = impl;
#endif
}

void B0Seeding2::process(const Input& input, const Output& output) const {
#if B0SEEDING2_HAS_SEEDING2
  using namespace b0_triplet_seeding;

  const auto [hits]                = input;
  auto [track_seeds, track_params] = output;
  const auto& impl                 = *m_impl;
  const auto& tcfg                 = m_cfg.triplet;

  const auto& gctx          = m_geoSvc->getActsGeometryContext();
  const auto& mctx          = m_geoSvc->getActsMagneticFieldContext();
  const auto field          = m_geoSvc->getFieldProvider();
  auto fieldCache           = field->makeCache(mctx);
  const FieldLookup fieldAt = [&](const Acts::Vector3& position) -> Acts::Vector3 {
    const auto value = field->getField(position, fieldCache);
    return value.ok() ? *value : Acts::Vector3::Zero();
  };

  // The same hits, order and units as B0TripletSeeding, with the field at every hit
  std::vector<Point> points;
  for (auto& station : groupStations(*hits, tcfg.stationGap)) {
    for (auto& point : station) {
      point.field = fieldAt(point.position);
      if (!point.field.isZero()) {
        point.fieldDirection = point.field.normalized();
      }
      points.push_back(point);
    }
  }

  // Space points in a frame, with the indices into points
  auto fill = [&](const Acts::RotationMatrix3& rotation, Acts::SpacePointContainer2& container) {
    container.clear();
    container.reserve(points.size());
    for (std::uint32_t i = 0; i < points.size(); ++i) {
      const Acts::Vector3 primed = rotation * points[i].position;
      const double r             = std::hypot(primed.x(), primed.y());
      if (r < m_cfg.rMin || r > m_cfg.rMax || primed.z() < m_cfg.zMin || primed.z() > m_cfg.zMax) {
        continue;
      }
      auto sp  = container.createSpacePoint();
      sp.xy()  = {static_cast<float>(primed.x()), static_cast<float>(primed.y())};
      sp.zr()  = {static_cast<float>(primed.z()), static_cast<float>(r)};
      sp.phi() = static_cast<float>(std::atan2(primed.y(), primed.x()));
      sp.varianceZ() =
          static_cast<float>(points[i].variance + m_cfg.extraSigmaZ * m_cfg.extraSigmaZ);
      sp.varianceR()     = static_cast<float>(points[i].variance);
      sp.copyFromIndex() = i;
    }
  };
  const auto columns = Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
                       Acts::SpacePointColumns::Phi | Acts::SpacePointColumns::VarianceZ |
                       Acts::SpacePointColumns::VarianceR | Acts::SpacePointColumns::CopyFromIndex;
  Acts::SpacePointContainer2 fixedSpacePoints(columns);
  fill(impl.fixedRotation, fixedSpacePoints);
  Acts::SpacePointContainer2 localSpacePoints(columns);

  auto filterConfig           = impl.filterConfig;
  filterConfig.experimentCuts = std::make_shared<TimeWindowCuts>(tcfg, points);
  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter seedFilter(filterConfig, filterState, filterCache, *impl.actsLogger);
  Acts::TripletSeeder::Cache seederCache;

  struct Found {
    float quality;
    std::array<std::uint32_t, 3> points;
  };
  std::vector<Found> found;
  std::vector<Acts::SpacePointIndex2> bottoms;
  std::vector<Acts::SpacePointIndex2> tops;
  std::size_t actsSeeds = 0;

  // Every hit in the middle range is a middle candidate, in the order of the hits
  for (Acts::SpacePointIndex2 m = 0; m < fixedSpacePoints.size(); ++m) {
    const auto fixedMiddle = fixedSpacePoints.at(m).asConst();
    const float rMFixed    = fixedMiddle.zr()[1];
    if (rMFixed < m_cfg.rMinMiddle || rMFixed > m_cfg.rMaxMiddle) {
      continue;
    }
    const std::uint32_t middleIndex = fixedMiddle.copyFromIndex();
    const Point& middlePoint        = points[middleIndex];

    const Acts::SpacePointContainer2* spacePoints = &fixedSpacePoints;
    const Finders* finders                        = &impl.fixed;
    Acts::SpacePointIndex2 middle                 = m;
    if (impl.model == FieldModel::PerMiddle) {
      const double bLocal = middlePoint.field.norm() / Acts::UnitConstants::T;
      if (!(bLocal > 0.)) {
        continue;
      }
      const auto bin =
          std::min(static_cast<std::size_t>(bLocal / m_cfg.bFieldBinWidth), impl.binned.size() - 1);
      finders = &impl.binned[bin];
      fill(frameRotation(middlePoint.field, impl.beamAxis), localSpacePoints);
      spacePoints      = &localSpacePoints;
      bool foundMiddle = false;
      for (Acts::SpacePointIndex2 i = 0; i < localSpacePoints.size(); ++i) {
        if (localSpacePoints.at(i).asConst().copyFromIndex() == middleIndex) {
          middle      = i;
          foundMiddle = true;
          break;
        }
      }
      if (!foundMiddle) {
        continue;
      }
    }
    const auto spM = spacePoints->at(middle);
    const float rM = spM.zr()[1];

    // Doublet candidates in r', with the hit-time window to the middle hit (the doublet stage
    // has no time cut)
    bottoms.clear();
    tops.clear();
    for (Acts::SpacePointIndex2 o = 0; o < spacePoints->size(); ++o) {
      if (o == middle) {
        continue;
      }
      const auto other        = spacePoints->at(o);
      const Point& otherPoint = points[other.copyFromIndex()];
      if (other.zr()[1] < rM) {
        if (inTimeWindow(tcfg, otherPoint, middlePoint)) {
          bottoms.push_back(o);
        }
      } else if (inTimeWindow(tcfg, middlePoint, otherPoint)) {
        tops.push_back(o);
      }
    }
    if (bottoms.empty() || tops.empty()) {
      continue;
    }
    auto bottomSps = spacePoints->subset(std::span<const Acts::SpacePointIndex2>(bottoms));
    auto topSps    = spacePoints->subset(std::span<const Acts::SpacePointIndex2>(tops));
    Acts::SeedContainer2 seeds;
    seeds.assignSpacePointContainer(*spacePoints);
    impl.seeder.createSeedsFromGroup(seederCache, *finders->bottom, *finders->top,
                                     *finders->triplet, seedFilter, *spacePoints, bottomSps, spM,
                                     topSps, seeds);
    for (const auto& seed : seeds) {
      const auto& indices = seed.spacePointIndices();
      found.push_back({.quality = seed.quality(),
                       .points  = {spacePoints->at(indices[0]).copyFromIndex(),
                                   spacePoints->at(indices[1]).copyFromIndex(),
                                   spacePoints->at(indices[2]).copyFromIndex()}});
      ++actsSeeds;
    }
  }

  // Rank by the Acts seed weight, ties by the hit indices, and write up to maxSeeds
  std::ranges::sort(found, [&](const Found& lhs, const Found& rhs) {
    auto key = [&](const Found& f) {
      return std::tuple{-f.quality, points[f.points[0]].index, points[f.points[1]].index,
                        points[f.points[2]].index};
    };
    return key(lhs) < key(rhs);
  });
  const Acts::Vector3 scalarField = m_cfg.bFieldInZ * Acts::UnitConstants::T * impl.fixedFieldAxis;
  std::size_t written             = 0;
  for (const auto& f : found) {
    if (written >= tcfg.maxSeeds) {
      break;
    }
    const Point& a = points[f.points[0]];
    const Point& b = points[f.points[1]];
    const Point& c = points[f.points[2]];
    Acts::Vector3 bField;
    if (impl.model == FieldModel::Scalar) {
      bField = scalarField;
    } else if (tcfg.fieldSamples > 0) {
      bField = fieldIntegral(a, b, c, tcfg.fieldSamples, fieldAt).field;
    } else {
      bField = b.field;
    }
    if (bField.isZero()) {
      continue;
    }
    const auto seed = makeSeed(tcfg, gctx, a, b, c, bField);
    if (!seed) {
      continue;
    }
    // The quality written is the Acts seed weight
    const Candidate candidate{
        .residual = -static_cast<double>(f.quality), .points = {&a, &b, &c}, .field = bField};
    writeSeed(*seed, candidate, *hits, *track_seeds, *track_params);
    ++written;
  }

  debug("{} hits, {} space points, {} Acts seeds, {} seeds", points.size(), fixedSpacePoints.size(),
        actsSeeds, written);
#else
  (void)input;
  (void)output;
#endif
}

} // namespace eicrecon
