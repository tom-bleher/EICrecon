// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include "B0TelescopeTracking.h"
#include "B0WeakPrior.h"
#include "extensions/spdlog/SpdlogToActs.h"
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#if Acts_VERSION_MAJOR >= 46
#include <Acts/EventData/BoundTrackParameters.hpp>
#else
#include <Acts/EventData/GenericBoundTrackParameters.hpp>
#endif
#include <Acts/EventData/ProxyAccessor.hpp>
#include <Acts/Propagator/ActorList.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilterExtensions.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StandardAborters.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/TrackFinding/TrackStateCreator.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#if __has_include(<Acts/TrackFitting/MbfSmoother.hpp>)
#include <Acts/TrackFitting/MbfSmoother.hpp>
#else
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#endif
#include <Acts/Utilities/TrackHelpers.hpp>
#include <ActsExamples/EventData/GeometryContainers.hpp>
#include <ActsExamples/EventData/IndexSourceLink.hpp>
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace eicrecon {
namespace {
using Container = ActsExamples::TrackContainer;
using Track = Container::TrackProxy;
using Transport = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
using ExtrapolationOptions = Transport::Options<Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>;
using Fitter = Acts::KalmanFitter<Transport, Acts::VectorMultiTrajectory>;
using BoundCovariance = std::decay_t<decltype(std::declval<Acts::BoundTrackParameters>().covariance().value())>;
const std::array<double, 6> units{Acts::UnitConstants::mm, Acts::UnitConstants::mm, 1, 1,
                                  1 / Acts::UnitConstants::GeV, Acts::UnitConstants::ns};
template <typename State>
bool isMeasurement(const State& state) {
  const auto f = state.typeFlags();
#if Acts_VERSION_MAJOR >= 45
  return f.isMeasurement() && !f.isOutlier() && !f.isHole();
#else
  return f.test(Acts::TrackStateFlag::MeasurementFlag) && !f.test(Acts::TrackStateFlag::OutlierFlag) &&
         !f.test(Acts::TrackStateFlag::HoleFlag);
#endif
}
std::vector<std::size_t> selectedIndices(const Track& track) {
  std::vector<std::size_t> result;
  for (const auto& state : track.trackStatesReversed()) if (isMeasurement(state))
    result.push_back(state.getUncalibratedSourceLink().template get<ActsExamples::IndexSourceLink>().index());
  std::sort(result.begin(), result.end());
  return result;
}
struct Calibrator {
  const edm4eic::Measurement2DCollection& measurements;
  const Acts::TrackingGeometry& geometry;
  const Acts::Surface* surface(const Acts::SourceLink& link) const {
    return geometry.findSurface(link.get<ActsExamples::IndexSourceLink>().geometryId());
  }
  void calibrate(const Acts::GeometryContext&, const Acts::CalibrationContext&,
                 const Acts::SourceLink& link, Acts::VectorMultiTrajectory::TrackStateProxy state) const {
    const auto index = link.get<ActsExamples::IndexSourceLink>().index();
    if (index >= measurements.size()) throw std::out_of_range("B0 source-link index");
    const auto m = measurements[index];
    const auto c = m.getCovariance();
    Acts::Vector2 local(m.getLoc().a * units[0], m.getLoc().b * units[1]);
    Eigen::Matrix2d cov;
    cov << c.xx, c.xy, c.xy, c.yy; cov *= units[0] * units[0];
    state.setUncalibratedSourceLink(Acts::SourceLink{link});
    state.allocateCalibrated(local, cov);
    state.setProjectorSubspaceIndices(std::array<std::uint8_t, 2>{0, 1});
  }
};
}

void B0TelescopeTracking::init() {
  m_provider = algorithms::ActsSvc::instance().acts_geometry_provider();
  if (!m_provider) throw std::runtime_error("B0 ACTS geometry unavailable");
  m_geometry = std::make_unique<b0::TelescopeGeometry>(*m_provider, m_cfg.stationZGap);
  CKFTracking::makeParticleHypothesis(m_cfg.particleHypothesisPdg);
  b0::weakPrior(m_cfg.weakPriorVariances, m_cfg.weakPriorScale);
  if (m_cfg.minStations < 3 || m_cfg.maxBranchesPerSurface == 0 || !(m_cfg.chi2Cut > 0) ||
      !(m_cfg.initializationStepMm > 0) || !std::isfinite(m_cfg.initializationStepMm) ||
      !(m_cfg.promptMaxD0Mm > 0) || !(m_cfg.promptMaxZ0Mm > 0) ||
      !std::isfinite(m_cfg.chi2Cut) || !std::isfinite(m_cfg.promptMaxD0Mm) ||
      !std::isfinite(m_cfg.promptMaxZ0Mm))
    throw std::invalid_argument("Invalid B0 telescope fitting configuration");
  m_logger = Acts::getDefaultLogger(std::string(name()),
      SpdlogToActsLevel(static_cast<spdlog::level::level_enum>(level())));
  m_finder = CKFTracking::makeCKFTrackingFunction(m_provider->trackingGeometry(), m_provider->getFieldProvider(), *m_logger);
  if (!m_cfg.diagnosticsFile.empty()) {
    std::ofstream out(m_cfg.diagnosticsFile);
    if (!out) throw std::runtime_error("Cannot initialize B0 diagnostic file");
    out << "{\"type\":\"configuration\",\"schema_version\":1,\"independent_refit\":"
        << (m_cfg.doFinalWeakPriorRefit || m_cfg.fitSeedMeasurementsDirectly ? "true" : "false")
        << ",\"direct_fit\":" << (m_cfg.fitSeedMeasurementsDirectly ? "true" : "false")
        << ",\"weak_prior_scale\":" << m_cfg.weakPriorScale << "}\n";
  }
}

void B0TelescopeTracking::process(const Input& input, const Output& output) const {
  const auto startTime = std::chrono::steady_clock::now();
  const auto [headers, seeds, measurements] = input;
  auto [outputStates, outputTracks] = output;
  if (headers->empty()) throw std::runtime_error("B0 diagnostics require EventHeader");
  const auto event = (*headers)[0].getEventNumber();
  const auto& gctx = m_provider->getActsGeometryContext();
  const auto& mctx = m_provider->getActsMagneticFieldContext();
  const auto& cctx = m_provider->getActsCalibrationContext();
  const auto field = m_provider->getFieldProvider();
  auto makeContainer = [] {
    Container c(std::make_shared<Acts::VectorTrackContainer>(), std::make_shared<Acts::VectorMultiTrajectory>());
    c.addColumn<unsigned int>("seed");
    c.addColumn<unsigned int>("b0IndependentRefit");
    return c;
  };
  auto tracks = makeContainer(), foundTracks = makeContainer(), fittedTracks = makeContainer();
  Acts::ProxyAccessor<unsigned int> seedColumn("seed"), refitColumn("b0IndependentRefit");
  Calibrator calibrator{*measurements, *m_provider->trackingGeometry()};
  ActsExamples::GeometryIdMultiset<ActsExamples::IndexSourceLink> links;
  for (std::size_t i = 0; i < measurements->size(); ++i)
    links.emplace(Acts::GeometryIdentifier{(*measurements)[i].getSurface()}, i);
  ActsExamples::IndexSourceLinkAccessor accessor; accessor.container = &links;
  Acts::GainMatrixUpdater updater;
  const Acts::MeasurementSelector::Config selectorConfig{{Acts::GeometryIdentifier(),
      {.etaBins = {}, .chi2CutOff = {m_cfg.chi2Cut}, .numMeasurementsCutOff = {m_cfg.maxBranchesPerSurface}}}};
  Acts::MeasurementSelector selector(selectorConfig);
  using Creator = Acts::TrackStateCreator<ActsExamples::IndexSourceLinkAccessor::Iterator, Container>;
  Creator creator;
  creator.sourceLinkAccessor.connect<&ActsExamples::IndexSourceLinkAccessor::range>(&accessor);
  creator.calibrator.connect<&Calibrator::calibrate>(&calibrator);
  creator.measurementSelector.connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(&selector);
  Acts::CombinatorialKalmanFilterExtensions<Container> findExtensions;
  findExtensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);
  findExtensions.createTrackStates.connect<&Creator::createTrackStates>(&creator);
  Acts::PropagatorPlainOptions propagation(gctx, mctx); propagation.maxSteps = 10000;
  CKFTracking::TrackFinderOptions findOptions(gctx, mctx, cctx, findExtensions, propagation);
  Acts::Navigator::Config nav{m_provider->trackingGeometry()};
  nav.resolveSensitive = true; nav.resolveMaterial = true; nav.resolvePassive = false;
  Transport transport(Acts::EigenStepper<>(field), Acts::Navigator(nav, m_logger->cloneWithSuffix("Navigator")),
                       m_logger->cloneWithSuffix("Propagator"));
  Fitter fitter(Transport(Acts::EigenStepper<>(field), Acts::Navigator(nav)), m_logger->cloneWithSuffix("Refit"));
#if __has_include(<Acts/TrackFitting/MbfSmoother.hpp>)
  using Smoother = Acts::MbfSmoother;
#else
  using Smoother = Acts::GainMatrixSmoother;
#endif
  Smoother smoother;
  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> fitExtensions;
  fitExtensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);
  fitExtensions.smoother.connect<&Smoother::operator()<Acts::VectorMultiTrajectory>>(&smoother);
  fitExtensions.calibrator.connect<&Calibrator::calibrate>(&calibrator);
  fitExtensions.surfaceAccessor.connect<&Calibrator::surface>(&calibrator);
  ExtrapolationOptions extrapolation(gctx, mctx);
  auto origin = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
  std::size_t invalidSeeds = 0, findingFailures = 0, refitFailures = 0, changedSets = 0, stationFailures = 0;
  std::size_t smoothingFailures = 0, promptFailures = 0;
  std::ostringstream records; records << std::setprecision(17);

  for (std::size_t iseed = 0; iseed < seeds->size(); ++iseed) {
    const auto seed = (*seeds)[iseed]; const auto p = seed.getParams();
    const auto* reference = m_geometry->surface(p.getSurface());
    if (!reference) { ++invalidSeeds; continue; }
    Acts::BoundVector values;
    values << p.getLoc().a, p.getLoc().b, p.getPhi(), p.getTheta(), p.getQOverP(), p.getTime();
    BoundCovariance cov;
    for (int i = 0; i < 6; ++i) {
      values(i) *= units[i];
      for (int j = 0; j < 6; ++j) cov(i, j) = p.getCovariance()(i, j) * units[i] * units[j];
    }
    if (!values.allFinite() || !cov.allFinite() || cov.llt().info() != Eigen::Success) { ++invalidSeeds; continue; }
    const int pdg = m_cfg.useSeedParticleHypothesis ? std::abs(p.getPdg()) : m_cfg.particleHypothesisPdg;
    const auto fitHypothesis = CKFTracking::makeParticleHypothesis(pdg);
    Acts::BoundTrackParameters seedParameters(reference->getSharedPtr(), values, cov, fitHypothesis);
    // Numerical initialization just upstream of the first sensor; no material
    // update here. The final fit traverses each selected measurement once.
    auto transform = reference->transform(gctx);
    transform.translation() -= seedParameters.direction() * m_cfg.initializationStepMm * Acts::UnitConstants::mm;
    auto upstream = Acts::Surface::makeShared<Acts::PlaneSurface>(transform);
    using Initializer = Acts::Propagator<Acts::EigenStepper<>>;
    Initializer initializer{Acts::EigenStepper<>(field)};
    Initializer::Options<> initialOptions(gctx, mctx);
    initialOptions.direction = Acts::Direction::Backward(); initialOptions.maxSteps = 1000;
    const auto moved = initializer.propagate(seedParameters, *upstream, initialOptions);
    if (!moved.ok() || !moved->endParameters) { ++invalidSeeds; continue; }
    const auto initial = moved->endParameters.value();
    foundTracks.clear(); fittedTracks.clear();

    auto emit = [&](Track track, bool independent) {
      if (!track.parameters().allFinite() || !track.covariance().allFinite() ||
          track.covariance().llt().info() != Eigen::Success) { ++refitFailures; return; }
      const auto indices = selectedIndices(track);
      std::set<std::size_t> stations;
      for (auto index : indices) {
        if (index >= measurements->size()) { ++stationFailures; return; }
        const auto station = m_geometry->stations().stationForSurface((*measurements)[index].getSurface());
        if (!station) { ++stationFailures; return; }
        stations.insert(*station);
      }
      if (stations.size() < m_cfg.minStations) { ++stationFailures; return; }
      if (m_cfg.promptSelection) {
        Acts::BoundTrackParameters local(track.referenceSurface().getSharedPtr(), track.parameters(), track.covariance(), fitHypothesis);
        auto options = extrapolation; options.direction = Acts::Direction::Backward();
        const auto atOrigin = transport.propagate(local, *origin, options);
        if (!atOrigin.ok() || !atOrigin->endParameters ||
            std::abs(atOrigin->endParameters->parameters()(0)) > m_cfg.promptMaxD0Mm * units[0] ||
            std::abs(atOrigin->endParameters->parameters()(1)) > m_cfg.promptMaxZ0Mm * units[1]) {
          ++promptFailures; return;
        }
      }
      seedColumn(track) = static_cast<unsigned int>(iseed); refitColumn(track) = independent ? 1U : 0U;
      auto accepted = tracks.makeTrack(); accepted.copyFrom(track);
      if (!m_cfg.diagnosticsFile.empty()) {
        for (const auto& state : track.trackStatesReversed()) if (isMeasurement(state)) {
          const auto index = state.getUncalibratedSourceLink().template get<ActsExamples::IndexSourceLink>().index();
          const auto m = (*measurements)[index]; const auto c = m.getCovariance();
          const Eigen::Vector2d residual = Eigen::Vector2d(m.getLoc().a, m.getLoc().b) - state.predicted().head<2>() / units[0];
          Eigen::Matrix2d innovation = state.predictedCovariance().topLeftCorner<2, 2>() / (units[0] * units[0]);
          innovation(0, 0) += c.xx; innovation(1, 1) += c.yy; innovation(0, 1) += c.xy; innovation(1, 0) += c.xy;
          if (!residual.allFinite() || !innovation.allFinite()) continue;
          records << "{\"type\":\"innovation\",\"event\":" << event << ",\"seed\":" << iseed
                  << ",\"track\":" << accepted.index() << ",\"surface\":\"" << m.getSurface()
                  << "\",\"measurement\":" << index << ",\"residual_mm\":[" << residual.x() << ',' << residual.y()
                  << "],\"covariance_mm2\":[" << innovation(0, 0) << ',' << innovation(0, 1) << ',' << innovation(1, 1) << "]}\n";
        }
      }
    };
    auto refit = [&](std::vector<std::size_t> indices) {
      std::sort(indices.begin(), indices.end());
      if (std::adjacent_find(indices.begin(), indices.end()) != indices.end() || indices.size() < 3) { ++changedSets; return; }
      std::set<std::uint64_t> surfaces;
      std::vector<Acts::SourceLink> selected;
      for (auto index : indices) {
        if (index >= measurements->size() || !surfaces.insert((*measurements)[index].getSurface()).second) { ++changedSets; return; }
        selected.emplace_back(ActsExamples::IndexSourceLink(Acts::GeometryIdentifier{(*measurements)[index].getSurface()}, index));
      }
      const auto weak = b0::weakPrior(m_cfg.weakPriorVariances, m_cfg.weakPriorScale);
      BoundCovariance prior;
      for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) prior(i, j) = weak(i, j) * units[i] * units[j];
      Acts::BoundTrackParameters start(initial.referenceSurface().getSharedPtr(), initial.parameters(), prior, fitHypothesis);
      Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> options(gctx, mctx, cctx, fitExtensions, propagation, reference);
      options.multipleScattering = true; options.energyLoss = true;
      options.referenceSurfaceStrategy = decltype(options.referenceSurfaceStrategy)::first;
      const auto result = fitter.fit(selected.begin(), selected.end(), start, options, fittedTracks);
      if (!result.ok()) { ++refitFailures; debug("B0 independent refit failed: {}", result.error().message()); return; }
      auto track = result.value();
      if (selectedIndices(track) != indices) { ++changedSets; return; }
      emit(track, true);
    };
    if (m_cfg.fitSeedMeasurementsDirectly) {
      std::set<std::pair<int, int>> seedHits;
      for (const auto& h : seed.getHits()) seedHits.emplace(h.getObjectID().collectionID, h.getObjectID().index);
      std::vector<std::size_t> indices;
      for (std::size_t i = 0; i < measurements->size(); ++i) {
        const auto m = (*measurements)[i];
        bool contained = m.hits_size() > 0;
        for (const auto& h : m.getHits()) contained &= seedHits.count({h.getObjectID().collectionID, h.getObjectID().index}) != 0;
        if (contained) indices.push_back(i);
      }
      refit(std::move(indices));
    } else {
      auto found = (*m_finder)(initial, findOptions, foundTracks);
      if (!found.ok()) { ++findingFailures; continue; }
      if (found.value().empty()) ++findingFailures;
      for (auto& track : found.value()) {
        if (m_cfg.doFinalWeakPriorRefit) { refit(selectedIndices(track)); continue; }
        if (!Acts::smoothTrack(gctx, track, *m_logger).ok() ||
            !Acts::extrapolateTrackToReferenceSurface(track, *reference, transport, extrapolation,
                                                       Acts::TrackExtrapolationStrategy::first, *m_logger).ok()) {
          ++smoothingFailures; continue;
        }
        emit(track, false);
      }
    }
  }
  const double elapsedMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - startTime).count();
  records << "{\"type\":\"event\",\"elapsed_ms\":" << elapsedMs << ",\"event\":" << event << ",\"seeds\":" << seeds->size()
          << ",\"tracks\":" << tracks.size() << ",\"invalid_seeds\":" << invalidSeeds
          << ",\"finding_failures\":" << findingFailures << ",\"refit_failures\":" << refitFailures
          << ",\"changed_measurement_sets\":" << changedSets << ",\"station_failures\":" << stationFailures
          << ",\"smoothing_failures\":" << smoothingFailures << ",\"prompt_selection_failures\":" << promptFailures << "}\n";
  if (!m_cfg.diagnosticsFile.empty()) {
    std::lock_guard lock(m_diagnosticsMutex);
    std::ofstream file(m_cfg.diagnosticsFile, std::ios::app);
    if (!(file << records.str())) throw std::runtime_error("Cannot append B0 diagnostics");
  }
  debug("B0 telescope tracks={} refitFailures={} changedSets={} invalidSeeds={}", tracks.size(), refitFailures, changedSets, invalidSeeds);
  *outputStates = new Acts::ConstVectorMultiTrajectory(std::move(tracks.trackStateContainer()));
  *outputTracks = new Acts::ConstVectorTrackContainer(std::move(tracks.container()));
}
} // namespace eicrecon
