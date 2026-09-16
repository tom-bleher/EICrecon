// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#pragma once
#include "algorithms/tracking/B0TelescopeSeeding.h"
#include "algorithms/tracking/B0TelescopeTruthSeeding.h"
#include "algorithms/tracking/B0TelescopeTracking.h"
#include "algorithms/tracking/B0TelescopeTrackExport.h"
#include "extensions/jana/JOmniFactory.h"
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

namespace eicrecon {

class B0TelescopeSeeding_factory : public JOmniFactory<B0TelescopeSeeding_factory, B0TelescopeSeedingConfig> {
  std::unique_ptr<B0TelescopeSeeding> m_algo;
  PodioInput<edm4hep::EventHeader> m_header{this};
  PodioInput<edm4eic::Measurement2D> m_measurements{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds{this};
  PodioOutput<edm4eic::TrackParameters> m_parameters{this};
  ParameterRef<double> p_stationZGap{this, "stationZGap", config().stationZGap, "B0 station gap [mm]"};
  ParameterRef<bool> p_fieldIntegral{this, "fieldIntegral", config().search.fit.fieldIntegral, "Use experimental field-integral estimator"};
  ParameterRef<unsigned> p_fitIterations{this, "fitIterations", config().search.fit.iterations, "Seed linearization iterations"};
  ParameterRef<unsigned> p_integrationSteps{this, "integrationSteps", config().search.fit.integrationSteps, "RK4 steps per seed interval"};
  ParameterRef<std::size_t> p_minStations{this, "minStations", config().search.minStations, "Minimum physical stations"};
  ParameterRef<std::size_t> p_maxTrials{this, "maxTrials", config().search.maxTrials, "Hard endpoint/candidate/extension work budget"};
  ParameterRef<std::size_t> p_maxCandidates{this, "maxCandidates", config().search.maxCandidates, "Maximum retained candidates"};
  ParameterRef<double> p_minMomentum{this, "minMomentum", config().search.minMomentum, "Minimum candidate momentum [GeV]"};
  ParameterRef<double> p_maxSlope{this, "maxSlope", config().search.maxSlope, "Maximum transverse slope"};
  ParameterRef<double> p_fieldBoundTesla{this, "fieldBoundTesla", config().search.fieldBoundTesla, "Explicit preselection field bound [T]"};
  ParameterRef<double> p_roadWidthMm{this, "roadWidthMm", config().search.roadWidthMm, "Search road width [mm]"};
  ParameterRef<double> p_maxChi2PerDof{this, "maxChi2PerDof", config().search.maxChi2PerDof, "Seed fit chi2/ndf limit"};
  ParameterRef<double> p_extensionChi2{this, "extensionChi2", config().search.extensionChi2, "Extension road chi2 limit"};
  ParameterRef<bool> p_useTiming{this, "useTiming", config().search.useTiming, "Enable optional time-of-flight preselection"};
  ParameterRef<double> p_timingNSigma{this, "timingNSigma", config().search.timingNSigma, "Timing window in sigma"};
  ParameterRef<double> p_timingModelSigmaNs{this, "timingModelSigmaNs", config().search.timingModelSigmaNs, "Timing/path-model allowance [ns]"};
  ParameterRef<double> p_massGeV{this, "massGeV", config().search.massGeV, "Timing mass hypothesis [GeV]"};
  ParameterRef<double> p_curvatureSignificance{this, "curvatureSignificance", config().curvatureSignificance, "Resolved charge significance threshold"};
  ParameterRef<double> p_proposalMaxMomentum{this, "proposalMaxMomentum", config().proposalMaxMomentum, "Finite proposal momentum cap [GeV]"};
  ParameterRef<int> p_forceCharge{this, "forceCharge", config().forceCharge, "0 infer, +/-1 diagnostic hypothesis"};
  ParameterRef<bool> p_bothCharges{this, "bothCharges", config().bothCharges, "Always emit both charge proposals"};
  ParameterRef<double> p_seedPositionSigma{this, "seedPositionSigma", config().seedPositionSigma, "Finding-only position floor [mm]"};
  ParameterRef<double> p_seedAngularSigma{this, "seedAngularSigma", config().seedAngularSigma, "Finding-only slope floor"};
  ParameterRef<double> p_seedQOverPSigma{this, "seedQOverPSigma", config().seedQOverPSigma, "Finding-only q/p floor [GeV^-1]"};
  ParameterRef<std::string> p_stationMapFile{this, "stationMapFile", config().stationMapFile, "Optional station manifest path"};
  ParameterRef<std::string> p_diagnosticsFile{this, "diagnosticsFile", config().diagnosticsFile, "Optional seed JSONL path"};
public:
  void Configure() {
    m_algo = std::make_unique<B0TelescopeSeeding>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }
  void Process(int32_t, uint64_t) {
    m_algo->process({m_header(), m_measurements()}, {m_seeds().get(), m_parameters().get()});
  }
};

class B0TelescopeTruthSeeding_factory : public JOmniFactory<B0TelescopeTruthSeeding_factory, B0TelescopeTruthSeedingConfig> {
  std::unique_ptr<B0TelescopeTruthSeeding> m_algo;
  PodioInput<edm4hep::EventHeader> m_header{this};
  PodioInput<edm4hep::SimTrackerHit> m_simHits{this};
  PodioInput<edm4eic::Measurement2D> m_measurements{this};
  PodioInput<edm4eic::MCRecoTrackerHitAssociation> m_rawAssociations{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds{this};
  PodioOutput<edm4eic::TrackParameters> m_parameters{this};
  ParameterRef<std::size_t> p_minStations{this, "minStations", config().minStations, "Minimum physical crossings"};
  ParameterRef<double> p_stationZGap{this, "stationZGap", config().stationZGap, "B0 station gap [mm]"};
  ParameterRef<double> p_minMomentumGeV{this, "minMomentumGeV", config().minMomentumGeV, "Minimum truth seed momentum [GeV]"};
  ParameterRef<double> p_relativeMomentumSmear{this, "relativeMomentumSmear", config().relativeMomentumSmear, "Log-momentum smear sigma"};
  ParameterRef<double> p_positionSigmaMm{this, "positionSigmaMm", config().positionSigmaMm, "Truth finding position sigma [mm]"};
  ParameterRef<double> p_slopeSigma{this, "slopeSigma", config().slopeSigma, "Truth finding slope sigma"};
  ParameterRef<double> p_qOverPSigma{this, "qOverPSigma", config().qOverPSigma, "Truth finding q/p sigma [GeV^-1]"};
  ParameterRef<double> p_timeSigmaNs{this, "timeSigmaNs", config().timeSigmaNs, "Truth finding time sigma [ns]"};
  ParameterRef<double> p_measurementPurityMin{this, "measurementPurityMin", config().measurementPurityMin, "Strict truth matching purity threshold"};
  ParameterRef<double> p_maxProjectionMm{this, "maxProjectionMm", config().maxProjectionMm, "Maximum local tangent binding distance [mm]"};
  ParameterRef<std::string> p_referenceFile{this, "referenceFile", config().referenceFile, "Optional local truth reference JSONL path"};
public:
  void Configure() {
    m_algo = std::make_unique<B0TelescopeTruthSeeding>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }
  void Process(int32_t, uint64_t) {
    m_algo->process({m_header(), m_simHits(), m_measurements(), m_rawAssociations()},
                    {m_seeds().get(), m_parameters().get()});
  }
};

class B0TelescopeTracking_factory : public JOmniFactory<B0TelescopeTracking_factory, B0TelescopeTrackingConfig> {
  std::unique_ptr<B0TelescopeTracking> m_algo;
  PodioInput<edm4hep::EventHeader> m_header{this};
  PodioInput<edm4eic::TrackSeed> m_seeds{this};
  PodioInput<edm4eic::Measurement2D> m_measurements{this};
  Output<Acts::ConstVectorMultiTrajectory> m_states{this};
  Output<Acts::ConstVectorTrackContainer> m_tracks{this};
  ParameterRef<bool> p_doFinalWeakPriorRefit{this, "doFinalWeakPriorRefit", config().doFinalWeakPriorRefit, "Independent final selected-hit fit"};
  ParameterRef<bool> p_fitSeedMeasurementsDirectly{this, "fitSeedMeasurementsDirectly", config().fitSeedMeasurementsDirectly, "Fit candidate hits without CKF finding"};
  ParameterRef<int> p_particleHypothesisPdg{this, "particleHypothesisPdg", config().particleHypothesisPdg, "Absolute mass-hypothesis PDG, not PID"};
  ParameterRef<bool> p_useSeedParticleHypothesis{this, "useSeedParticleHypothesis", config().useSeedParticleHypothesis, "Use truth seed PDG for validation only"};
  ParameterRef<std::size_t> p_minStations{this, "minStations", config().minStations, "Minimum fitted physical stations"};
  ParameterRef<double> p_stationZGap{this, "stationZGap", config().stationZGap, "B0 station gap [mm]"};
  ParameterRef<double> p_chi2Cut{this, "chi2Cut", config().chi2Cut, "CKF measurement chi2 limit"};
  ParameterRef<std::size_t> p_maxBranchesPerSurface{this, "maxBranchesPerSurface", config().maxBranchesPerSurface, "CKF branch limit per surface"};
  ParameterRef<double> p_weakPriorScale{this, "weakPriorScale", config().weakPriorScale, "Scale fresh absolute weak variances"};
  ParameterRef<double> p_initializationStepMm{this, "initializationStepMm", config().initializationStepMm, "Numerical starting offset upstream [mm]"};
  ParameterRef<bool> p_promptSelection{this, "promptSelection", config().promptSelection, "Post-fit prompt selection, not a beamspot prior"};
  ParameterRef<double> p_promptMaxD0Mm{this, "promptMaxD0Mm", config().promptMaxD0Mm, "Prompt d0 limit [mm]"};
  ParameterRef<double> p_promptMaxZ0Mm{this, "promptMaxZ0Mm", config().promptMaxZ0Mm, "Prompt z0 limit [mm]"};
  ParameterRef<std::string> p_diagnosticsFile{this, "diagnosticsFile", config().diagnosticsFile, "Optional per-chain diagnostics JSONL path"};
  ParameterRef<double> p_weakLoc0{this, "weakPriorLoc0Variance", config().weakPriorVariances[0], "Absolute weak-prior variance [mm^2]"};
  ParameterRef<double> p_weakLoc1{this, "weakPriorLoc1Variance", config().weakPriorVariances[1], "Absolute weak-prior variance [mm^2]"};
  ParameterRef<double> p_weakPhi{this, "weakPriorPhiVariance", config().weakPriorVariances[2], "Absolute weak-prior variance [rad^2]"};
  ParameterRef<double> p_weakTheta{this, "weakPriorThetaVariance", config().weakPriorVariances[3], "Absolute weak-prior variance [rad^2]"};
  ParameterRef<double> p_weakQOverP{this, "weakPriorQOverPVariance", config().weakPriorVariances[4], "Absolute weak-prior variance [GeV^-2]"};
  ParameterRef<double> p_weakTime{this, "weakPriorTimeVariance", config().weakPriorVariances[5], "Absolute weak-prior variance [ns^2]"};
public:
  void Configure() {
    m_algo = std::make_unique<B0TelescopeTracking>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }
  void Process(int32_t, uint64_t) {
    m_algo->process({m_header(), m_seeds(), m_measurements()},
                    {&m_states().emplace_back(), &m_tracks().emplace_back()});
  }
};

class B0TelescopeTrackExport_factory : public JOmniFactory<B0TelescopeTrackExport_factory, NoConfig> {
  std::unique_ptr<B0TelescopeTrackExport> m_algo;
  PodioInput<edm4eic::Measurement2D> m_measurements{this};
  PodioInput<edm4eic::TrackSeed> m_seeds{this};
  Input<Acts::ConstVectorMultiTrajectory> m_states{this};
  Input<Acts::ConstVectorTrackContainer> m_tracks{this};
  PodioInput<edm4eic::MCRecoTrackerHitAssociation> m_rawAssociations{this};
  PodioOutput<edm4eic::Trajectory> m_trajectories{this};
  PodioOutput<edm4eic::TrackParameters> m_parameters{this};
  PodioOutput<edm4eic::Track> m_outputTracks{this};
  PodioOutput<edm4eic::MCRecoTrackParticleLink> m_links{this};
  PodioOutput<edm4eic::MCRecoTrackParticleAssociation> m_associations{this};
public:
  void Configure() {
    m_algo = std::make_unique<B0TelescopeTrackExport>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }
  void Process(int32_t, uint64_t) {
    const auto states = m_states();
    if (states.empty() || !states.front()) throw std::runtime_error("Missing B0 states");
    const auto tracks = m_tracks();
    if (tracks.empty() || !tracks.front()) throw std::runtime_error("Missing B0 tracks");
    m_algo->process({m_measurements(), m_seeds(), states.front(), tracks.front(), m_rawAssociations()},
                    {m_trajectories().get(), m_parameters().get(), m_outputTracks().get(), m_links().get(), m_associations().get()});
  }
};
} // namespace eicrecon
