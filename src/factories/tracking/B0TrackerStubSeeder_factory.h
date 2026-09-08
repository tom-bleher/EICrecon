// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <cstdint>
#include <edm4hep/EventHeaderCollection.h>
#include <memory>
#include <stdexcept>
#include <string>

#include "algorithms/tracking/B0CandidateReplay.h"
#include "algorithms/tracking/B0TrackerStubSeeder.h"
#include "algorithms/tracking/B0TrackerStubSeederConfig.h"
#include "extensions/jana/JOmniFactory.h"

namespace eicrecon {

class B0TrackerStubSeeder_factory
    : public JOmniFactory<B0TrackerStubSeeder_factory, B0TrackerStubSeederConfig> {

private:
  using AlgoT = eicrecon::B0TrackerStubSeeder;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4hep::EventHeader> m_headers_input{this};
  PodioInput<edm4eic::TrackerHit> m_hits_input{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds_output{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};

  ParameterRef<float> m_minAbsFieldY{this, "minAbsFieldY", config().minAbsFieldY,
                                     "minimum |By| from the ACTS field provider [T]"};
  ParameterRef<float> m_zFieldEntrance{
      this, "zFieldEntrance", config().zFieldEntrance,
      "ion-frame z of B0pf field entrance [mm]; <=0 uses B0PF geometry"};
  ParameterRef<float> m_stationZGap{this, "stationZGap", config().stationZGap,
                                    "ion-frame z gap [mm] that starts a new B0 station"};
  ParameterRef<unsigned int> m_minStations{this, "minStations", config().minStations,
                                           "minimum distinct stations per seed"};
  ParameterRef<float> m_yRoadWidth{this, "yRoadWidth", config().yRoadWidth,
                                   "half-width [mm] of the non-bend-plane preselection road"};
  ParameterRef<float> m_maxYBeamlineResidual{
      this, "maxYBeamlineResidual", config().maxYBeamlineResidual,
      "max endpoint residual [mm] relative to a beamline ray"};
  ParameterRef<unsigned int> m_maxCombinations{this, "maxCombinations", config().maxCombinations,
                                               "cap on hit combinations per event"};
  ParameterRef<unsigned int> m_maxSeeds{this, "maxSeeds", config().maxSeeds,
                                        "maximum seeds per event"};
  ParameterRef<unsigned int> m_maxSharedHits{this, "maxSharedHits", config().maxSharedHits,
                                             "max shared hits between accepted seeds"};
  ParameterRef<float> m_sharedHitDistance{
      this, "sharedHitDistance", config().sharedHitDistance,
      "same-station hits closer than this [mm] count as shared"};
  ParameterRef<int> m_charge{this, "charge", config().charge,
                             "0=infer charge, +/-1=force a diagnostic charge hypothesis"};
  ParameterRef<bool> m_testBothCharges{this, "testBothCharges", config().testBothCharges,
                                       "emit both charge hypotheses for each candidate"};
  ParameterRef<float> m_minCurvatureSignificance{
      this, "minCurvatureSignificance", config().minCurvatureSignificance,
      "minimum |q/p| significance for a resolved charge sign"};
  ParameterRef<bool> m_constrainToBeamline{this, "constrainToBeamline",
                                           config().constrainToBeamline,
                                           "use origin->field-entrance chord for the direction"};
  ParameterRef<float> m_pMin{this, "pMin", config().pMin, "min accepted fitted momentum [GeV]"};
  ParameterRef<float> m_maxYResidual{this, "maxYResidual", config().maxYResidual,
                                     "max RMS y(z) compatibility residual [mm]"};
  ParameterRef<float> m_maxXResidual{this, "maxXResidual", config().maxXResidual,
                                     "max RMS x(z) compatibility residual [mm]"};
  ParameterRef<float> m_maxAbsTransverseSlope{this, "maxAbsTransverseSlope",
                                              config().maxAbsTransverseSlope,
                                              "max transverse slope at first B0 station"};
  ParameterRef<float> m_locaVariance{this, "locaVariance", config().locaVariance,
                                     "fallback seed loc0 variance"};
  ParameterRef<float> m_locbVariance{this, "locbVariance", config().locbVariance,
                                     "fallback seed loc1 variance"};
  ParameterRef<float> m_phiVariance{this, "phiVariance", config().phiVariance,
                                    "fallback seed phi variance"};
  ParameterRef<float> m_thetaVariance{this, "thetaVariance", config().thetaVariance,
                                      "fallback seed theta variance"};
  ParameterRef<float> m_qOverPVariance{this, "qOverPVariance", config().qOverPVariance,
                                       "fallback seed q/p variance"};
  ParameterRef<float> m_timeVariance{this, "timeVariance", config().timeVariance,
                                     "seed time variance [ns^2]"};
  ParameterRef<float> m_beamSpotSizeX{this, "beamSpotSizeX", config().beamSpotSizeX,
                                      "interaction-vertex sigma along lab x [mm]"};
  ParameterRef<float> m_beamSpotSizeY{this, "beamSpotSizeY", config().beamSpotSizeY,
                                      "interaction-vertex sigma along lab y [mm]"};
  ParameterRef<float> m_beamSpotSizeZ{this, "beamSpotSizeZ", config().beamSpotSizeZ,
                                      "interaction-vertex sigma along lab z [mm]"};
  ParameterRef<float> m_scatteringScale{
      this, "scatteringScale", config().scatteringScale,
      "scattering angle coefficient [GeV rad], sigma=scale*|q/p|"};
  ParameterRef<float> m_qOverPRelativeUncertainty{this, "qOverPRelativeUncertainty",
                                                  config().qOverPRelativeUncertainty,
                                                  "relative q/p scattering/model uncertainty"};
  ParameterRef<std::string> m_candidateFile{
      this, "candidateFile", config().candidateFile,
      "JSON of external RecHit candidates; empty keeps combinatorial seeding"};

  B0ReplayFile m_replay;

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
    if (!config().candidateFile.empty()) {
      m_replay = loadB0ReplayFile(config().candidateFile);
    }
  }

  void ChangeRun(int32_t /* run_number */) {}

  void Process(int32_t /* run_number */, uint64_t /* event_number */) {
    if (config().candidateFile.empty()) {
      m_algo->process({m_hits_input()}, {m_seeds_output().get(), m_trackparams_output().get()});
      return;
    }
    const auto headers = m_headers_input();
    if (headers->size() != 1) {
      throw std::runtime_error(
          "B0TrackerStubSeeder: expected one EventHeader for candidate replay");
    }
    const std::int64_t run   = static_cast<std::int64_t>(headers->at(0).getRunNumber());
    const std::int64_t event = static_cast<std::int64_t>(headers->at(0).getEventNumber());
    const auto it            = m_replay.events.find({run, event});
    if (it == m_replay.events.end()) {
      throw std::runtime_error("B0TrackerStubSeeder: no replay candidates for run/event (" +
                               std::to_string(run) + ", " + std::to_string(event) + ")");
    }
    auto bound = bindReplayCandidates(*m_hits_input(), it->second.candidates);
    m_algo->process({m_hits_input()}, {m_seeds_output().get(), m_trackparams_output().get()},
                    bound);
  }
};

} // namespace eicrecon
