// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <memory>

#include "algorithms/tracking/B0TrackerStubSeeder.h"
#include "algorithms/tracking/B0TrackerStubSeederConfig.h"
#include "extensions/jana/JOmniFactory.h"

namespace eicrecon {

class B0TrackerStubSeeder_factory
    : public JOmniFactory<B0TrackerStubSeeder_factory, B0TrackerStubSeederConfig> {

private:
  using AlgoT = eicrecon::B0TrackerStubSeeder;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::TrackerHit> m_hits_input{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds_output{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};

  ParameterRef<float> m_minAbsFieldY{this, "minAbsFieldY", config().minAbsFieldY,
                                     "minimum |By| from the ACTS field provider [T]"};
  ParameterRef<float> m_zFieldEntrance{
      this, "zFieldEntrance", config().zFieldEntrance,
      "ion-frame z of B0pf field entrance [mm]; <=0 uses B0PF geometry"};
  ParameterRef<unsigned int> m_fieldSamples{this, "fieldSamples", config().fieldSamples,
                                            "minimum ACTS field mesh points per candidate"};
  ParameterRef<unsigned int> m_fieldFitIterations{this, "fieldFitIterations",
                                                  config().fieldFitIterations,
                                                  "field-integral fit iterations per candidate"};
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
                                     "seed time variance"};
  ParameterRef<float> m_seedTime{this, "seedTime", config().seedTime,
                                 "seed time prior [ns] (no B0 timing digitization yet)"};
  ParameterRef<float> m_phiModelVariance{this, "phiModelVariance", config().phiModelVariance,
                                         "calibrated residual phi variance"};
  ParameterRef<float> m_phiQOverPScale{this, "phiQOverPScale", config().phiQOverPScale,
                                       "q/p-scaled phi model uncertainty [GeV rad]"};
  ParameterRef<float> m_thetaModelVariance{this, "thetaModelVariance", config().thetaModelVariance,
                                           "calibrated residual theta variance"};
  ParameterRef<float> m_thetaQOverPScale{this, "thetaQOverPScale", config().thetaQOverPScale,
                                         "q/p-scaled theta model uncertainty [GeV rad]"};
  ParameterRef<float> m_qOverPModelVariance{this, "qOverPModelVariance",
                                            config().qOverPModelVariance,
                                            "calibrated residual q/p variance"};
  ParameterRef<float> m_qOverPRelativeUncertainty{this, "qOverPRelativeUncertainty",
                                                  config().qOverPRelativeUncertainty,
                                                  "relative q/p trajectory-model uncertainty"};
  ParameterRef<float> m_fieldRelativeUncertainty{this, "fieldRelativeUncertainty",
                                                 config().fieldRelativeUncertainty,
                                                 "external relative field-integral uncertainty"};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void ChangeRun(int32_t /* run_number */) {}

  void Process(int32_t /* run_number */, uint64_t /* event_number */) {
    m_algo->process({m_hits_input()}, {m_seeds_output().get(), m_trackparams_output().get()});
  }
};

} // namespace eicrecon
