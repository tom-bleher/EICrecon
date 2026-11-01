// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

#include "algorithms/tracking/B0Seeding2.h"
#include "algorithms/tracking/B0Seeding2Config.h"
#include "extensions/jana/JOmniFactory.h"
#include "services/geometry/acts/ACTSGeo_service.h"

namespace eicrecon {

class B0Seeding2_factory : public JOmniFactory<B0Seeding2_factory, B0Seeding2Config> {

private:
  using AlgoT = eicrecon::B0Seeding2;
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::TrackerHit> m_hits_input{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds_output{this};
  PodioOutput<edm4eic::TrackParameters> m_trackparams_output{this};

  // Declares the dependency on the ACTS geometry and field (see B0TripletSeeding_factory)
  Service<ACTSGeo_service> m_ACTSGeoSvc{this};

  // Shared with B0TripletSeeding (same names and meaning)
  ParameterRef<double> m_stationGap{this, "stationGap", config().triplet.stationGap,
                                    "z gap that starts a new station, mm"};
  ParameterRef<double> m_minMomentum{this, "minMomentum", config().triplet.minMomentum,
                                     "min seed momentum in the seed field, GeV"};
  ParameterRef<double> m_anchorDistance{this, "anchorDistance", config().triplet.anchorDistance,
                                        "seed anchor upstream of the first hit, mm"};
  ParameterRef<std::size_t> m_maxSeeds{this, "maxSeeds", config().triplet.maxSeeds,
                                       "max seeds per event"};
  ParameterRef<double> m_timeWindowSigmas{this, "timeWindowSigmas",
                                          config().triplet.timeWindowSigmas,
                                          "hit-time window in sigmas, 0 off"};
  ParameterRef<double> m_timeWindowTolerance{this, "timeWindowTolerance",
                                             config().triplet.timeWindowTolerance,
                                             "extra hit-time window, ns"};
  ParameterRef<double> m_timeErrorFloor{this, "timeErrorFloor", config().triplet.timeErrorFloor,
                                        "min hit time error, ns"};
  ParameterRef<double> m_slowParticleMass{this, "slowParticleMass",
                                          config().triplet.slowParticleMass,
                                          "mass of the slowest seeded particle, GeV"};
  ParameterRef<std::size_t> m_fieldSamples{this, "fieldSamples", config().triplet.fieldSamples,
                                           "field samples for the effective field, 0 middle hit"};
  ParameterRef<double> m_positionError{this, "positionError", config().triplet.positionError,
                                       "seed position error, mm"};
  ParameterRef<double> m_angleError{this, "angleError", config().triplet.angleError,
                                    "seed angle error, rad"};
  ParameterRef<double> m_qOverPRelativeError{this, "qOverPRelativeError",
                                             config().triplet.qOverPRelativeError,
                                             "relative q/p error floor"};
  ParameterRef<double> m_timeError{this, "timeError", config().triplet.timeError,
                                   "seed time error, ns"};
  ParameterRef<bool> m_massFromTimeOfFlight{
      this, "massFromTimeOfFlight", config().triplet.massFromTimeOfFlight,
      "choose the proton or pion hypothesis of a seed from the times of its hits"};
  ParameterRef<double> m_pionSigmas{
      this, "pionSigmas", config().triplet.pionSigmas,
      "a seed is a pion if its 1/beta is this many errors below that of a proton"};
  ParameterRef<double> m_maxInverseBetaError{
      this, "maxInverseBetaError", config().triplet.maxInverseBetaError,
      "no hypothesis is written if the error on 1/beta exceeds this"};

  // Seeding2
  ParameterRef<std::string> m_fieldModel{this, "fieldModel", config().fieldModel,
                                         "scalar, effective or perMiddle"};
  ParameterRef<double> m_bFieldInZ{this, "bFieldInZ", config().bFieldInZ,
                                   "field of the fixed-frame finders, tesla"};
  ParameterRef<double> m_bFieldBinWidth{this, "bFieldBinWidth", config().bFieldBinWidth,
                                        "|B| bin width of the per-middle finders, tesla"};
  ParameterRef<double> m_rMin{this, "rMin", config().rMin, "mm"};
  ParameterRef<double> m_rMax{this, "rMax", config().rMax, "mm"};
  ParameterRef<double> m_zMin{this, "zMin", config().zMin, "mm"};
  ParameterRef<double> m_zMax{this, "zMax", config().zMax, "mm"};
  ParameterRef<double> m_rMinMiddle{this, "rMinMiddle", config().rMinMiddle, "mm"};
  ParameterRef<double> m_rMaxMiddle{this, "rMaxMiddle", config().rMaxMiddle, "mm"};
  ParameterRef<double> m_extraSigmaZ{this, "extraSigmaZ", config().extraSigmaZ,
                                     "extra hit error along z', mm"};
  ParameterRef<double> m_deltaRMinBottom{this, "deltaRMinBottom", config().deltaRMinBottom, "mm"};
  ParameterRef<double> m_deltaRMaxBottom{this, "deltaRMaxBottom", config().deltaRMaxBottom, "mm"};
  ParameterRef<double> m_deltaRMinTop{this, "deltaRMinTop", config().deltaRMinTop, "mm"};
  ParameterRef<double> m_deltaRMaxTop{this, "deltaRMaxTop", config().deltaRMaxTop, "mm"};
  ParameterRef<double> m_deltaZMin{this, "deltaZMin", config().deltaZMin, "mm"};
  ParameterRef<double> m_deltaZMax{this, "deltaZMax", config().deltaZMax, "mm"};
  ParameterRef<double> m_collisionRegionMin{this, "collisionRegionMin", config().collisionRegionMin,
                                            "mm"};
  ParameterRef<double> m_collisionRegionMax{this, "collisionRegionMax", config().collisionRegionMax,
                                            "mm"};
  ParameterRef<double> m_cotThetaMax{this, "cotThetaMax", config().cotThetaMax, ""};
  ParameterRef<double> m_impactMax{this, "impactMax", config().impactMax, "mm"};
  ParameterRef<double> m_minPt{this, "minPt", config().minPt, "GeV, transverse to z'"};
  ParameterRef<double> m_helixCutTolerance{this, "helixCutTolerance", config().helixCutTolerance,
                                           ""};
  ParameterRef<double> m_sigmaScattering{this, "sigmaScattering", config().sigmaScattering, ""};
  ParameterRef<double> m_radLengthPerSeed{this, "radLengthPerSeed", config().radLengthPerSeed, ""};
  ParameterRef<bool> m_sortedByCotTheta{this, "sortedByCotTheta", config().sortedByCotTheta, ""};
  ParameterRef<double> m_deltaInvHelixDiameter{this, "deltaInvHelixDiameter",
                                               config().deltaInvHelixDiameter, "1/mm"};
  ParameterRef<double> m_filterDeltaRMin{this, "filterDeltaRMin", config().filterDeltaRMin, "mm"};
  ParameterRef<double> m_compatSeedWeight{this, "compatSeedWeight", config().compatSeedWeight, ""};
  ParameterRef<std::size_t> m_compatSeedLimit{this, "compatSeedLimit", config().compatSeedLimit,
                                              ""};
  ParameterRef<double> m_impactWeightFactor{this, "impactWeightFactor", config().impactWeightFactor,
                                            ""};
  ParameterRef<double> m_zOriginWeightFactor{this, "zOriginWeightFactor",
                                             config().zOriginWeightFactor, ""};
  ParameterRef<std::size_t> m_maxSeedsPerSpM{this, "maxSeedsPerSpM", config().maxSeedsPerSpM, ""};
  ParameterRef<std::size_t> m_maxSeedsPerSpMConf{this, "maxSeedsPerSpMConf",
                                                 config().maxSeedsPerSpMConf, ""};
  ParameterRef<std::size_t> m_maxQualitySeedsPerSpMConf{this, "maxQualitySeedsPerSpMConf",
                                                        config().maxQualitySeedsPerSpMConf, ""};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void Process(int32_t /* run_number */, uint64_t /* event_number */) override {
    m_algo->process({m_hits_input()}, {m_seeds_output().get(), m_trackparams_output().get()});
  }
};

} // namespace eicrecon
