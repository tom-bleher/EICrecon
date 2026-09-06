// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#pragma once

#include "algorithms/tracking/B0TrackerActsSeeding.h"
#include "extensions/jana/JOmniFactory.h"

namespace eicrecon {
class B0TrackerActsSeeding_factory
    : public JOmniFactory<B0TrackerActsSeeding_factory, B0TrackerActsSeedingConfig> {
private:
  PodioInput<edm4eic::TrackSeed> m_candidates{this};
  PodioInput<edm4eic::Measurement2D> m_measurements{this};
  PodioOutput<edm4eic::TrackSeed> m_seeds{this};
  PodioOutput<edm4eic::TrackParameters> m_parameters{this};
  ParameterRef<bool> m_useMultipoint{this, "useMultipoint", config().useMultipoint,
                                     "use ACTS multipoint; false copies the original stub seeds"};
  ParameterRef<double> m_inflation{this, "covarianceInflation", config().covarianceInflation,
                                   "inflate the correlated seed prior before CKF updates"};
  ParameterRef<bool> m_diagonal{this, "diagonalCovariance", config().diagonalCovariance,
                                "discard seed off-diagonal covariance before inflation"};
  ParameterRef<double> m_scattering{this, "scatteringScale", config().scatteringScale,
                                    "first-surface angular model sigma coefficient [GeV rad]"};
  ParameterRef<double> m_qop{this, "qOverPRelativeUncertainty", config().qOverPRelativeUncertainty,
                             "relative q/p model uncertainty"};
  ParameterRef<double> m_time{this, "timeVariance", config().timeVariance, "time variance [ns^2]"};
  ParameterRef<double> m_field{this, "minAbsField", config().minAbsField, "minimum field norm [T]"};
  ParameterRef<double> m_significance{this, "minCurvatureSignificance",
                                      config().minCurvatureSignificance,
                                      "minimum q/p significance to resolve the charge"};
  ParameterRef<int> m_charge{this, "charge", config().charge, "0 infer charge; +/-1 force charge"};
  ParameterRef<bool> m_both{this, "testBothCharges", config().testBothCharges,
                            "emit both charge hypotheses"};
  ParameterRef<unsigned int> m_maxSeeds{this, "maxSeeds", config().maxSeeds,
                                        "maximum emitted seeds"};
  std::unique_ptr<B0TrackerActsSeeding> m_algorithm;

public:
  void Configure() {
    m_algorithm = std::make_unique<B0TrackerActsSeeding>(GetPrefix());
    m_algorithm->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algorithm->applyConfig(config());
    m_algorithm->init();
  }
  void Process(int32_t /* run_number */, uint64_t /* event_number */) override {
    m_algorithm->process({m_candidates(), m_measurements()},
                         {m_seeds().get(), m_parameters().get()});
  }
};
} // namespace eicrecon
