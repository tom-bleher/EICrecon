// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include <JANA/JApplication.h>
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/tracking/B0TelescopeFactories.h"
#include "factories/tracking/AmbiguitySolver_factory.h"
#include <stdexcept>
#include <string>

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  using namespace eicrecon;
  std::string response = "effective";
  app->SetDefaultParameter("b0_telescope:response", response,
      "Input response: effective (legacy 70 um proxy) or aclgad (requires b0_aclgad plugin/physical profile)");
  if (response != "effective" && response != "aclgad")
    throw std::invalid_argument("b0_telescope:response must be effective or aclgad");
  const std::string measurements = response == "aclgad" ? "B0ACLGADMeasurements" : "B0TrackerMeasurements";
  const std::string associations = response == "aclgad" ? "B0ACLGADRawHitAssociations" : "B0TrackerRawHitAssociations";
  app->Add(new JOmniFactoryGeneratorT<B0TelescopeSeeding_factory>(
      "B0TelescopeSeeds", {"EventHeader", measurements},
      {"B0TelescopeSeeds", "B0TelescopeSeedParameters"}, {}, app));
  app->Add(new JOmniFactoryGeneratorT<B0TelescopeTruthSeeding_factory>(
      "B0TelescopeTruthSeeds", {"EventHeader", "B0TrackerHits", measurements, associations},
      {"B0TelescopeTruthSeeds", "B0TelescopeTruthSeedParameters"}, {}, app));
  for (const bool truth : {false, true}) {
    for (const bool direct : {false, true}) {
      const std::string stem = std::string("B0Telescope") + (truth ? "Truth" : "") + (direct ? "Direct" : "CKF");
      const std::string seeds = truth ? "B0TelescopeTruthSeeds" : "B0TelescopeSeeds";
      B0TelescopeTrackingConfig config;
      config.doFinalWeakPriorRefit = true;
      config.fitSeedMeasurementsDirectly = direct;
      config.useSeedParticleHypothesis = truth;
      app->Add(new JOmniFactoryGeneratorT<B0TelescopeTracking_factory>(
          stem, {"EventHeader", seeds, measurements},
          {stem + "StatesUnfiltered", stem + "ActsTracksUnfiltered"}, config, app));
      app->Add(new JOmniFactoryGeneratorT<AmbiguitySolver_factory>(
          stem + "Ambiguity", {stem + "StatesUnfiltered", stem + "ActsTracksUnfiltered"},
          {stem + "States", stem + "ActsTracks"}, {.n_measurements_min = 3}, app));
      for (const bool unfiltered : {false, true}) {
        const std::string suffix = unfiltered ? "Unfiltered" : "";
        app->Add(new JOmniFactoryGeneratorT<B0TelescopeTrackExport_factory>(
            stem + "Tracks" + suffix,
            {measurements, seeds, stem + "States" + suffix,
             stem + "ActsTracks" + suffix, associations},
            {stem + "Trajectories" + suffix, stem + "Parameters" + suffix, stem + "Tracks" + suffix,
             stem + "Links" + suffix, stem + "Associations" + suffix}, {}, app));
      }
    }
  }
}
}
