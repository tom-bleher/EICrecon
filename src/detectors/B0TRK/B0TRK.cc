// Copyright 2022, Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <string>
#include <vector>
#include <edm4eic/unit_system.h>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // Digitization
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "B0TrackerRawHits", {"EventHeader", "B0TrackerHits"},
      {"B0TrackerRawHits", "B0TrackerRawHitLinks", "B0TrackerRawHitAssociations"},
      {
          // The realistic B0 sensor is 50 um thick (B0TrackerSensorThickness), so
          // a minimum-ionising proton deposits only ~14 keV most probably and 10 %
          // of hits fall below the 10 keV inherited from the far-forward Roman
          // pots. Use the central silicon value: on a 41 GeV proton gun this
          // takes the stub-seeded efficiency from 79 % to 85 % and the fraction
          // of four-measurement tracks from 72 % to 91 %, at unchanged fake rate.
          .threshold      = 0.54 * dd4hep::keV,
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  // Convert raw digitized hits into hits with geometry info (ready for tracking)
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "B0TrackerRecHits", {"B0TrackerRawHits"}, {"B0TrackerRecHits"},
      {
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));
}
} // extern "C"
