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

#include "algorithms/digi/SiliconChargeSharingConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/digi/SiliconChargeSharing_factory.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/tracking/LGADHitClustering_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // Share each deposit among the 500 um AC-LGAD pads, with the sharing width of the TOF sensors.
  // The cut is on each pad's share of one deposit. TOF uses 6 keV, which is too high here: a MIP
  // deposits only about 12 keV in the 50 um B0 sensors, so most pads would be cut. Shares of at
  // least 0.5 keV also keep a non-zero charge in the integer-keV raw hits
  app->Add(new JOmniFactoryGeneratorT<SiliconChargeSharing_factory>(
      "B0TrackerSharedHits", {"B0TrackerHits"}, {"B0TrackerSharedHits"},
      {
          .sigma_mode     = SiliconChargeSharingConfig::ESigmaMode::rel,
          .sigma_sharingx = 0.5,
          .sigma_sharingy = 0.5,
          .min_edep       = 0.5 * dd4hep::keV,
          .readout        = "B0TrackerHits",
      },
      app));

  // Digitization. The threshold acts on each shared piece before summing, so it stays at zero
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "B0TrackerRawHits", {"EventHeader", "B0TrackerSharedHits"},
      {"B0TrackerRawHits", "B0TrackerRawHitLinks", "B0TrackerRawHitAssociations"},
      {
          .threshold      = 0.0,
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  // Convert raw digitized hits into pad hits with geometry info
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "B0TrackerRecHits", {"B0TrackerRawHits"}, {"B0TrackerRecHits"},
      {
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  // Merge touching pads into one charge-weighted measurement per particle crossing
  app->Add(new JOmniFactoryGeneratorT<LGADHitClustering_factory>(
      "B0TrackerMeasurements", {"B0TrackerRecHits"}, {"B0TrackerMeasurements"},
      {
          .readout = "B0TrackerHits",
          .useAve  = true,
      },
      app));
}
} // extern "C"
