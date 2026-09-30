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
#include "factories/digi/SiliconChargeSharing_factory.h"
#include "factories/digi/SiliconTrackerDigi_factory.h"
#include "factories/tracking/SiliconPadClustering_factory.h"
#include "factories/tracking/TrackerHitReconstruction_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // AC-LGAD pads share each deposit with their neighbours (electrode size of the BNL sensor
  // tested on EICROC0: 100 um metal on 500 um pitch)
  app->Add(new JOmniFactoryGeneratorT<SiliconChargeSharing_factory>(
      "B0TrackerSharedHits", {"B0TrackerHits"}, {"B0TrackerSharedHits"},
      {
          .min_edep    = 0,
          .readout     = "B0TrackerHits",
          .model       = SiliconChargeSharingConfig::EModel::trapezoid,
          .electrode_x = 0.1 * dd4hep::mm,
          .electrode_y = 0.1 * dd4hep::mm,
      },
      app));

  // Digitization: noise and threshold act on the summed pad signal, on the deposited-energy
  // scale (MIP MPV 12.3 keV in 50 um Si). Noise is 1.5% of the MPV and the pad threshold 4 sigma;
  // the relative term is an effective 2-pad resolution floor. Tuned so that 150 um electrodes
  // reproduce HPK 500 um-pitch pixel test-beam results (Dutta et al., NIM A (2025) 170224)
  app->Add(new JOmniFactoryGeneratorT<SiliconTrackerDigi_factory>(
      "B0TrackerRawHits", {"EventHeader", "B0TrackerSharedHits"},
      {"B0TrackerRawHits", "B0TrackerRawHitLinks", "B0TrackerRawHitAssociations"},
      {
          .threshold          = 0.74 * dd4hep::keV,
          .timeResolution     = 30 * edm4eic::unit::ps,
          .thresholdOnCellSum = true,
          .noise              = 0.18 * dd4hep::keV,
          .relativeNoise      = 0.25,
      },
      app));

  // Convert raw digitized hits into pad hits with geometry info
  app->Add(new JOmniFactoryGeneratorT<TrackerHitReconstruction_factory>(
      "B0TrackerRecHits", {"B0TrackerRawHits"}, {"B0TrackerRecHits"},
      {
          .timeResolution = 30 * edm4eic::unit::ps,
      },
      app));

  // Merge the pads of each particle crossing into one cluster hit for tracking. Resolutions are
  // the core residual widths of the digitization above for 100 um electrodes (41 GeV protons)
  app->Add(new JOmniFactoryGeneratorT<SiliconPadClustering_factory>(
      "B0TrackerClusterHits", {"B0TrackerRecHits"}, {"B0TrackerClusterHits"},
      {
          .readout               = "B0TrackerHits",
          .electrode_x           = 0.1 * dd4hep::mm,
          .electrode_y           = 0.1 * dd4hep::mm,
          .single_pad_resolution = 0.052 * dd4hep::mm,
          .shared_pad_resolution = 0.030 * dd4hep::mm,
      },
      app));
}
} // extern "C"
