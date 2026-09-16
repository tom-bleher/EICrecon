// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include <JANA/JApplication.h>
#include "B0ACLGAD_factory.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  app->Add(new JOmniFactoryGeneratorT<eicrecon::B0ACLGAD_factory>(
      "B0ACLGADResponse", {"EventHeader", "B0TrackerHits"},
      {"B0ACLGADRawHits", "B0ACLGADRawHitLinks", "B0ACLGADRawHitAssociations",
       "B0ACLGADChannelHits", "B0ACLGADMeasurements"}, {}, app));
}
}
