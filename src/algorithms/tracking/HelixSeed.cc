// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include "HelixSeed.h"

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Utilities/MathHelpers.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>
#include <edm4eic/Cov6f.h>
#include <edm4hep/Vector2f.h>
#include <edm4hep/Vector3f.h>
#include <cmath>

namespace eicrecon {

bool addHelixSeed(const HelixSeedPrior& prior, const std::array<std::uint32_t, 3>& triplet,
                  float quality, const edm4eic::TrackerHitCollection& hits,
                  std::span<const Acts::Vector3> positions, const Acts::Vector3& field,
                  edm4eic::TrackSeedCollection& seeds,
                  edm4eic::TrackParametersCollection& parameters) {
  const Acts::Vector3& first = positions[triplet[0]];
  // Acts 45.2 added the start time; the hit time is set below
#if Acts_VERSION_MAJOR > 45 || (Acts_VERSION_MAJOR == 45 && Acts_VERSION_MINOR >= 2)
  const Acts::FreeVector helix = Acts::estimateTrackParamsFromSeed(first, 0., positions[triplet[1]],
                                                                   positions[triplet[2]], field);
#else
  const Acts::FreeVector helix =
      Acts::estimateTrackParamsFromSeed(first, positions[triplet[1]], positions[triplet[2]], field);
#endif
  const double qOverP = helix[Acts::eFreeQOverP] * Acts::UnitConstants::GeV / edm4eic::unit::GeV;
  if (!helix.allFinite() || qOverP == 0.) {
    return false;
  }

  // CKFTracking rebuilds the perigee along global z from the seed's point, so the seed is at
  // loc0 = loc1 = 0. Forward tracks run nearly parallel to it, so the loc1 and phi errors are
  // scaled up.
  const Acts::Vector3 direction = helix.segment<3>(Acts::eFreeDir0);
  const double theta            = Acts::VectorHelpers::theta(direction);
  const double sinTheta         = std::sin(theta);
  // The perigee is not defined for a track along its axis
  if (!(sinTheta > 0.)) {
    return false;
  }

  edm4eic::Cov6f cov;
  cov(0, 0) = static_cast<float>(Acts::square(prior.positionError));
  cov(1, 1) = static_cast<float>(Acts::square(prior.positionError / sinTheta));
  cov(2, 2) = static_cast<float>(Acts::square(prior.angleError / sinTheta));
  cov(3, 3) = static_cast<float>(Acts::square(prior.angleError));
  cov(4, 4) = static_cast<float>(Acts::square(prior.qOverPRelativeError * qOverP));
  cov(5, 5) = static_cast<float>(Acts::square(prior.timeError));

  auto params = parameters.create();
  params.setType(-1); // seed
  params.setLoc({0.F, 0.F});
  params.setPhi(static_cast<float>(Acts::VectorHelpers::phi(direction)));
  params.setTheta(static_cast<float>(theta));
  params.setQOverP(static_cast<float>(qOverP));
  params.setTime(hits[triplet[0]].getTime());
  params.setCovariance(cov);

  const Acts::Vector3 anchor =
      (first - prior.anchorDistance * direction) / Acts::UnitConstants::mm * edm4eic::unit::mm;
  auto seed = seeds.create();
  seed.setPerigee({static_cast<float>(anchor.x()), static_cast<float>(anchor.y()),
                   static_cast<float>(anchor.z())});
  seed.setQuality(quality);
  seed.setParams(params);
  for (const auto i : triplet) {
    seed.addToHits(hits[i]);
  }
  return true;
}

} // namespace eicrecon
