// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Tom Bleher

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <edm4eic/Cov6f.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/Vector3f.h>
#include <podio/RelationRange.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <initializer_list>

#include "algorithms/tracking/B0TrackSeeding.h"
#include "algorithms/tracking/B0TrackSeedingConfig.h"

// A uniform field across the beam, tilted as it is far off the midplane of the B0 dipole: a frame
// along any other direction would bend these tracks in the slope the Acts finders take as straight
TEST_CASE("B0TrackSeeding finds helices in a field off the vertical", "[B0TrackSeeding]") {
  const double tilt = 0.5; // rad
  const Acts::Vector3 fieldVector =
      1.3 * Acts::UnitConstants::T * Acts::Vector3(std::sin(tilt), std::cos(tilt), 0.);
  const Acts::ConstantBField field(fieldVector);
  const Acts::Vector3 along  = fieldVector.normalized();
  const Acts::Vector3 across = Acts::Vector3::UnitZ().cross(along); // the force for q > 0

  eicrecon::B0TrackSeeding algo("B0TrackSeeding");
  for (const double q : {+1., -1.}) {
    for (const double p : {5., 20., 100.}) {
      // Hits on four stations, on a track starting at z = 5800 mm along z and 10 mrad along the field
      edm4eic::TrackerHitCollection hits;
      const double radius = p / fieldVector.norm();
      for (const double z : {5850., 6272., 6601., 6880.}) {
        const double dz = z - 5800.;
        const Acts::Vector3 position =
            q * (radius - std::sqrt(radius * radius - dz * dz)) * across + 0.01 * dz * along +
            z * Acts::Vector3::UnitZ();
        hits.create(0ULL,
                    edm4hep::Vector3f{static_cast<float>(position.x()),
                                      static_cast<float>(position.y()),
                                      static_cast<float>(position.z())},
                    edm4eic::CovDiag3f{4.e-4F, 4.e-4F, 0.F}, 3.F, 0.F, 0.F, 0.F);
      }

      edm4eic::TrackSeedCollection seeds;
      edm4eic::TrackParametersCollection parameters;
      algo.findSeeds(hits, field, Acts::MagneticFieldContext{}, seeds, parameters);

      // Stations 1-2-3, 1-2-4, 1-3-4 and 2-3-4
      REQUIRE(seeds.size() == 4);
      for (const auto seed : seeds) {
        CHECK(seed.getParams().getQOverP() == Catch::Approx(q / p).epsilon(0.02));
        CHECK(seed.getParams().getTime() == Catch::Approx(3.));
        CHECK(std::sqrt(seed.getParams().getCovariance()(4, 4)) ==
              Catch::Approx(algo.getConfig().qOverPRelativeError / p).epsilon(0.02));
        CHECK(seed.getPerigee().z < seed.getHits()[0].getPosition().z);
      }
    }
  }
}
