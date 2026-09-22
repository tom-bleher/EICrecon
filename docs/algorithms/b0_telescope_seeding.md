# Experimental B0 telescope seeding

This branch adds a small, opt-in specialization of the modular ACTS seeding
framework. It is intended for low-occupancy commissioning of the realistic B0
tracker, not as a validated tracking-efficiency result.

## Architecture

`B0TrackerRecHits` -> four physical station groups -> `Acts::TripletSeeder`
-> ACTS three-point parameter estimator -> ACTS backward field propagation
-> the existing B0 CKF, ambiguity solver and EDM conversion.

The telescope policies derive from `Acts::DoubletSeedFinder`,
`Acts::TripletSeedFinder`, and `Acts::ITripletSeedFilter`.
`createSeedsFromGroups` receives the station triples (1,2,3), (1,2,4),
(1,3,4), and (2,3,4). Radial-order and cot(theta)-sorting fast paths are
explicitly disabled. The built-in solenoidal compatibility cuts and filter
are not used. No cylindrical grid or KD-tree is needed for this first version.

The detector rotation comes from `B0Tracker_rotation`. Readout layer IDs are
decoded with DD4hep, not bit shifts or clustering of whichever hits happened to
arrive. The default explicit contract is layers (1,2), (3,4), (5,6), (7,8) for
four physical stations. `layersPerStation=1` selects the ideal four-layer
readout. Keep the individual front/back measurements; they are not averaged.

Doublets use a forward separation and loose Cartesian slopes. Triplets use
Cartesian sagitta windows, wider in the bending plane. There is no interaction-
point or production-vertex cut. The stored seed quality is minus the normalized
sagitta score, not a fitted chi-square. Exact duplicate triplets are removed;
nearby front/back alternatives are left for the existing ambiguity solver.

Parameter estimation is separate from seed finding. ACTS' conformal three-point
estimator receives the actual three-dimensional field vector sampled at the
middle hit. It estimates a state at a plane through the first hit. ACTS then
propagates that state and its covariance backward through the field provider to
the origin-centred perigee required by the current EICrecon CKF. This reference
surface does NOT constrain the vertex to the origin. In particular, a tangent
measured in the magnet is not simply relabelled as a direction at the IP.

## Build and run

Use an eic-shell installation with ACTS >= 45.3. The adapter supports the
Seeding2 names in ACTS 45.3/46 and Seeding names in ACTS 47, and detects the
mutable-view versus returned-cursor interface signatures. Older ACTS builds
keep the existing seeder; requesting telescope mode fails explicitly.

Build and source the realistic detector fork first. Reconstruct using the same
geometry and a material map generated for that geometry. Then, in this checkout:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build --target install -j4
source install/bin/eicrecon-this.sh

# DETECTOR_PATH must refer to the built realistic detector, not the stock install.
export DETECTOR_CONFIG=epic_ip6_extended

eicrecon \
  -Ptracking:B0Seeding=telescope \
  -Ptracking:B0TrackerSeeds:layersPerStation=2 \
  -Pacts:LayerEnvelopeZ=5 \
  -Ppodio:output_file=b0-telescope.edm4eic.root \
  input.edm4hep.root
```

The output names remain `B0TrackerSeeds`, `B0TrackerSeedParameters`,
`B0TrackerCKFTracks`, etc. The normal output configuration preserves their
relations; when trimming output collections, also write the referenced hit,
parameter, measurement and track collections. `tracking:B0Seeding=standard`
(the default) retains the existing seed finder. The truth-seeded chain is
unchanged and can be used as a control on the same events.

Telescope mode sets only the B0 hit-seeded CKF and ambiguity solver minimum
measurement counts to three. It does not widen the CKF chi-square cut, change
its particle hypothesis, add a station selector inside CKF, or modify central
tracking. Three fitted measurements are not necessarily three physical
stations: final track-level station coverage still needs to be checked in
validation. This branch retains the independent geometry-clearance hook from
the realistic-geometry bring-up; its longitudinal envelope defaults to 5 mm.

## Tests

The geometry kernel can be tested without JANA, ROOT, DD4hep or ACTS:

```sh
cmake -S src/tests/b0_telescope_math_test -B build-b0-math
cmake --build build-b0-math
ctest --test-dir build-b0-math --output-on-failure
```

The ACTS integration cases are in the usual algorithms test target:

```sh
cmake --build build --target algorithms_test -j4
./build/src/tests/algorithms_test/algorithms_test '[b0telescope]'
```

They exercise all four station triples, missing stations, front/back grouping,
decreasing radius, displaced tracks, rejection of crossed track combinations,
event-local state, and the parameter estimator with both charge signs in a
uniform transverse field. The standalone geometry tests were run with GCC 14,
AddressSanitizer and UndefinedBehaviorSanitizer: 2,048 checks passed. The full
ACTS/EICrecon build and detector reconstruction were not available in the
implementation workspace and have not been certified by those numerical tests.

## Scope and next validation

The defaults are deliberately loose, uncalibrated commissioning settings.
`maxCombinations` rejects an entire busy event with a warning rather than taking
an input-order-biased prefix. Seed caps use deterministic ranking; a global cap
emits a warning. Account for rejected events and capped seeds in efficiencies.
The per-middle cap is per station triple, as its parameter description states.

The local helix estimate assumes an approximately uniform field over the seed;
subsequent propagation uses the full field. Initial covariance is a configurable
diagonal approximation at the first-hit plane, not a tuned pull model. Seed
transport is field-only; material effects and mass dependence are handled by the
existing pion-hypothesis CKF. No MC information is used to form telescope seeds.
Zero-field, degenerate, unresolved-curvature, failed-transport and invalid-
covariance estimates are skipped and counted. There is no invented fixed
momentum or charge-sign fallback.

Validate clean positive/negative pion and proton samples against truth-seeded
CKF on the same geometry first. Record each stage separately: reconstructible
physical stations, compatible triplets, successful parameter estimates,
successful CKF fits, reference extrapolation, and ambiguity selection. Then test
missing stations, mixed tracks, displaced production and background occupancy.
Compare momentum/direction residuals, fake/duplicate rates, wall time and output
identity at one versus multiple threads. Backward transport to a perigee and the
unchanged downstream extrapolator can still reject displaced or unusual tracks;
absence of an IP cut is not a claim of complete secondary-particle acceptance.
