# ACTS geometry/material-map runtime contract

The ACTS geometry service normally remains permissive: it uses the `material-map`
constant declared by the loaded DD4hep geometry when present, otherwise the
historical generic fallback, and still permits `-Pacts:MaterialMap=...` overrides.

For detector/physics validation this can be too permissive. A map can load while
belonging to an older geometry, and a detector can have approach surfaces that
received no mapped material. The optional settings below make those assumptions
machine-checkable without hard-coding B0 into the generic ACTS service.

## Parameters

- `acts:RequireGeometryMaterialMap=false` (default): when true, the loaded DD4hep
  geometry must declare a `material-map`; the selected map and the geometry-
  declared map are hashed and must have identical SHA-256 content. A relocated
  byte-identical copy is therefore accepted, while an old/different override is
  rejected before event processing.
- `acts:ExpectedMaterialMapSHA256=` (default: geometry constant
  `material-map-sha256` when available): optional canonical content pin. The
  selected map must match it; with `RequireGeometryMaterialMap=true`, the
  geometry-declared map must match it as well.
- `acts:RequireMaterialCoverage=false` (default): enables strict final ACTS
  surface-material coverage checks.
- `acts:RequiredMaterialDetectorConstants=`: comma-separated DD4hep constant
  names identifying detector IDs whose ACTS **approach surfaces** must all have
  material. For example, the B0 validation chain uses
  `B0Tracker_Station_1_ID`. These constants are resolved to the same 8-bit
  detector identifier stored in `GeometryIdentifier::extra()` by the existing
  DD4hep-to-ACTS geometry-ID hook.

The service logs the geometry-declared path, final selected path, whether the
selection came from geometry/fallback/override, and hashes when a hash check is
requested. Existing jobs retain their previous behavior unless a strict option
is enabled.

## B0 validation

The B0 benchmark enables the strict geometry/map match and B0 approach-surface
coverage by default when it performs reconstruction. For an intentional expert
diagnostic with a mismatched map, set:

```sh
B0_STRICT_MATERIAL_MAP=0 \
src/benchmarks/reconstruction/b0_tracking/run_b0_benchmark.sh \
  --sim input.edm4hep.root --out results/diagnostic
```

This opt-out is recorded in provenance and should not be used for quoted B0
tracking performance.

To additionally pin a known validated artifact by content:

```sh
export ACTS_MATERIAL_MAP_SHA256=<64-hex-sha256>
```

For the realistic IP6 B0 geometry, the geometry-side material validation remains
the authority for producing/certifying a map. Runtime equality and surface
coverage are guards against accidental mismatch; they do **not** replace Geant4
material closure, bin-occupancy, navigation, or engineering-tolerance validation.

## Coverage semantics

Strict coverage is intentionally detector-scoped and based on final ACTS
`GeometryIdentifier`s, not detector-name string matching. Every final ACTS
surface with:

1. `approach() != 0`, and
2. `extra()` equal to one of the explicitly requested detector IDs

must have non-null `surfaceMaterial()`. Initialization fails if no approach
surfaces are found for a requested detector ID or if any selected surface lacks
material.

This is stricter than the historical diagnostic, which only reported a layer as
covered when *some* approach surface in the layer received material. The old
whole-geometry diagnostic remains for non-strict runs.

## Standalone contract test

The content-identity helper has no ACTS/DD4hep dependency:

```sh
cmake -S src/tests/material_map_contract -B build/material-map-contract \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/material-map-contract
ctest --test-dir build/material-map-contract --output-on-failure
```

The test checks standard SHA-256 vectors, relocated identical content, mismatched
content, missing geometry declarations, expected-hash normalization and failure
cases. Full service compilation and detector-backed coverage checks still need
the supported eic-shell stack.
