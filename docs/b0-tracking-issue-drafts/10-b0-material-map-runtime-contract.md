# B0: enforce a runtime contract between loaded geometry and ACTS material map

## Motivation

The realistic B0 ePIC configurations explicitly select a geometry-matched material map, but EICrecon still allows a different `acts:MaterialMap` to be supplied and falls back to the generic `calibrations/materials-map.cbor` when the geometry does not declare one.

`EpicJsonMaterialDecorator` currently reports layers for which no approach surface received material, but the check is diagnostic only: it logs a critical message rather than failing initialization. It also records success at the layer level when **any** approach surface has material, which is weaker than demonstrating complete B0 coverage.

This means a reconstruction job can start successfully even when the map/geometry pairing has not been validated for the realistic B0 telescope.

Companion geometry-side work lives in `tom-bleher/epic` and validates map provenance, surface coverage and material closure.

## Goal

Make the geometry/material-map relationship explicit, machine-checkable and optionally strict so that a production B0 benchmark cannot silently use an unvalidated map.

The core capability should remain detector-agnostic. B0-specific policy should be supplied by configuration/provenance rather than hard-coded detector-name checks in ACTS services.

## Proposed design

### 1. Expose declared vs selected material map

At ACTS service initialization, record separately:

- the material map declared by the loaded DD4hep geometry/configuration;
- the final material map selected after JANA parameter overrides;
- whether the generic fallback was used.

Log these values at info level and include them in benchmark provenance.

### 2. Add an opt-in strict geometry/map policy

Provide a reusable parameter such as `acts:RequireGeometryMaterialMap` (name illustrative) that can reject unintended map overrides for configurations that declare a validated map.

Path equality alone is not a sufficient long-term identity check because the same artifact may exist at another local path. Prefer a provenance identifier/content hash when available; a path check can be an initial conservative guard.

### 3. Strengthen material-coverage validation

When strict validation is requested, fail if required mapped approach surfaces are missing material rather than emitting only a log message.

The mechanism should support an explicit required-surface/detector set supplied by geometry or validation metadata so that central-detector surfaces that intentionally differ do not cause false failures.

### 4. Record immutable provenance

B0 benchmark output should record at least:

- ePIC geometry/configuration identity;
- EICrecon revision;
- selected material-map path and content hash;
- expected map/provenance identifier from the geometry-side validation;
- ACTS version;
- whether strict coverage checks passed.

## Acceptance criteria

- [ ] EICrecon distinguishes geometry-declared, user-selected and fallback material maps.
- [ ] A strict mode can reject a mismatched map before event processing.
- [ ] Strict coverage validation can fail on missing required B0 approach-surface material.
- [ ] The capability is generic/default-safe and is not keyed to string matching on `B0Tracker`.
- [ ] B0 benchmark manifests record the selected map hash and geometry identity.
- [ ] A deliberately mismatched/old map is rejected in the strict B0 validation configuration.

## Validation

Test at least:

1. matched realistic B0 geometry + validated map: pass;
2. explicit old/generic map with strict mode: fail with actionable diagnostics;
3. no geometry-declared map in ordinary configurations: preserve current default behavior unless strict policy is explicitly enabled;
4. local copy of the validated artifact: accepted when its content/provenance identity matches.

## Non-goals

- Preventing expert diagnostic runs with an intentionally mismatched map when strict mode is disabled.
- Encoding B0-specific logic inside the generic ACTS geometry service.
- Replacing the geometry-side material-mapping validation with a runtime filename check.
