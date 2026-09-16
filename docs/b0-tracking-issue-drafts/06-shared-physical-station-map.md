# B0 tracking: use one shared physical-station map in seeding, CKF selection, and benchmarks

## Motivation

The realistic B0 geometry has eight ACTS front/back layers but four physical stations. Correct reconstruction therefore depends on consistently answering the question: **which ACTS surface/hit belongs to which physical station?**

The current branch implements this concept in multiple places:

- `B0TrackerStubSeeder` discovers B0 surfaces, clusters them in ion-frame `z`, and maps volume IDs to station indices;
- `CKFTracking` separately scans B0 sensitive surfaces and builds its own geometry-ID-to-station mapping for `numB0StationsMin`;
- the offline benchmark independently clusters station positions/hits with the same nominal 50 mm gap.

These implementations currently aim to encode the same physical concept, but duplicated discovery/selection rules are a correctness and maintenance risk as the geometry changes.

Relevant code:

- `src/algorithms/tracking/B0TrackerStubSeeder.cc`
- `src/algorithms/tracking/CKFTracking.cc`
- `src/algorithms/tracking/CKFTrackingConfig.h`
- `src/benchmarks/reconstruction/b0_tracking/analyze_b0_tracking.py`

Geometry context: four stations separated by roughly 225–276 mm, with front/back sensor planes only a few mm apart.

## Goal

Create a single authoritative B0 physical-station mapping abstraction and reuse it everywhere that station multiplicity or station identity matters.

## Proposed design

Introduce a small geometry utility, for example:

```text
B0StationMap
  - discover B0 sensitive ACTS surfaces from the geometry
  - compute ion-frame z for each surface
  - cluster front/back surfaces into physical stations
  - map GeometryIdentifier / volume ID / surface to station index
  - expose station z intervals/means
```

The exact API can differ, but the key property is that the seeder and CKF should not have independent definitions of the same mapping.

### Discovery robustness

Prefer stable detector/geometry identifiers over string path matching where available. If path/name matching remains necessary, isolate it inside the shared helper and test it against both:

- the older four-layer/one-plane-per-station geometry;
- the realistic eight-layer front/back geometry.

### Configuration

Keep the clustering gap configurable, but define it once. Validate that the chosen value is safely between:

- maximum within-station front/back separation;
- minimum inter-station separation.

Emit a clear initialization error/warning if the discovered geometry does not produce the expected station topology for a B0 chain that requires it.

## Acceptance criteria

- [ ] Seeder and CKF use the same in-process station-map implementation.
- [ ] The station gap/configuration is not duplicated across independent algorithms without a shared source.
- [ ] Unit tests cover 4-plane and 8-plane geometries using representative z values.
- [ ] Front/back surfaces of one station map to the same physical station.
- [ ] Neighboring physical stations never merge for the current geometry.
- [ ] Unknown/unmapped B0 surfaces are observable through counters/logging and do not silently alter station multiplicity.
- [ ] Benchmark provenance records the station-map definition used by reconstruction.
- [ ] Where practical, export station identity or enough metadata for the offline benchmark to use the reconstruction definition rather than reimplementing it.

## Non-goals

- Changing the detector geometry.
- Requiring exactly two measurements per station.
- Treating front and back as statistically identical measurements; this issue only centralizes physical-station identity.
