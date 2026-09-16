# Opt-in B0 telescope reconstruction

This parallel path does not replace `B0TrackerSeeds`, modify central tracking, change detector pitch, or merge local-plane tracks into `CombinedTracks`. Load `b0_telescope` and request its output collections; standard digitization/tracking still supplies `B0TrackerMeasurements`.

## Four comparison chains

| Stem | Finding | Final fit |
|---|---|---|
| `B0TelescopeCKF` | local data seeds + ACTS CKF | independent weak-prior KF |
| `B0TelescopeDirect` | local data candidates + sensor extension | independent weak-prior KF |
| `B0TelescopeTruthCKF` | simulated-crossing local seeds + CKF | independent weak-prior KF |
| `B0TelescopeTruthDirect` | truth-associated sensor measurements | independent weak-prior KF |

Each stem produces `Tracks`, `Parameters`, `Trajectories`, `Links`, `Associations`; append `Unfiltered` for pre-ambiguity outputs. ACTS containers remain transient. Data seeds are `B0TelescopeSeeds`/`B0TelescopeSeedParameters`; truth seeds are `B0TelescopeTruthSeeds`/`B0TelescopeTruthSeedParameters`. Request an explicit collection list when timing performance: four chains intentionally cost more than one.

## Reference surfaces

`TrackParameters.surface` is the actual sensitive-plane geometry ID. The six bound parameters are `(loc0,loc1,phi,theta,q/p,time)` at that plane. The dedicated exporter preserves it and transforms the full covariance into Cartesian coordinates.

These separately named `Track` collections contain B0-plane position/momentum, **not vertex momentum**. Do not feed them unchanged to vertex-based analyses. Exact-zero q/p retains meaningful bound parameters but is marked `Track.type=-1` with a nonphysical zero Cartesian momentum placeholder; use the bound parameters instead.

`promptSelection=false` is the default local/displaced mode and does not require propagation to the IP. Setting it true applies d0/z0 cuts after full-field backward propagation. This is a post-fit **selection**, not a calibrated external beam-spot likelihood. No beam-spot-constrained fit is claimed.

## Independent fitting

Finding uses the hit-derived seed and explicit search floors. Final fitting uses its transported mean but a fresh absolute diagonal covariance, **not a scaled seed posterior**. It assimilates the selected measurements once with ACTS field/material transport. Changed measurement sets and fit failures are counted and rejected, with no silent fallback to seed-influenced covariance.

The reusable `doFinalWeakPriorRefit` option defaults off; this opt-in plugin explicitly enables it. Direct fits always use fresh initialization. No behavior is selected implicitly by a B0 station-count flag.

Scan `weakPriorScale=1,10,100` and starting means before quoting covariance coverage. Finite weak priors still need convergence checks. Synthetic information-accounting tests are not proof of detector-level calibrated pulls. Early predicted innovations may retain initialization dependence.

## Configuration

JANA parameters use `<plugin>:<factory-tag>:<parameter>`, for example:

```sh
-Pb0_telescope:B0TelescopeSeeds:fieldIntegral=true
-Pb0_telescope:B0TelescopeSeeds:useTiming=true
-Pb0_telescope:B0TelescopeSeeds:stationMapFile=/absolute/path/stations.json
-Pb0_telescope:B0TelescopeSeeds:diagnosticsFile=/absolute/path/seeds.jsonl
-Pb0_telescope:B0TelescopeTruthSeeds:referenceFile=/absolute/path/truth.jsonl
-Pb0_telescope:B0TelescopeCKF:diagnosticsFile=/absolute/path/ckf.jsonl
-Pb0_telescope:B0TelescopeCKF:weakPriorScale=10
-Pb0_telescope:B0TelescopeCKF:promptSelection=true
```

Check emitted names with the executable's parameter listing. Each writer needs a distinct filename and an existing parent directory. JSONL event writes are mutex-protected; event ordering is not guaranteed. Do not share filenames across processes.

Seeder parameters include `minStations`, `maxTrials`, `maxCandidates`, `minMomentum` [GeV], `maxSlope`, `fieldBoundTesla`, `roadWidthMm`, `maxChi2PerDof`, `extensionChi2`, `curvatureSignificance`, `proposalMaxMomentum` [GeV], `forceCharge`, `bothCharges`. Search floors `seedPositionSigma` [mm], `seedAngularSigma` [slope], `seedQOverPSigma` [1/GeV] are not detector resolutions. See [NUMERICS.md](NUMERICS.md) for approximation and budget definitions.

`particleHypothesisPdg` is a fit-mass hypothesis, not reconstructed PID. Realistic chains default to protons; configure other samples explicitly. Truth chains use seed PDG if `useSeedParticleHypothesis=true`; disable it for a strictly identical-mass comparison. Timing has a separate `massGeV` and model allowance. The field-integral and timing options default off.

Absolute weak variances: `weakPriorLoc0Variance`, `weakPriorLoc1Variance` [mm^2], `weakPriorPhiVariance`, `weakPriorThetaVariance` [rad^2], `weakPriorQOverPVariance` [GeV^-2], `weakPriorTimeVariance` [ns^2]. They are numerical initialization settings, not measured detector properties.

## Truth reference and diagnostics

The truth seeder uses actual B0 crossings and the same surface/station map as the realistic seeder and fitter. It does not reuse central production-vertex or eta cuts. Supported charged species are e, mu, pi, K and p. Local simulated-hit momentum/time and deterministic positive momentum smearing initialize truth seeds. Truth relations never enter the realistic finder.

The truth JSONL contains particle identities, station counts, unsmeared local states, seed identities and unmapped/invalid counters. A tangent projects a simulated hit only onto its own thin sensor plane; `maxProjectionMm` bounds that approximation (default 0.1 mm). It is not a long-range straight-line extrapolation. Check this approximation against the precision needed for the analysis.

Truth-associated direct fits use strict majority purity and at most one measurement per surface. Truth seeds can exist despite digitization losses. Reconstructibility therefore remains a crossing-based denominator, independent of seed/fit success. Record supported species, momentum thresholds and missing truth-reference states separately.

The station manifest uses string uint64 surface IDs, physical station, ACTS layer, centre and normal. Seeder counters include endpoint/candidate work and truncation. Tracking JSONL records stage failures and predicted residuals/innovation covariance by measurement. These are separate from the legacy aggregate counters.

Before promotion, require supported-stack builds, detector-backed smoke runs of all four chains, Geant4/ACTS surface/material comparison, prior/charge scans, displaced and mixed-species samples, background overlays, close pairs, one-sided coverage, module edges and missing stations. Measure efficiency, fakes, duplicates and runtime tails; no gain is claimed merely because a new algorithm exists.

References: [ACTS fitting](https://acts.readthedocs.io/en/latest/core/reconstruction/track_fitting.html), [upstream B0 discussion](https://github.com/eic/EICrecon/issues/2930), and fork PR #3's review drafts.
