# B0 tracking benchmark

Pinned evaluation for the B0 stub-seeded CKF chain. This is **prompt
forward-proton reconstruction for solenoid-off `epic_ip6_extended`**, not a
general B0 tracker.

Reconstruction behaviour is unchanged by this benchmark. Failure counters are
additive only.

## Scope relative to central tracking

The B0 chain is a separate factory graph (`B0TrackerSeeds` → `B0TrackerCKF*` →
combined collections). Central seeding and CKF stay on their own instances.

| Work | Where | Central tracking |
|---|---|---|
| PR1 counters + this benchmark | B0 seeder; observe-only increments in shared `CKFTracking` | Unchanged |
| PR2 seed families | `B0TrackerStubSeeder` only | Unchanged |
| PR3 independent final refit | Shared `CKFTracking` **capability**, B0-only behavior | Must remain unchanged |
| PR4 cut scans | B0 factory config (`chi2CutOff`, seeder residuals) | Unchanged |
| PR5 occupancy tests | Benchmarks | Unchanged |

PR1 may increment counters inside `CKFTracking.cc` only when
`numB0StationsMin > 0`. That gate is observability, not a physics change:
central CKF keeps `numB0StationsMin = 0` and the accept/reject branches are
the same either way.

**PR3 must not become a global CKF change.** Do not infer “this is B0” from
`numB0StationsMin` in order to change fitted parameters. Add an explicit
config flag, default **off**:

```text
CKFTrackingConfig::doFinalWeakPriorRefit = false   // central, and the default
```

Enable it only on the B0 factory instances in `src/global/tracking/tracking.cc`.
Central remains `seed → CKF → smooth/extrapolate`. B0 becomes
`seed → CKF finding → weak-prior refit of selected measurements → smooth/extrapolate`.
That is a shared capability with a B0-only behavior change, not a change to
ACTS reconstruction everywhere.

## One command

From an eic-shell with this EICrecon install sourced:

```bash
src/benchmarks/reconstruction/b0_tracking/run_b0_benchmark.sh \
  --reco reco.edm4eic.root \
  --out results/b0_tracking
```

To reconstruct first (still no seeder/CKF parameter changes):

```bash
export B0_TRACKING_COUNTERS_FILE=/tmp/b0_counters.json   # optional; the run script sets this
src/benchmarks/reconstruction/b0_tracking/run_b0_benchmark.sh \
  --sim sim.edm4hep.root \
  --out results/b0_tracking
```

Outputs:

- `report.md` — efficiencies, residual summaries, failure counts
- `summary.json`
- `provenance.json` — SHAs, geometry, material map hash, matching definitions
- `b0_tracking_validation.pdf`
- `counters.json` when `B0_TRACKING_COUNTERS_FILE` is set (also dumped at eicrecon exit)

Self-test of the geometry-free helpers:

```bash
src/benchmarks/reconstruction/b0_tracking/run_b0_benchmark.sh --self-test
```

## Frozen baseline

Fill every field in `provenance.json` from the machine that produces a result.
The intended pin is:

| Item | Value |
|---|---|
| EICrecon | `feat/b0-tracking` at the SHA recorded by the run script |
| Geometry | local `epic` install, `DETECTOR_CONFIG=epic_ip6_extended` |
| Compact | `$DETECTOR_PATH/$DETECTOR_CONFIG.xml` |
| Field | generated `epic_ip6_extended.xml` includes `compact/fields/beamline_18x275.xml`. `B0PF_Bmax` and `B0PF_GradientMax` match `beamline_5x41.xml`; downstream magnets do not. There is no solenoid in this configuration. |
| Material map | `calibrations/materials-map-ip6-extended-cbc0605e892b.cbor` (override with `ACTS_MATERIAL_MAP`) |
| Sample | prompt protons, 5×41, \(p=8\)–\(41\) GeV, ion-frame \(\theta=4\)–\(22\) mrad |
| Reference surface | origin perigee `(0,0,0)` |
| Physical station | ion-frame \(z\) single-linkage clustering, gap \(50\,\mathrm{mm}\) |
| Reconstructible | `generatorStatus==1`, \(\lvert\mathrm{PDG}\rvert=2212\), hits in \(\ge 3\) physical B0 stations |
| Seed match | majority of seed rec hits → raw-hit associations → `MCParticle` |
| Track match | `B0TrackerCKFTrackAssociations` with weight \(\ge 0.5\) |

### Suggested gun (npsim)

```bash
npsim --compactFile $DETECTOR_PATH/$DETECTOR_CONFIG.xml \
  --numberOfEvents 3000 \
  --enableGun \
  --gun.particle proton \
  --gun.distribution uniform \
  --gun.momentumMin "8*GeV" --gun.momentumMax "41*GeV" \
  --gun.thetaMin "4*mrad" --gun.thetaMax "22*mrad" \
  --outputFile b0_protons.edm4hep.root
```

Use the ion-beam direction appropriate to IP6 (25 mrad crossing). Do not treat a
lab-frame \(\theta\) gun as identical to the ion-frame 4–22 mrad sample quoted
below until that is verified in `provenance.json`.

## Efficiency definitions

A truth proton is reconstructible independently of whether reconstruction
succeeded.

\[
\epsilon_{\mathrm{seed}}
=
\frac{N(\text{reconstructible with a compatible seed})}
     {N(\text{reconstructible})}
\]

\[
\epsilon_{\mathrm{track}\mid\mathrm{seed}}
=
\frac{N(\text{seeded reconstructible with a matched stub-seeded track})}
     {N(\text{seeded reconstructible})}
\]

\[
\epsilon_{\mathrm{total}}
=
\frac{N(\text{reconstructible with a matched stub-seeded track})}
     {N(\text{reconstructible})}
\]

Fake rate: reconstructed tracks with no reconstructible match, over all
reconstructed B0 tracks. Duplicate rate: extra matched tracks beyond one per
truth particle, over reconstructible truth.

Final-track pulls and

\[
D^2=(\hat a-a_{\mathrm{truth}})^T C^{-1}(\hat a-a_{\mathrm{truth}})
\]

are **diagnostics**. Do not quote them as physics results until the independent
final refit exists.

Per-station Kalman measurement residuals are not yet in the PDF: the persisted
EDM has \(\chi^2/\mathrm{ndf}\) and hit associations, not fitted residuals on
each B0 surface. Use \(\chi^2/\mathrm{ndf}\) until those states are exported.

## Historical numbers (not the benchmark)

`B0TrackerStubSeederConfig.h` comments quote a 3000-proton gun study:

| Quantity | Quoted value |
|---|---:|
| Seed direction | 13–24 μrad |
| Seed \(q/p\) | 2.4–4.3 % |
| Seed pull widths | 0.85–1.0 |
| Seeded four-station events | 765 |
| Stub-seeded reconstructed | 698 |
| Truth-seeded reconstructed | 723 |

Those yields are **conditional on an already-seeded four-station sample**. They
are not \(\epsilon_{\mathrm{total}}\) as defined above. Treat them as historical
until this command reproduces them on the pinned sample; then record the new
numbers in the run’s `report.md`.

## Failure counters

Set `B0_TRACKING_COUNTERS_FILE` before `eicrecon`. At exit the process writes
JSON with seeder cuts, including `overlapRejectedSubset` vs
`overlapRejectedUnrelated`, plus stub and truth CKF/ambiguity counts. Summing
`overlapRejected*` is the planned `overlapRejected` total.

These counters must not change track parameters or event acceptance.
