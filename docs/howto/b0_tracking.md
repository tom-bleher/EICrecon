# B0 tracking with ACTS 47.7

Use the same `epic_ip6_extended` geometry and 5×41 GeV optics in simulation and
reconstruction. This reduced configuration contains B0 tracking, not the central
tracker. Its material map is specific to the ACTS geometry identifiers produced
by that configuration.

## Select the reconstruction installation explicitly

An ACTS overlay does not select an EICrecon checkout automatically. Inside the
matching eic-shell, activate the geometry, ACTS, and an EICrecon installation
built against that ACTS version:

```bash
source /path/to/epic/install/bin/thisepic.sh epic_ip6_extended
source /path/to/acts-47.7/install/bin/this_acts.sh
source /path/to/EICrecon/install/bin/eicrecon-this.sh
export JANA_PLUGIN_PATH="$EICrecon_ROOT/lib/EICrecon/plugins:/opt/local/lib/JANA/plugins"
command -v eicrecon
printf '%s\n' "$EICrecon_ROOT" "$DETECTOR_PATH/$DETECTOR_CONFIG.xml" "$JANA_PLUGIN_PATH"
```

Put these selections in the overlay's startup script as well. Keep its
`CMAKE_PREFIX_PATH` consistent with the selected reconstruction and ACTS prefixes.
The explicit plugin path excludes plugins from other EICrecon/ACTS builds. An
executable's configure-time version string can lag incremental source builds;
record the source commit and actual loaded library paths with validation runs.

## Simulation and material

Build the detector with `-DEPIC_BUILD_DDG4_PLUGINS=ON` when using its
`Geant4TrackerWeightedPlacementAction`. Select this action explicitly in the B0
simulation steering; the presence of its library does not change npsim's default
action. The action keeps accumulated deposits distinct between reused sensor
placements. Use `.edm4hep.root` simulation output.

The geometry's `material-map` constant selects the reduced-geometry map. Run
where calibration and field resources resolve, or provide an absolute
`-Pacts:MaterialMap=/path/to/map.cbor`. B0 startup checks require usable material
on the approach surfaces marked for mapping by the geometry; a wrong or
incomplete map must fail the job. Vacuum bins are valid and do not themselves
indicate a missing map. Material recording/conversion without a material
decorator is a separate workflow.

## B0 particle hypothesis

Both B0 CKF chains default to the proton hypothesis, PDG 2212, for forward-proton
studies. Central tracking retains its pion default. These are fit mass hypotheses,
not particle identification; the charge sign comes from the seed. Override the
B0 defaults for another known species, for example pions:

```text
-Ptracking:B0TrackerCKFTrajectories:ParticleHypothesisPdg=211
-Ptracking:B0TrackerCKFTruthSeededTrajectories:ParticleHypothesisPdg=211
```

Before production, run a small matched simulation/reconstruction sample. Verify
measurement surfaces and local coordinates, material/navigation coverage over B0
acceptance, physical-station counts, persisted relations, and residual/pull
closure. A successful startup or positive covariance alone does not validate
tracking resolution or uncertainty coverage.
