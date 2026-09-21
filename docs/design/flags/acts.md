## ACTS flags


### Logging

The `acts:LogLevel` sets log level for most operations. Currently, some may require source code changes to update the verbosity.

### Material map

Material map in JSon format can be loaded with **acts:MaterialMap** flag:

```yaml
acts:MaterialMap=/path/to/file/material.cbp
```

The default value for MaterialMap `calibrations/materials-map.cbor`.
When EICRecon runs, DD4Hep downloads calibrations to the current running directory
including material map to `calibrations/materials-map.cbor`.

### Geometry envelopes

The DD4hep to ACTS conversion pads the bounds of every subdetector tracking
volume beyond its outermost layers. The pads are exposed as

```yaml
acts:LayerEnvelopeR=1   # mm, radial pad
acts:LayerEnvelopeZ=5   # mm, longitudinal pad
```

They are per-volume pads, not per-layer envelopes (those come from the DD4hep
`envelope_r/z_min/max` layer parameters). The z pad must cover `r*tan(tilt)`
for layers tilted off the beam axis, such as the B0 tracker at 25 mrad, or the
tilted sensors sit behind the volume entry point and the navigator skips them
without recording a hole. It must also stay below half the z gap between
neighbouring tracking volumes, otherwise geometry conversion fails with
"Misconfiguration in volume building".
