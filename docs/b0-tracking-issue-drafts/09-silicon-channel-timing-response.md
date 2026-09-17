# SiliconTrackerDigi: define channel-level timing semantics before using B0 timing for pattern recognition

## Motivation

The shared `SiliconTrackerDigi` currently smears each `SimTrackerHit` time independently and assigns the raw channel the minimum smeared time among contributions in that cell.

For a multi-contribution channel this is not equivalent to a single detector/electronics timing measurement. Even when two contributions have the same true time and independent Gaussian jitter with width `sigma_t`, the expected minimum is biased early:

```text
E[min(epsilon_1, epsilon_2)] = -sigma_t / sqrt(pi)
```

For the B0 configuration (`timeResolution = 30 ps`) two simultaneous contributions therefore produce an illustrative bias of about `-16.9 ps` purely from the aggregation rule. The bias changes with contribution multiplicity.

This matters before the B0 candidate builder enables time-of-flight compatibility by default. Timing should reject incompatible trajectories, not encode an occupancy-dependent implementation artifact.

Related work:

- PR #8 fixes the separate channel-threshold / truth-association inconsistency.
- The ePIC geometry follow-up proposes a realistic AC-LGAD response model with charge sharing, gain/noise, threshold and timing.

## Goal

Give silicon channels an explicit, documented timing-response contract. The generic digitizer should not accidentally derive a channel timestamp by taking an order statistic of independently smeared deposits unless that behavior represents a deliberate discriminator model.

## Proposed staged approach

### 1. Define generic channel timing semantics

For the shared silicon digitizer, choose and document a channel-level rule. A minimal detector-agnostic model could be:

1. collect all contributions belonging to the channel/readout window;
2. determine a deterministic channel time from unsmeared contributions (for example first-arrival time when that is the intended discriminator approximation);
3. smear the **channel measurement once** by the configured resolution.

Do not silently mix independent per-contribution jitter with a single reported channel uncertainty.

### 2. Keep detector-specific effects out of the generic fallback

AC-LGAD-specific response can add configurable effects such as:

- amplitude-dependent time walk;
- neighboring-electrode pulse combination;
- gain/noise dependence;
- threshold crossing / discriminator model;
- common clock or event-time uncertainty.

The shared fallback should remain simple and explicit rather than pretending to be the final B0 electronics model.

### 3. Propagate the meaning of the timestamp

Document whether `RawTrackerHit.timeStamp` represents:

- first threshold crossing;
- pulse time estimator;
- earliest physical deposit;
- another channel observable.

`TrackerHitReconstruction` should attach a covariance/resolution consistent with that definition.

## Tests

- [ ] One contribution: configured resolution is preserved.
- [ ] Two simultaneous contributions do not acquire an order-statistic bias unless deliberately configured.
- [ ] Contribution ordering does not change the result for a fixed RNG seed.
- [ ] Separated contributions follow the documented discriminator/aggregation rule.
- [ ] Channel multiplicity does not silently change the meaning of the reported uncertainty.
- [ ] B0 timing-candidate tests validate efficiency/background rejection with the selected response model.

## Acceptance criteria

- [ ] Channel timing semantics are explicit in code/configuration documentation.
- [ ] The generic implementation applies timing smearing at the correct aggregation level.
- [ ] B0's 30 ps configuration has a documented physical interpretation.
- [ ] Timing-based B0 candidate rejection remains default-off until detector-backed validation demonstrates calibrated behavior.

## Non-goals

- Requiring a full waveform simulation for the generic silicon digitizer.
- Encoding one AC-LGAD design as universal silicon-detector behavior.
