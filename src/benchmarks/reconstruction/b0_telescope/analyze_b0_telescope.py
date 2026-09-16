#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Surface-consistent B0 evaluation; truth matching is not an analysis selection.

The pure evaluation functions accept dictionaries for dependency-light tests.
PODIO is imported only by the ROOT adapter. No simulation performance is baked in.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

DEFAULT_CHAINS = ("B0TelescopeCKF", "B0TelescopeDirect", "B0TelescopeTruthCKF", "B0TelescopeTruthDirect")
ParticleKey = tuple[int, int, int]


def strict_json(text: str) -> Any:
    def invalid(value: str) -> None:
        raise ValueError(f"Non-finite JSON number: {value}")
    return json.loads(text, parse_constant=invalid)


def read_jsonl(path: Path) -> Iterable[dict]:
    with path.open(encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                row = strict_json(line)
                if not isinstance(row, dict):
                    raise ValueError("Expected a JSON object")
                yield row
            except (ValueError, TypeError) as error:
                raise ValueError(f"{path}:{number}: {error}") from error


def particle_key(event: int, identity: Iterable[int]) -> ParticleKey:
    parts = tuple(identity)
    if len(parts) != 2 or any(not isinstance(i, int) for i in parts):
        raise ValueError("Particle identity must contain collection ID and index")
    return int(event), parts[0], parts[1]


@dataclass
class TruthReference:
    particles: dict[ParticleKey, dict] = field(default_factory=dict)
    states: dict[tuple[ParticleKey, str], dict] = field(default_factory=dict)
    seeds: dict[tuple[int, int], ParticleKey] = field(default_factory=dict)
    events: set[int] = field(default_factory=set)
    configuration: dict = field(default_factory=dict)
    counters: Counter = field(default_factory=Counter)
    event_counters: dict[int, dict] = field(default_factory=dict)

    @classmethod
    def read(cls, path: Path) -> "TruthReference":
        result = cls()
        for row in read_jsonl(path):
            kind = row.get("type")
            if kind == "configuration":
                if result.configuration or row.get("schema_version") != 1:
                    raise ValueError("Missing/duplicate/unsupported truth configuration")
                result.configuration = row
                continue
            event = int(row["event"])
            if kind == "event":
                if event in result.events:
                    raise ValueError(f"Duplicate event ID {event}; evaluate separate runs separately")
                result.events.add(event)
                result.event_counters[event] = {k: v for k, v in row.items() if k not in ("type", "event")}
                result.counters.update(result.event_counters[event])
            elif kind == "particle":
                key = particle_key(event, row["particle_id"])
                if key in result.particles:
                    raise ValueError(f"Duplicate truth particle {key}")
                result.particles[key] = row
            elif kind == "truth_state":
                key = (particle_key(event, row["particle_id"]), str(row["surface"]))
                if key in result.states:
                    raise ValueError(f"Duplicate truth surface state {key}")
                state = np.asarray(row["state"], dtype=float)
                if state.shape != (6,) or not np.isfinite(state).all():
                    raise ValueError("Invalid local truth state")
                result.states[key] = row
            elif kind == "truth_seed":
                key = (event, int(row["seed"]))
                if key in result.seeds:
                    raise ValueError(f"Duplicate truth seed {key}")
                result.seeds[key] = particle_key(event, row["particle_id"])
            else:
                raise ValueError(f"Unknown truth record type {kind!r}")
        if not result.configuration:
            raise ValueError("Truth reference has no schema/projection configuration")
        return result


@dataclass(frozen=True)
class Selection:
    pdgs: frozenset[int] = frozenset()
    min_stations: int = 3
    min_production_momentum: float = 0.0
    match_threshold: float = 0.5

    def __post_init__(self) -> None:
        if self.min_stations < 3 or not math.isfinite(self.min_production_momentum) or self.min_production_momentum < 0:
            raise ValueError("Invalid reconstructibility selection")
        if not (0.5 <= self.match_threshold < 1):
            raise ValueError("Strict-majority threshold must be in [0.5,1)")

    def accepts(self, row: dict) -> bool:
        return (bool(row["supported"]) and (not self.pdgs or row["pdg"] in self.pdgs)
                and row["n_stations"] >= self.min_stations
                and row["p_production_gev"] >= self.min_production_momentum)


def match_particle(matches: list[dict], threshold: float) -> tuple[tuple[int, int] | None, str]:
    """Return a strict majority; ties and contradictory associations are not resolved by order."""
    weights: dict[tuple[int, int], float] = defaultdict(float)
    for item in matches:
        weight = float(item["weight"])
        if not math.isfinite(weight) or not (0 <= weight <= 1):
            raise ValueError("Association weight is not a finite probability")
        weights[tuple(item["particle_id"])] += weight
    winners = [identity for identity, weight in weights.items() if weight > threshold]
    if len(winners) == 1:
        return winners[0], "matched"
    return None, "ambiguous" if weights else "unassociated"


def rate(successes: int, total: int) -> dict:
    """Binomial rate with a 95% Wilson interval, undefined rather than zero for no denominator."""
    if not 0 <= successes <= total:
        raise ValueError("Invalid binomial counts")
    if total == 0:
        return {"numerator": successes, "denominator": total, "value": None, "low95": None, "high95": None}
    z = 1.959963984540054
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return {"numerator": successes, "denominator": total, "value": p,
            "low95": max(0.0, centre - half), "high95": min(1.0, centre + half)}


def residual_diagnostics(estimate: dict, truth: dict) -> dict:
    if str(estimate["surface"]) != str(truth["surface"]):
        raise ValueError("Cannot compare bound parameters on different surfaces")
    a = np.asarray(estimate["state"], dtype=float)
    b = np.asarray(truth["state"], dtype=float)
    cov = np.asarray(estimate["covariance"], dtype=float)
    if a.shape not in ((5,), (6,)) or b.shape not in ((5,), (6,)) or cov.shape not in ((5, 5), (6, 6)):
        raise ValueError("Expected five spatial bound parameters and their covariance")
    a, b, cov = a[:5], b[:5], cov[:5, :5]
    if not np.isfinite(a).all() or not np.isfinite(b).all() or not np.isfinite(cov).all():
        raise ValueError("Non-finite parameters/covariance")
    if not np.allclose(cov, cov.T, rtol=1e-6, atol=1e-12):
        raise ValueError("Asymmetric covariance")
    cov = (cov + cov.T) * 0.5
    try:
        lower = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as error:
        raise ValueError("Non-positive-definite covariance; inverse alone is insufficient") from error
    delta = a - b
    delta[2] = (delta[2] + math.pi) % (2 * math.pi) - math.pi
    whitened = np.linalg.solve(lower, delta)
    pulls = delta / np.sqrt(np.diag(cov))
    q, truth_q = float(a[4]), float(b[4])
    return {"residual": delta.tolist(), "pull": pulls.tolist(), "d2": float(whitened @ whitened),
            "wrong_charge": bool(q * truth_q < 0), "unresolved_charge": q == 0,
            "relative_momentum_residual": None if q == 0 or truth_q == 0 else abs(truth_q / q) - 1}


def quantiles(values: list[float]) -> dict:
    finite = np.asarray([x for x in values if x is not None and math.isfinite(x)], dtype=float)
    if not finite.size:
        return {"n": 0}
    return {"n": int(finite.size), "mean": float(finite.mean()), "std": float(finite.std()),
            "q025": float(np.quantile(finite, .025)), "q16": float(np.quantile(finite, .16)),
            "median": float(np.median(finite)), "q84": float(np.quantile(finite, .84)),
            "q975": float(np.quantile(finite, .975)), "p95": float(np.quantile(finite, .95)),
            "p99": float(np.quantile(finite, .99))}


def evaluate(events: Iterable[dict], truth: TruthReference, selection: Selection,
             chains: tuple[str, ...] = DEFAULT_CHAINS) -> tuple[dict, list[dict]]:
    seen: set[int] = set()
    target: set[ParticleKey] = set()
    seeded = {name: set() for name in chains}
    reconstructed = {name: set() for name in chains}
    counters = {name: Counter() for name in chains}
    output_rows: list[dict] = []
    occupancy: dict[int, int] = {}
    by_event: dict[int, list[tuple[ParticleKey, dict]]] = defaultdict(list)
    for key, row in truth.particles.items():
        by_event[key[0]].append((key, row))
    for event in events:
        number = int(event["event"])
        if number in seen:
            raise ValueError(f"Duplicate event {number}; do not silently combine run-local IDs")
        if number not in truth.events:
            raise ValueError(f"ROOT event {number} is absent from the truth sidecar")
        seen.add(number)
        occupancy[number] = int(event["n_hits"])
        current = {key for key, row in by_event[number] if selection.accepts(row)}
        target.update(current)
        for name in chains:
            if name not in event["chains"]:
                raise ValueError(f"Required reconstruction chain {name} missing")
            chain = event["chains"][name]
            count = counters[name]
            count["events"] += 1
            count["seeds"] += len(chain["seeds"])
            for seed in chain["seeds"]:
                identity, _ = match_particle(seed["matches"], selection.match_threshold)
                if identity is not None:
                    key = particle_key(number, identity)
                    if key in current:
                        seeded[name].add(key)
            matched_tracks: dict[ParticleKey, list[dict]] = defaultdict(list)
            for track in chain["tracks"]:
                count["tracks"] += 1
                identity, reason = match_particle(track["matches"], selection.match_threshold)
                if identity is None:
                    count["fake_" + reason] += 1
                    continue
                key = particle_key(number, identity)
                if key not in current:
                    count["matched_outside_target"] += 1
                    continue
                matched_tracks[key].append(track)
            for key, alternatives in matched_tracks.items():
                reconstructed[name].add(key)
                count["duplicate_tracks"] += len(alternatives) - 1
                best = min(alternatives, key=lambda t: (-t["n_measurements"],
                    t["chi2"] / max(1, t["ndf"]), tuple(t["id"])))
                estimate = best.get("parameters")
                if estimate is None:
                    count["missing_parameters"] += 1
                    continue
                reference = truth.states.get((key, str(estimate["surface"])))
                if reference is None:
                    count["missing_same_surface_truth"] += 1
                    continue
                try:
                    diagnostic = residual_diagnostics(estimate, reference)
                except ValueError:
                    count["invalid_covariance_or_parameters"] += 1
                    continue
                particle = truth.particles[key]
                output_rows.append({"chain": name, "event": number, "particle_id": list(key[1:]),
                    "track_id": best["id"], "measurement_ids": best.get("measurement_ids", []),
                    "surface": str(estimate["surface"]), "state": estimate["state"],
                    "covariance": estimate["covariance"], "pdg": particle["pdg"],
                    "n_stations": particle["n_stations"], "p_production_gev": particle["p_production_gev"],
                    "vertex_z_mm": particle["vertex_mm"][2], "n_hits": occupancy[number], **diagnostic})
    if not seen:
        raise ValueError("No events were evaluated")
    truth_counts: Counter = Counter()
    for number in seen:
        truth_counts.update(truth.event_counters.get(number, {}))
    report = {"schema_version": 1, "events": len(seen), "target_particles": len(target),
        "selection": {"signed_pdgs": sorted(selection.pdgs), "empty_pdg_list": "all supported charged species",
            "min_stations": selection.min_stations, "min_production_momentum_gev": selection.min_production_momentum,
            "match_rule": f"strict purity > {selection.match_threshold}"},
        "truth_binding": truth.configuration, "truth_stage_counters": dict(truth_counts),
        "reference": "actual sensitive plane, five spatial bound parameters", "chains": {}}
    bins = {"p_production_gev": [0, 3, 8, 15, 25, 41, 80, 150, 300, math.inf],
            "vertex_z_mm": [-math.inf, -200, 200, 1000, 3000, 5000, 6000, 7000, math.inf],
            "n_stations": [2.5, 3.5, 4.5, math.inf], "n_hits": [0, 8, 16, 32, 64, 128, 512, math.inf]}
    for name in chains:
        count = counters[name]
        reco = reconstructed[name]
        seeds = seeded[name]
        rows = [row for row in output_rows if row["chain"] == name]
        fakes = count["fake_unassociated"] + count["fake_ambiguous"]
        summary = {"counts": dict(count), "seed_efficiency": rate(len(seeds), len(target)),
            "track_given_seed": rate(len(reco & seeds), len(seeds)), "total_efficiency": rate(len(reco), len(target)),
            "recovered_without_matching_seed": len(reco - seeds), "fake_fraction": rate(fakes, count["tracks"]),
            "duplicate_fraction_of_tracks": rate(count["duplicate_tracks"], count["tracks"]),
            "matched_non_target_fraction": rate(count["matched_outside_target"], count["tracks"]),
            "relative_momentum_residual": quantiles([r["relative_momentum_residual"] for r in rows]),
            "d2": quantiles([r["d2"] for r in rows]),
            "d2_coverage95": rate(sum(r["d2"] <= 11.070497693516351 for r in rows), len(rows)),
            "wrong_charge": rate(sum(r["wrong_charge"] for r in rows), len(rows)),
            "unresolved_charge": sum(r["unresolved_charge"] for r in rows),
            "pulls": {label: quantiles([r["pull"][i] for r in rows])
                      for i, label in enumerate(("loc0", "loc1", "phi", "theta", "qOverP"))}, "binned_efficiency": {}}
        for variable, edges in bins.items():
            summaries = []
            for low, high in zip(edges, edges[1:]):
                def value(key: ParticleKey) -> float:
                    row = truth.particles[key]
                    return (occupancy[key[0]] if variable == "n_hits" else row["vertex_mm"][2]
                            if variable == "vertex_z_mm" else row[variable])
                members = {key for key in target if low <= value(key) < high}
                summaries.append({"low": low if math.isfinite(low) else None,
                    "high": high if math.isfinite(high) else None, "seed": rate(len(members & seeds), len(members)),
                    "track": rate(len(members & reco), len(members))})
            summary["binned_efficiency"][variable] = summaries
        report["chains"][name] = summary
    return report, output_rows


def _object_id(obj: Any) -> tuple[int, int]:
    identity = obj.getObjectID()
    return int(identity.collectionID), int(identity.index)


def _required(frame: Any, name: str) -> Any:
    try:
        collection = frame.get(name)
    except Exception as error:
        raise ValueError(f"Missing persisted collection {name}") from error
    if collection is None:
        raise ValueError(f"Missing persisted collection {name}")
    return collection


def _bound(parameters: Any) -> dict:
    local = parameters.getLoc()
    cov = parameters.getCovariance()
    return {"surface": str(parameters.getSurface()),
        "state": [float(local.a), float(local.b), float(parameters.getPhi()), float(parameters.getTheta()),
                  float(parameters.getQOverP()), float(parameters.getTime())],
        "covariance": [[float(cov(i, j)) for j in range(6)] for i in range(6)]}


def podio_events(path: Path, truth: TruthReference, chains: tuple[str, ...]) -> Iterable[dict]:
    from podio.reading import get_reader
    reader = get_reader(str(path))
    for frame in reader.get("events"):
        headers = _required(frame, "EventHeader")
        if len(headers) != 1:
            raise ValueError("Exactly one EventHeader required")
        number = int(headers[0].getEventNumber())
        raw_votes: dict[tuple[int, int], dict[tuple[int, int], float]] = defaultdict(lambda: defaultdict(float))
        for a in _required(frame, "B0TrackerRawHitAssociations"):
            w = float(a.getWeight())
            if not math.isfinite(w) or w < 0:
                raise ValueError("Invalid raw-hit truth weight")
            particle = a.getSimHit().getParticle()
            if particle.isAvailable() and w > 0:
                raw_votes[_object_id(a.getRawHit())][_object_id(particle)] += w
        result = {"event": number, "n_hits": len(_required(frame, "B0TrackerMeasurements")), "chains": {}}
        for stem in chains:
            truth_chain = stem.startswith("B0TelescopeTruth")
            seed_name = "B0TelescopeTruthSeeds" if truth_chain else "B0TelescopeSeeds"
            seeds = []
            for index, seed in enumerate(_required(frame, seed_name)):
                if truth_chain:
                    identity = truth.seeds.get((number, index))
                    if identity is None:
                        raise ValueError(f"Persisted truth seed {number}:{index} missing from sidecar")
                    matches = [{"particle_id": identity[1:], "weight": 1.0}]
                else:
                    hits = {_object_id(h): h for h in seed.getHits()}
                    votes: dict[tuple[int, int], float] = defaultdict(float)
                    for h in hits.values():
                        contributions = raw_votes.get(_object_id(h.getRawHit()), {})
                        total = sum(contributions.values())
                        if total > 0:
                            for identity, weight in contributions.items():
                                votes[identity] += weight / total
                    matches = [{"particle_id": identity, "weight": weight / len(hits)}
                               for identity, weight in votes.items()] if hits else []
                seeds.append({"id": _object_id(seed), "matches": matches})
            matches_by_track: dict[tuple[int, int], list[dict]] = defaultdict(list)
            for a in _required(frame, stem + "Associations"):
                matches_by_track[_object_id(a.getRec())].append({"particle_id": _object_id(a.getSim()), "weight": float(a.getWeight())})
            tracks = []
            for t in _required(frame, stem + "Tracks"):
                parameters = list(t.getTrajectory().getTrackParameters())
                if len(parameters) > 1:
                    raise ValueError("Local exporter must persist one unambiguous reference state per track")
                identity = _object_id(t)
                measurements = list(t.getMeasurements())
                tracks.append({"id": identity, "n_measurements": len(measurements),
                    "measurement_ids": sorted(_object_id(m) for m in measurements),
                    "chi2": float(t.getChi2()), "ndf": int(t.getNdf()),
                    "parameters": _bound(parameters[0]) if parameters else None, "matches": matches_by_track[identity]})
            result["chains"][stem] = {"seeds": seeds, "tracks": tracks}
        yield result


def diagnostic_summary(paths: list[Path]) -> dict:
    result = {}
    for path in paths:
        count: Counter = Counter()
        elapsed = []
        configurations = []
        events = set()
        for row in read_jsonl(path):
            if row["type"] == "configuration":
                configurations.append(row)
            elif row["type"] in ("event", "seed_event"):
                event = int(row["event"])
                if event in events:
                    raise ValueError(f"Duplicate diagnostic event in {path}")
                events.add(event)
                for key, value in row.items():
                    if key not in ("type", "event", "elapsed_ms") and isinstance(value, (int, float)):
                        count[key] += value
                if "elapsed_ms" in row:
                    elapsed.append(float(row["elapsed_ms"]))
        result[str(path)] = {"configuration": configurations, "events": len(events), "totals": dict(count),
                             "elapsed_ms": quantiles(elapsed)}
    return result


def write_report(report: dict, rows: list[dict], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (out / "matched_tracks.jsonl").open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, allow_nan=False) + "\n")
    def percent(value: dict) -> str:
        return "undefined" if value["value"] is None else f'{100 * value["value"]:.2f}% ({value["numerator"]}/{value["denominator"]})'
    lines = ["# B0 local tracking evaluation", "", f'Events: {report["events"]}; target particles: {report["target_particles"]}.',
        "", "Reference: actual sensor plane; five spatial bound parameters. Matching uses strict majority purity.",
        "Real tracks outside the selected species/phase space are **not fakes**. Missing/invalid reference states are counted, not assigned zero residual.",
        "", "| Chain | Seed efficiency | Track given seed | Total efficiency | Fake fraction |", "|---|---:|---:|---:|---:|"]
    for name, values in report["chains"].items():
        lines.append(f'| {name} | {percent(values["seed_efficiency"])} | {percent(values["track_given_seed"])} | {percent(values["total_efficiency"])} | {percent(values["fake_fraction"])} |')
    lines.extend(["", "The JSON report contains Wilson intervals, duplicate/non-target rates, charge errors, residual/pull quantiles, D² diagnostics, binned efficiencies and rejection counters.",
        "Do not infer calibrated coverage from a refit flag. Validate prior independence, truth binding, material, selection effects and model assumptions first.",
        "The seed/conditional efficiencies need not multiply to the total if CKF reconstructs a particle from a seed matched to another particle; such cases are counted explicitly.", ""])
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reco", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--chain", action="append", default=[])
    parser.add_argument("--pdg", type=int, action="append", default=[], help="signed PDG; repeat, empty means all supported")
    parser.add_argument("--min-stations", type=int, default=3)
    parser.add_argument("--min-production-momentum", type=float, default=0)
    parser.add_argument("--match-threshold", type=float, default=.5)
    parser.add_argument("--diagnostics", action="append", type=Path, default=[])
    parser.add_argument("--provenance", type=Path)
    args = parser.parse_args()
    chains = tuple(args.chain or DEFAULT_CHAINS)
    if len(set(chains)) != len(chains):
        parser.error("Duplicate chain names")
    truth = TruthReference.read(args.truth)
    selection = Selection(frozenset(args.pdg), args.min_stations, args.min_production_momentum, args.match_threshold)
    report, rows = evaluate(podio_events(args.reco, truth, chains), truth, selection, chains)
    report["diagnostics"] = diagnostic_summary(args.diagnostics)
    if args.provenance:
        from provenance import verify_outputs
        report["provenance"] = strict_json(args.provenance.read_text(encoding="utf-8"))
        verify_outputs(report["provenance"], args.reco, args.truth)
    write_report(report, rows, args.out)
    print(args.out / "report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
