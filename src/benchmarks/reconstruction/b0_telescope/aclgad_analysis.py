#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""AC-LGAD cluster-aware adapter for the existing B0 surface-consistent benchmark."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path


def measurement_votes(hits: list[dict], weights: list[float], raw_votes: dict) -> dict:
    """Mirror B0MeasurementTruth.h: noise stays in the denominator."""
    if not hits:
        return {}
    if weights and len(weights) != len(hits):
        raise ValueError("Measurement hits/weights length mismatch")
    if len({tuple(h['id']) for h in hits}) != len(hits):
        raise ValueError("Duplicate measurement channel")
    weights = weights or [1.0] * len(hits)
    if any(not math.isfinite(w) or w < 0 for w in weights):
        raise ValueError("Invalid measurement weight")
    total = sum(weights)
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Invalid measurement weight sum")
    result = defaultdict(float)
    for hit, weight in zip(hits, weights):
        votes = raw_votes.get(tuple(hit['raw']), {})
        if any(not math.isfinite(w) or w < 0 for w in votes.values()):
            raise ValueError("Invalid raw truth weight")
        norm = max(1.0, sum(votes.values()))
        if not math.isfinite(norm):
            raise ValueError("Invalid raw truth sum")
        for particle, fraction in votes.items():
            result[particle] += weight / total * fraction / norm
    return dict(result)


def seed_votes(hit_ids: list[tuple], measurements: dict, hit_to_measurement: dict, raw_votes: dict) -> dict:
    """Average physical measurements, not the number of pads in each cluster."""
    ids = set(hit_ids)
    if len(ids) != len(hit_ids):
        raise ValueError("Duplicate seed hit")
    if any(identity not in hit_to_measurement for identity in ids):
        raise ValueError("Seed hit has no measurement")
    selected = {hit_to_measurement[identity] for identity in ids}
    votes = defaultdict(float)
    for identity in selected:
        m = measurements[identity]
        if not {tuple(h['id']) for h in m['hits']} <= ids:
            raise ValueError("Seed contains only part of a cluster")
        for particle, weight in measurement_votes(m['hits'], m['weights'], raw_votes).items():
            votes[particle] += weight / len(selected)
    return dict(votes)


def events(path: Path, truth, chains: tuple[str, ...]):
    from podio.reading import get_reader
    from analyze_b0_telescope import _required, _object_id, _bound
    reader = get_reader(str(path))
    for frame in reader.get('events'):
        headers = _required(frame, 'EventHeader')
        if len(headers) != 1:
            raise ValueError('Exactly one EventHeader required')
        number = int(headers[0].getEventNumber())
        raw = defaultdict(lambda: defaultdict(float))
        for association in _required(frame, 'B0ACLGADRawHitAssociations'):
            weight = float(association.getWeight())
            if not math.isfinite(weight) or weight < 0:
                raise ValueError('Invalid raw association weight')
            particle = association.getSimHit().getParticle()
            if particle.isAvailable() and weight > 0:
                raw[_object_id(association.getRawHit())][_object_id(particle)] += weight
        measurements, reverse = {}, {}
        for m in _required(frame, 'B0ACLGADMeasurements'):
            identity = _object_id(m)
            hits = [{'id': _object_id(h), 'raw': _object_id(h.getRawHit())} for h in m.getHits()]
            for h in hits:
                if h['id'] in reverse:
                    raise ValueError('A channel belongs to multiple AC-LGAD measurements')
                reverse[h['id']] = identity
            cov = m.getCovariance()
            measurements[identity] = {'hits': hits, 'weights': [float(w) for w in m.getWeights()],
                'surface': str(m.getSurface()), 'state': [float(m.getLoc().a), float(m.getLoc().b), float(m.getTime())],
                'covariance': [[float(cov.xx), float(cov.xy), float(cov.xz)],
                               [float(cov.xy), float(cov.yy), float(cov.yz)],
                               [float(cov.xz), float(cov.yz), float(cov.zz)]]}
        result = {'event': number, 'n_hits': len(measurements), 'chains': {},
                  'response_measurements': list(measurements.values()), 'raw_votes': dict(raw),
                  'n_channels': len(_required(frame, 'B0ACLGADRawHits'))}
        for stem in chains:
            truth_chain = stem.startswith('B0TelescopeTruth')
            seed_name = 'B0TelescopeTruthSeeds' if truth_chain else 'B0TelescopeSeeds'
            seeds = []
            for index, seed in enumerate(_required(frame, seed_name)):
                if truth_chain:
                    identity = truth.seeds.get((number, index))
                    if identity is None:
                        raise ValueError('Truth seed is absent from sidecar')
                    matches = [{'particle_id': identity[1:], 'weight': 1.0}]
                else:
                    votes = seed_votes([_object_id(h) for h in seed.getHits()], measurements, reverse, raw)
                    matches = [{'particle_id': identity, 'weight': weight} for identity, weight in votes.items()]
                seeds.append({'id': _object_id(seed), 'matches': matches})
            matches_by_track = defaultdict(list)
            for association in _required(frame, stem + 'Associations'):
                matches_by_track[_object_id(association.getRec())].append(
                    {'particle_id': _object_id(association.getSim()), 'weight': float(association.getWeight())})
            tracks = []
            for track in _required(frame, stem + 'Tracks'):
                parameters = list(track.getTrajectory().getTrackParameters())
                if len(parameters) != 1:
                    raise ValueError('Expected one unambiguous local reference state')
                ms = list(track.getMeasurements())
                tracks.append({'id': _object_id(track), 'n_measurements': len(ms),
                    'measurement_ids': sorted(_object_id(m) for m in ms),
                    'chi2': float(track.getChi2()), 'ndf': int(track.getNdf()),
                    'parameters': _bound(parameters[0]), 'matches': matches_by_track[_object_id(track)]})
            result['chains'][stem] = {'seeds': seeds, 'tracks': tracks}
        yield result


def response_rows(event: dict, truth):
    import numpy as np
    from analyze_b0_telescope import match_particle, particle_key
    for measurement in event['response_measurements']:
        votes = measurement_votes(measurement['hits'], measurement['weights'], event['raw_votes'])
        identity, status = match_particle([{'particle_id': key, 'weight': w} for key, w in votes.items()], .5)
        row = {'event': event['event'], 'surface': measurement['surface'],
               'channels': len(measurement['hits']), 'match_status': status}
        if identity is not None:
            key = particle_key(event['event'], identity)
            reference = truth.states.get((key, measurement['surface']))
            row['particle_id'] = list(identity)
            if reference is None:
                row['reference_status'] = 'missing'
            else:
                ref = np.asarray(reference['state'], dtype=float)[[0, 1, 5]]
                covariance = np.asarray(measurement['covariance'], dtype=float)
                estimate = np.asarray(measurement['state'], dtype=float)
                if not np.isfinite(covariance).all() or not np.isfinite(estimate).all():
                    raise ValueError('Non-finite AC-LGAD measurement')
                if not np.allclose(covariance, covariance.T, atol=1e-12):
                    raise ValueError('Asymmetric AC-LGAD covariance')
                lower = np.linalg.cholesky(covariance)
                delta = estimate - ref
                whitened = np.linalg.solve(lower, delta)
                row.update(reference_status='matched', residual_mm_mm_ns=delta.tolist(),
                    pull=(delta / np.sqrt(covariance.diagonal())).tolist(), d2=float(whitened @ whitened))
        yield row


def main() -> None:
    from analyze_b0_telescope import (DEFAULT_CHAINS, TruthReference, Selection, evaluate,
                                      write_report, diagnostic_summary)
    from provenance import verify_outputs
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reco', 'truth', 'out', 'provenance'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--chain', action='append', choices=DEFAULT_CHAINS, default=[])
    parser.add_argument('--pdg', type=int, action='append', default=[])
    parser.add_argument('--diagnostics', type=Path, action='append', default=[])
    args = parser.parse_args()
    provenance = json.loads(args.provenance.read_text())
    verify_outputs(provenance, args.reco, args.truth)
    truth = TruthReference.read(args.truth)
    chains = tuple(args.chain or DEFAULT_CHAINS)
    args.out.mkdir(parents=True, exist_ok=False)
    with (args.out / 'cluster_residuals.jsonl').open('w') as stream:
        def observed():
            for event in events(args.reco, truth, chains):
                for row in response_rows(event, truth):
                    stream.write(json.dumps(row, allow_nan=False) + '\n')
                yield event
        report, rows = evaluate(observed(), truth, Selection(frozenset(args.pdg)), chains)
    report['response_model'] = 'physical-pitch phenomenological AC-LGAD; calibration is external'
    report['provenance'] = provenance
    report['diagnostics'] = diagnostic_summary(args.diagnostics)
    write_report(report, rows, args.out)
    print(args.out / 'report.md')


if __name__ == '__main__':
    main()
