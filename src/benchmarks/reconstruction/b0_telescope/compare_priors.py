#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare fresh-prior scans only for identical hits and reference surfaces."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from analyze_b0_telescope import quantiles, read_jsonl, strict_json


def compare(first: list[dict], second: list[dict]) -> dict:
    def index(rows: list[dict]) -> dict:
        result = {}
        for row in rows:
            key = (row['chain'], row['event'], tuple(row['particle_id']))
            if key in result:
                raise ValueError('Duplicate representative track identity')
            result[key] = row
        return result
    a, b = index(first), index(second)
    changed_hits = changed_surface = 0
    shifts, changes = [], []
    for key in sorted(a.keys() & b.keys()):
        x, y = a[key], b[key]
        if sorted(map(tuple, x['measurement_ids'])) != sorted(map(tuple, y['measurement_ids'])):
            changed_hits += 1; continue
        if str(x['surface']) != str(y['surface']):
            changed_surface += 1; continue
        if not x['measurement_ids']:
            raise ValueError('Prior comparison needs explicit measurement identities')
        dx = np.asarray(y['state'], dtype=float)[:5] - np.asarray(x['state'], dtype=float)[:5]
        dx[2] = (dx[2] + math.pi) % (2 * math.pi) - math.pi
        c, d = np.asarray(x['covariance'], float)[:5, :5], np.asarray(y['covariance'], float)[:5, :5]
        if not all(np.isfinite(v).all() for v in (dx, c, d)):
            raise ValueError('Non-finite prior-scan data')
        lower = np.linalg.cholesky(c)
        np.linalg.cholesky(d)
        shifts.append(float(np.linalg.norm(np.linalg.solve(lower, dx))))
        relative = np.linalg.solve(lower, d - c)
        relative = np.linalg.solve(lower, relative.T).T
        changes.append(float(np.linalg.norm(relative, ord=2)))
    return {'common_particles': len(a.keys() & b.keys()), 'only_first': len(a.keys() - b.keys()),
        'only_second': len(b.keys() - a.keys()), 'changed_hits': changed_hits,
        'changed_reference_surface': changed_surface, 'identical_hit_set_tracks': len(shifts),
        'mean_shift_in_reference_covariance': quantiles(shifts),
        'relative_covariance_spectral_change': quantiles(changes),
        'max_mean_shift': max(shifts, default=None), 'max_covariance_change': max(changes, default=None)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('other', type=Path, nargs='+')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-mean-shift', type=float, default=.01)
    parser.add_argument('--max-covariance-change', type=float, default=.01)
    args = parser.parse_args()
    if any(not math.isfinite(v) or v <= 0 for v in (args.max_mean_shift, args.max_covariance_change)):
        parser.error('Acceptance thresholds must be finite and positive')
    first_summary = strict_json((args.reference / 'summary.json').read_text())
    def fingerprint(summary: dict) -> tuple:
        manifest = summary.get('provenance', {})
        if manifest.get('status') != 'completed':
            raise ValueError('Prior scans require completed run provenance')
        return tuple(manifest['inputs'][name]['sha256'] for name in ('sim', 'compact', 'calibration', 'material'))
    fingerprint_first = fingerprint(first_summary)
    first = list(read_jsonl(args.reference / 'matched_tracks.jsonl'))
    output = {}
    failed = False
    for directory in args.other:
        summary = strict_json((directory / 'summary.json').read_text())
        if fingerprint(summary) != fingerprint_first or summary.get('selection') != first_summary.get('selection'):
            raise ValueError('Input geometry/material/sample or analysis selections differ')
        result = compare(first, list(read_jsonl(directory / 'matched_tracks.jsonl')))
        passed = (result['identical_hit_set_tracks'] > 0 and result['max_mean_shift'] <= args.max_mean_shift
                  and result['max_covariance_change'] <= args.max_covariance_change)
        result['passes_requested_tolerances_on_common_hits'] = passed
        failed |= not passed
        output[str(directory)] = result
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({'reference': str(args.reference), 'comparisons': output,
        'thresholds': {'mean': args.max_mean_shift, 'covariance': args.max_covariance_change}},
        indent=2, allow_nan=False) + '\n')
    return int(failed)

if __name__ == '__main__':
    raise SystemExit(main())
