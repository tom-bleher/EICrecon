#!/usr/bin/env python3
"""Compare two check_b0_multipoint.py reports on identical truth denominators."""
import argparse
import json
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('baseline')
p.add_argument('candidate')
p.add_argument('--bootstrap-samples', type=int, default=2000)
p.add_argument('--bootstrap-seed', type=int, default=271828)
a = p.parse_args()
if a.bootstrap_samples < 1:
    p.error('--bootstrap-samples must be positive')
left, right = [json.load(open(path)) for path in [a.baseline, a.candidate]]
assert left.get('provenance') == right.get('provenance'), 'Analysis configurations differ'
l = {tuple(r['key']): r for r in left['per_particle']}
r = {tuple(x['key']): x for x in right['per_particle']}
assert l.keys() == r.keys(), 'Truth acceptance denominators differ; compare identical input events/configuration'
common = sorted(k for k in l if l[k]['success'] and r[k]['success'])
width = lambda values: float(0.5*np.diff(np.percentile(values,[16,84]))[0]) if len(values) else None
x = np.array([l[k]['qop_residual_percent'] for k in common])
y = np.array([r[k]['qop_residual_percent'] for k in common])
rng = np.random.default_rng(a.bootstrap_seed)
bootstrap = []
for _ in range(a.bootstrap_samples):
    if not common:
        break
    indices = rng.integers(0, len(common), len(common))
    bootstrap.append(width(y[indices]) - width(x[indices]))
result = {'paired_width_difference_percent': width(y)-width(x) if common else None,
          'paired_width_difference_bootstrap95_percent': np.percentile(bootstrap,[2.5,97.5]).tolist() if bootstrap else None,
          'bootstrap': {'samples': a.bootstrap_samples, 'seed': a.bootstrap_seed, 'unit': 'matched primary proton', 'paired': True},
          'baseline': a.baseline, 'candidate': a.candidate, 'accepted_primary_protons': len(l),
          'baseline_success': sum(x['success'] for x in l.values()),
          'candidate_success': sum(x['success'] for x in r.values()),
          'gained': sum(not l[k]['success'] and r[k]['success'] for k in l),
          'lost': sum(l[k]['success'] and not r[k]['success'] for k in l),
          'paired_success': len(common),
          'paired_tails_above_20_percent': {'baseline': int(np.count_nonzero(np.abs(x)>20)), 'candidate': int(np.count_nonzero(np.abs(y)>20))},
          'paired_wrong_charge': {'baseline': int(np.count_nonzero(x < -100)), 'candidate': int(np.count_nonzero(y < -100))},
          'paired_half68_percent': {'baseline': width(x), 'candidate': width(y)},
          'first_hit_validation': {'sensor_bound_seeds': right['statistics'].get('sensor_bound_seeds', 0),
                                   'invalid_surface_matches': len(right['invalid_seed_surface_matches']),
                                   'tracks_missing_start_hit': len(right['tracks_missing_start_hit'])},
          'baseline_measurement_counts': left['track_measurement_counts'],
          'candidate_measurement_counts': right['track_measurement_counts']}
print(json.dumps(result,indent=2))
