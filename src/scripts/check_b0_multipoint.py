#!/usr/bin/env python3
"""Check seed surface/first-hit relations and truth-matched final CKF output.

All events are read. No truth-nearest candidate selection: final duplicate tracks
are ordered by descending measurement count, then chi2, then track index.
Momentum widths are 0.5*(p84-p16), NumPy linear interpolation, no rounding in JSON.
"""
import argparse
from collections import Counter
import json
import sys

import awkward as ak
import numpy as np
import uproot

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file')
parser.add_argument('--prefix', default='B0TrackerCKF')
parser.add_argument('--crossing-angle', type=float, default=-0.025, help='Outgoing-ion axis rotation in radians')
args = parser.parse_args()
if not np.isfinite(args.crossing_angle):
    parser.error('--crossing-angle must be finite')
f = uproot.open(args.file)
t = f['events']
a = lambda name: ak.to_list(t[name].array())
p = args.prefix
run, event = a('EventHeader.runNumber'), a('EventHeader.eventNumber')
assert all(len(r) == len(e) == 1 for r, e in zip(run, event)), 'Require exactly one EventHeader per event'
keys = [(r[0], e[0]) for r, e in zip(run, event)]
assert len(set(keys)) == len(keys), 'Duplicate event identity'
surface = a('B0TrackerSeedParameters.surface')
seed_cov = a('B0TrackerSeedParameters.covariance.covariance[21]')
seed_values = [a('B0TrackerSeedParameters.'+field) for field in ['loc.a','loc.b','phi','theta','qOverP','time']]
seed_param = a('_B0TrackerSeeds_params.index')
seed_b, seed_e = a('B0TrackerSeeds.hits_begin'), a('B0TrackerSeeds.hits_end')
seed_hits = a('_B0TrackerSeeds_hits.index')
meas_b, meas_e = a('B0TrackerMeasurements.hits_begin'), a('B0TrackerMeasurements.hits_end')
meas_hits = a('_B0TrackerMeasurements_hits.index')
meas_surface = a('B0TrackerMeasurements.surface')
traj_seed = a('_' + p + 'Trajectories_seed.index')
track_traj = a('_' + p + 'Tracks_trajectory.index')
track_b, track_e = a(p + 'Tracks.measurements_begin'), a(p + 'Tracks.measurements_end')
track_meas = a('_' + p + 'Tracks_measurements.index')
track_q = a(p + 'Tracks.charge')
track_chi2 = a(p + 'Tracks.chi2')
track_p = [a(p + 'Tracks.momentum.' + c) for c in 'xyz']
rec, sim, weight = a('_' + p + 'TrackAssociations_rec.index'), a('_' + p + 'TrackAssociations_sim.index'), a(p + 'TrackAssociations.weight')
mc_status, mc_pdg, mc_charge = a('MCParticles.generatorStatus'), a('MCParticles.PDG'), a('MCParticles.charge')
mc_p = [a('MCParticles.momentum.' + c) for c in 'xyz']
hit_part = a('_B0TrackerHits_particle.index')
hit_pos = [a('B0TrackerHits.position.' + c) for c in 'xyz']
hit_mom = [a('B0TrackerHits.momentum.' + c) for c in 'xyz']
# Validate collection identity too: array positions alone are not PODIO identity.
meta = f['podio_metadata']['events___CollectionTypeInfo']
md = ak.to_list(meta.arrays())[0]
ids = dict(zip(md['events___CollectionTypeInfo.name'], md['events___CollectionTypeInfo.collectionID']))
for relation, target in [('B0TrackerSeeds_params', 'B0TrackerSeedParameters'),
                         ('B0TrackerSeeds_hits', 'B0TrackerRecHits'),
                         ('B0TrackerMeasurements_hits', 'B0TrackerRecHits'),
                         (p+'Trajectories_seed', 'B0TrackerSeeds'),
                         (p+'Tracks_trajectory', p+'Trajectories'),
                         (p+'Tracks_measurements', 'B0TrackerMeasurements'),
                         (p+'TrackAssociations_rec', p+'Tracks'),
                         (p+'TrackAssociations_sim', 'MCParticles'),
                         ('B0TrackerHits_particle', 'MCParticles')]:
    indices = a('_' + relation + '.index')
    collids = a('_' + relation + '.collectionID')
    assert all(i < 0 or c == ids[target] for ii, cc in zip(indices, collids) for i, c in zip(ii, cc)), relation

stats = Counter({name:0 for name in ['seeds','tracks','fake_tracks_purity_below_half','duplicates_all_truth_particles','duplicates_for_accepted_primary_protons','wrong_charge_quality_selected','tails_absolute_qop_error_above_20_percent']})
measurements = Counter()
bad_covariances = []
minimum_correlation_eigenvalue = None
missing_first, missing_surface, unmatched_surface, accepted, residuals, paired_rows = [], [], [], [], [], []
for ev, key in enumerate(keys):
    start_hits = {}
    for s, pi in enumerate(seed_param[ev]):
        stats['seeds'] += 1
        cov = np.zeros((6,6))
        k = 0
        for i in range(6):
            for j in range(i+1):
                cov[i,j] = cov[j,i] = seed_cov[ev][pi][k]
                k += 1
        values = np.array([field[ev][pi] for field in seed_values])
        if not np.isfinite(values).all() or not np.isfinite(cov).all() or np.any(np.diag(cov)<=0):
            bad_covariances.append([*key,s,'nonfinite_or_nonpositive_diagonal'])
        else:
            sigma = np.sqrt(np.diag(cov))
            smallest = float(np.linalg.eigvalsh(cov/np.outer(sigma,sigma))[0])
            minimum_correlation_eigenvalue = smallest if minimum_correlation_eigenvalue is None else min(minimum_correlation_eigenvalue,smallest)
            if smallest <= 0:
                bad_covariances.append([*key,s,smallest])
        surf = surface[ev][pi]
        if surf == 0:
            stats['legacy_origin_seeds'] += 1
            continue
        stats['sensor_bound_seeds'] += 1
        hits = set(seed_hits[ev][seed_b[ev][s]:seed_e[ev][s]])
        on_start = set()
        for m, ms in enumerate(meas_surface[ev]):
            if ms == surf:
                on_start.update(hits & set(meas_hits[ev][meas_b[ev][m]:meas_e[ev][m]]))
        if len(on_start) != 1:
            unmatched_surface.append([*key, s, surf, sorted(on_start)])
        else:
            start_hits[s] = next(iter(on_start))
    match = {}
    for tr, q in enumerate(track_q[ev]):
        stats['tracks'] += 1
        ms = track_meas[ev][track_b[ev][tr]:track_e[ev][tr]]
        measurements[len(ms)] += 1
        trajectory = track_traj[ev][tr]
        assert 0 <= trajectory < len(traj_seed[ev]), 'Dangling track trajectory'
        seed = traj_seed[ev][trajectory]
        assert 0 <= seed < len(seed_param[ev]), 'Dangling trajectory seed'
        hs = set()
        for m in ms:
            assert 0 <= m < len(meas_b[ev]), 'Dangling measurement'
            hs.update(meas_hits[ev][meas_b[ev][m]:meas_e[ev][m]])
        if seed in start_hits:
            stats['tracks_with_sensor_seed'] += 1
            start_surface = surface[ev][seed_param[ev][seed]]
            if any(meas_surface[ev][m] == start_surface for m in ms):
                stats['tracks_containing_start_surface'] += 1
            else:
                missing_surface.append([*key, tr, seed, start_surface])
            if start_hits[seed] not in hs:
                missing_first.append([*key, tr, seed, start_hits[seed], len(ms)])
            else:
                stats['tracks_containing_start_hit'] += 1
        assocs = [(w, mi) for ri, mi, w in zip(rec[ev], sim[ev], weight[ev]) if ri == tr]
        if not assocs or max(w for w, _ in assocs) < 0.5:
            stats['fake_tracks_purity_below_half'] += 1
            continue
        _, mc = max(assocs, key=lambda x: (x[0], -x[1]))
        match.setdefault(mc, []).append(tr)
    stats['duplicates_all_truth_particles'] += sum(len(tracks)-1 for tracks in match.values())
    for mc, (status, pdg) in enumerate(zip(mc_status[ev], mc_pdg[ev])):
        if status != 1 or pdg != 2212:
            continue
        truth_p = np.linalg.norm([c[ev][mc] for c in mc_p])
        clean_z = []
        for h, particle in enumerate(hit_part[ev]):
            if particle != mc:
                continue
            ph = np.linalg.norm([c[ev][h] for c in hit_mom])
            if ph <= 0.5*truth_p or hit_mom[2][ev][h] <= 0:
                continue
            clean_z.append(hit_pos[0][ev][h]*np.sin(args.crossing_angle) + hit_pos[2][ev][h]*np.cos(args.crossing_angle))
        z = np.sort(clean_z)
        stations = int(len(z)>0) + int(np.count_nonzero(np.diff(z) >= 50))
        if stations < 3:
            continue
        accepted.append([*key, mc])
        candidates = match.get(mc, [])
        if not candidates:
            paired_rows.append({'key': [*key, mc], 'success': False})
            continue
        stats['accepted_primary_protons_matched'] += 1
        stats['duplicates_for_accepted_primary_protons'] += len(candidates)-1
        tr = min(candidates, key=lambda i: (-(track_e[ev][i]-track_b[ev][i]), track_chi2[ev][i], i))
        reco_p = np.linalg.norm([c[ev][tr] for c in track_p])
        truth_qop = mc_charge[ev][mc]/truth_p
        reco_qop = track_q[ev][tr]/reco_p
        residual = 100*(reco_qop/truth_qop-1)
        if np.sign(reco_qop) != np.sign(truth_qop):
            stats['wrong_charge_quality_selected'] += 1
        if abs(residual) > 20:
            stats['tails_absolute_qop_error_above_20_percent'] += 1
        residuals.append(residual)
        paired_rows.append({'key': [*key, mc], 'success': True, 'track': tr,
                            'nmeasurements': track_e[ev][tr]-track_b[ev][tr],
                            'qop_residual_percent': residual})

r = np.asarray(residuals)
result = {'valid': not (unmatched_surface or bad_covariances or missing_surface),
          'provenance': {'crossing_angle_rad': args.crossing_angle, 'prefix': args.prefix,
                         'percentile_method': 'numpy linear', 'association_minimum_purity': 0.5,
                         'minimum_stations': 3, 'station_separation_mm': 50},
          'input': args.file, 'events': len(keys), 'statistics': dict(stats),
          'track_measurement_counts': dict(measurements),
          'invalid_seed_surface_matches': unmatched_surface,
          'invalid_seed_covariances': bad_covariances,
          'minimum_seed_correlation_eigenvalue': minimum_correlation_eigenvalue,
          'tracks_missing_start_hit': missing_first,
          'tracks_missing_start_surface': missing_surface,
          'denominator': 'Primary protons with >=3 stations separated by >=50 mm in outgoing-ion z, each retaining >half origin momentum and positive lab pz; match purity >=0.5; quality-selected final track.',
          'accepted_primary_protons': len(accepted),
          'qop_resolution': {'successes': len(r), 'half68_percent': float((np.percentile(r,84)-np.percentile(r,16))/2) if len(r) else None,
                             'median_percent': float(np.median(r)) if len(r) else None},
          'per_particle': paired_rows}
print(json.dumps(result, indent=2))

sys.exit(0 if result["valid"] else 1)
