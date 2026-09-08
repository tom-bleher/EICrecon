"""B0 final-track benchmark; run inside eic-shell (numpy, awkward, uproot).

Matching adapted from presentation17/scripts/b0_truth.py, retaining its
fractional direct-primary convention and relation/charge-closure checks.

Direct-primary truth and fractional measurement-count matching for B0 plots.

A measurement is one vote, shared equally between distinct threshold-passing linked SimTrackerHits.
Unknown raw-hit contributions receive an unknown share. Secondary, non-primary,
and unknown shares never count as direct primary signal, but remain in the
measurement denominator. Persisted association weights are not charge fractions.
"""
import math
from collections import defaultdict
from fractions import Fraction

SECONDARY_BIT = 1 << 30
MIN_PURITY = 0.5
MATCHING_DEFINITION = (
    'Fractional measurement-count direct-primary purity: each fitted measurement '
    'has one vote split equally among distinct threshold-passing SimTrackerHits '
    '(plus one unknown share per unlinked raw hit). Only quality bit30-clear '
    'hits linked to generatorStatus=1 PDG=2212 count for that primary; all '
    'secondary, other-particle and unknown shares remain in the denominator. '
    'This is not energy/charge purity; association and measurement weights are '
    'not used. Active contributors pass the explicitly supplied digitization threshold; '
    'their rounded charges must reproduce each linked raw hit. Purity >= 0.5. Each track chooses highest direct purity, ties '
    'lowest MC index. Each primary chooses highest direct purity, most fitted '
    'measurements, lowest chi2, then track index. No truth-nearest matching.'
)
ACCEPTANCE_DEFINITION = (
    'Generated status=1 PDG=2212 proton with quality bit30-clear simulated hits '
    'in at least three distinct physical B0 stations; no reconstructed-hit, '
    'hit-momentum or direction requirement.'
)
COMMON_FIELDS = [
    'MCParticles.generatorStatus', 'MCParticles.PDG',
    'B0TrackerHits.cellID', 'B0TrackerHits.quality', 'B0TrackerHits.eDep',
    '_B0TrackerHits_particle.index', '_B0TrackerHits_particle.collectionID',
    'B0TrackerRawHits.cellID', 'B0TrackerRawHits.charge', 'B0TrackerRecHits.cellID',
    '_B0TrackerRecHits_rawHit.index', '_B0TrackerRecHits_rawHit.collectionID',
    '_B0TrackerRawHitAssociations_rawHit.index',
    '_B0TrackerRawHitAssociations_rawHit.collectionID',
    '_B0TrackerRawHitAssociations_simHit.index',
    '_B0TrackerRawHitAssociations_simHit.collectionID',
    'B0TrackerMeasurements.hits_begin', 'B0TrackerMeasurements.hits_end',
    '_B0TrackerMeasurements_hits.index', '_B0TrackerMeasurements_hits.collectionID',
]


def track_fields(prefix):
    return [prefix + 'Tracks.' + name for name in
            ['chi2', 'measurements_begin', 'measurements_end']] + [
        '_' + prefix + 'Tracks_measurements.' + name for name in ['index', 'collectionID']]


def direct_primary(status, pdg, quality):
    return status == 1 and pdg == 2212 and not (int(quality) & SECONDARY_BIT)


def station(cell, layers_per_station):
    layer = (int(cell) >> 8) & 15
    maximum = 4 * layers_per_station
    if not 1 <= layer <= maximum:
        raise ValueError(f'Invalid B0 layer {layer} for {layers_per_station} layers/station')
    return (layer - 1) // layers_per_station + 1


def parallel(*arrays):
    if len({len(a) for a in arrays}) != 1:
        raise ValueError('Relation/data arrays have inconsistent lengths')
    return zip(*arrays)


def checked_index(index, collection, expected, size):
    if collection != expected or not 0 <= index < size:
        raise ValueError(f'Invalid relation ({collection},{index}); expected collection {expected}, size {size}')
    return index


def checked_slice(begin, end, size):
    if not 0 <= begin <= end <= size:
        raise ValueError(f'Invalid relation slice [{begin}:{end}], size {size}')
    return range(begin, end)


def primary_station_sets(get, ids, layers_per_station):
    """Return direct and ancestry station sets, including primaries with no hits."""
    status, pdg = get('MCParticles.generatorStatus'), get('MCParticles.PDG')
    primaries = {i for i, (s, p) in enumerate(parallel(status, pdg)) if s == 1 and p == 2212}
    direct, ancestry = ({i: set() for i in primaries} for _ in range(2))
    owners = []
    for cell, quality, owner, collection in parallel(
            get('B0TrackerHits.cellID'), get('B0TrackerHits.quality'),
            get('_B0TrackerHits_particle.index'), get('_B0TrackerHits_particle.collectionID')):
        checked_index(owner, collection, ids['MCParticles'], len(status))
        st = station(cell, layers_per_station)
        primary = owner if direct_primary(status[owner], pdg[owner], quality) else None
        owners.append(primary)
        if owner in primaries:
            ancestry[owner].add(st)
        if primary is not None:
            direct[primary].add(st)
    return direct, ancestry, owners


class EventTruth:
    def __init__(self, get, ids, layers_per_station, threshold_keV):
        self.get, self.ids = get, ids
        if not math.isfinite(threshold_keV) or threshold_keV <= 0:
            raise ValueError('Digitization threshold in keV must be finite and positive')
        self.threshold_keV = threshold_keV
        self.direct_stations, self.ancestry_stations, self.sim_primary = primary_station_sets(get, ids, layers_per_station)
        self.primaries = set(self.direct_stations)
        cells = get('B0TrackerHits.cellID')
        raw_cells, rec_cells = get('B0TrackerRawHits.cellID'), get('B0TrackerRecHits.cellID')
        self.raw_sim_all = defaultdict(set)
        for raw, rc, sim, sc in parallel(
                get('_B0TrackerRawHitAssociations_rawHit.index'),
                get('_B0TrackerRawHitAssociations_rawHit.collectionID'),
                get('_B0TrackerRawHitAssociations_simHit.index'),
                get('_B0TrackerRawHitAssociations_simHit.collectionID')):
            checked_index(raw, rc, ids['B0TrackerRawHits'], len(raw_cells))
            checked_index(sim, sc, ids['B0TrackerHits'], len(cells))
            if raw_cells[raw] != cells[sim]:
                raise ValueError('Raw/sim association crosses cell IDs')
            self.raw_sim_all[raw].add(sim)
        energy, charge = get('B0TrackerHits.eDep'), get('B0TrackerRawHits.charge')
        list(parallel(energy, cells))
        list(parallel(charge, raw_cells))
        if any(not math.isfinite(e) or e < 0 for e in energy):
            raise ValueError('Simulated energy deposits must be finite and nonnegative')
        self.raw_sim = defaultdict(set)
        self.inactive_links = 0
        self.unknown_raw_hits = len(raw_cells) - len(self.raw_sim_all)
        for raw, sims in self.raw_sim_all.items():
            accepted = {sim for sim in sims if energy[sim] >= threshold_keV * 1e-6}
            # SiliconTrackerDigi rounds each accepted deposit before summing it.
            expected = sum(math.floor(energy[sim] * 1e6 + .5) for sim in accepted)
            if expected != charge[raw]:
                raise ValueError(f'Raw charge closure failed at raw hit {raw}: '
                                 f'expected {expected}, stored {charge[raw]}, threshold {threshold_keV} keV')
            self.raw_sim[raw] = accepted
            self.inactive_links += len(sims) - len(accepted)
        self.rec_raw = []
        for cell, raw, collection in parallel(rec_cells, get('_B0TrackerRecHits_rawHit.index'),
                                               get('_B0TrackerRecHits_rawHit.collectionID')):
            checked_index(raw, collection, ids['B0TrackerRawHits'], len(raw_cells))
            if cell != raw_cells[raw]:
                raise ValueError('Reconstructed/raw hit relation crosses cell IDs')
            self.rec_raw.append(raw)
        self.measurement_shares = []
        links, collections = get('_B0TrackerMeasurements_hits.index'), get('_B0TrackerMeasurements_hits.collectionID')
        list(parallel(links, collections))
        for begin, end in parallel(get('B0TrackerMeasurements.hits_begin'), get('B0TrackerMeasurements.hits_end')):
            tokens = set()
            for i in checked_slice(begin, end, len(links)):
                rec = checked_index(links[i], collections[i], ids['B0TrackerRecHits'], len(rec_cells))
                raw = self.rec_raw[rec]
                sims = self.raw_sim[raw]
                tokens.update(('sim', sim) for sim in sims)
                if not sims:
                    tokens.add(('unknown_raw', raw))
            shares = defaultdict(Fraction)
            for kind, index in sorted(tokens):
                primary = self.sim_primary[index] if kind == 'sim' else None
                if primary is not None:
                    shares[primary] += Fraction(1, len(tokens))
            self.measurement_shares.append(dict(shares))

    def match(self, prefix, minimum=MIN_PURITY):
        if not 0 < minimum <= 1:
            raise ValueError('Purity threshold must be in (0,1]')
        g, ids = self.get, self.ids
        links, collections = g('_' + prefix + 'Tracks_measurements.index'), g('_' + prefix + 'Tracks_measurements.collectionID')
        list(parallel(links, collections))
        choices = {}
        tracks = []
        for tr, (begin, end, chi2) in enumerate(parallel(
                g(prefix + 'Tracks.measurements_begin'), g(prefix + 'Tracks.measurements_end'), g(prefix + 'Tracks.chi2'))):
            if not math.isfinite(chi2):
                raise ValueError('Track chi2 is nonfinite')
            measurements = [checked_index(links[i], collections[i], ids['B0TrackerMeasurements'], len(self.measurement_shares))
                            for i in checked_slice(begin, end, len(links))]
            if len(set(measurements)) != len(measurements):
                raise ValueError('Track repeats a fitted measurement')
            votes = defaultdict(Fraction)
            for measurement in measurements:
                for primary, fraction in self.measurement_shares[measurement].items():
                    votes[primary] += fraction
            fractions = {mc: count / len(measurements) for mc, count in votes.items()}
            row = dict(track=tr, measurements=sorted(measurements), nmeasurements=len(measurements), chi2=chi2, fractions={mc: float(value) for mc, value in fractions.items()},
                       primary=None, purity=0.0)
            if fractions:
                mc = min(fractions, key=lambda p: (-fractions[p], p))
                row.update(primary=mc, purity=float(fractions[mc]))
                if fractions[mc] >= minimum:
                    rank = (-fractions[mc], -len(measurements), chi2, tr)
                    if mc not in choices or rank < choices[mc][0]:
                        choices[mc] = (rank, row)
            tracks.append(row)
        return {mc: row for mc, (_, row) in choices.items()}, tracks


def perigee(vertex, momentum, charge):
    """Origin/lab-z perigee of a ray; assumes field-free truth-vertex region."""
    import numpy as np
    v, p = np.asarray(vertex, dtype=float), np.asarray(momentum, dtype=float)
    pt2, mag = p[:2] @ p[:2], np.linalg.norm(p)
    if not np.isfinite([*v, *p, charge]).all() or pt2 <= 0 or mag <= 0 or charge == 0:
        raise ValueError('Truth ray needs finite vertex/momentum, nonzero pt and charge')
    point = v - (v[:2] @ p[:2]) / pt2 * p
    phi = math.atan2(p[1], p[0])
    return np.array([-point[0] * math.sin(phi) + point[1] * math.cos(phi), point[2],
                     phi, math.atan2(math.sqrt(pt2), p[2]), charge / mag])


def covariance_metrics(packed, fitted, truth=None):
    """Never clamp/remove invalid matrices. Report their status explicitly."""
    import numpy as np
    if len(packed) != 21 or len(fitted) != 6:
        raise ValueError('Expected packed Cov6f and six bound parameters')
    cov = np.array([[packed[min(i, j) + max(i, j) * (max(i, j) + 1) // 2]
                     for j in range(6)] for i in range(6)])
    result = {'covariance_status': 'nonfinite', 'pulls': None, 'joint95': None,
              'mahalanobis2': None, 'qop_residual_percent': None}
    if truth is not None and math.isfinite(fitted[4]):
        result['qop_residual_percent'] = float(100 * (fitted[4] / truth[4] - 1))
    if not np.isfinite(cov).all() or not np.isfinite(fitted).all():
        return result
    result['covariance_status'] = 'nonpositive'
    if np.any(np.diag(cov) <= 0):
        return result
    scale = np.sqrt(np.diag(cov))
    correlation = cov / np.outer(scale, scale)
    eigenvalues = np.linalg.eigvalsh(correlation)
    result['min_normalized_eigenvalue'] = float(eigenvalues[0])
    if eigenvalues[0] <= 0:
        # Bound storage rounding element by element rather than rejecting small
        # positive eigenvalues with a global, arbitrary condition-number cut.
        stored = cov.astype(np.float32)
        above = np.nextafter(stored, np.float32(np.inf)).astype(float)
        below = np.nextafter(stored, np.float32(-np.inf)).astype(float)
        upward, downward = above - cov, cov - below
        upward = np.where(np.isfinite(upward), upward, downward)
        downward = np.where(np.isfinite(downward), downward, upward)
        bound = float(np.linalg.norm(.5 * np.maximum(upward, downward) /
                                     np.outer(scale, scale), 2))
        result['rounding_eigenvalue_bound'] = bound
        if eigenvalues[0] >= -bound:
            result['covariance_status'] = 'roundoff_limited'
        return result
    result['covariance_status'] = 'positive_definite'
    if truth is not None:
        residual = np.asarray(fitted[:5]) - truth
        residual[2] = math.atan2(math.sin(residual[2]), math.cos(residual[2]))
        pulls = residual / scale[:5]
        distance = float(pulls @ np.linalg.solve(correlation[:5, :5], pulls))
        result.update(pulls=pulls.tolist(), mahalanobis2=distance, joint95=distance <= 11.070497693516351,
                      qop_residual_percent=float(100 * residual[4] / truth[4]))
    return result


def fitted_parameters(get, ids, prefix, track, seed_sizes):
    tr = prefix + 'Tracks'
    trajectory = checked_index(get('_' + tr + '_trajectory.index')[track],
                               get('_' + tr + '_trajectory.collectionID')[track],
                               ids[prefix + 'Trajectories'], len(get(prefix + 'Trajectories.trackParameters_begin')))
    begin, end = (get(prefix + 'Trajectories.trackParameters_' + x)[trajectory] for x in ('begin', 'end'))
    refs = '_' + prefix + 'Trajectories_trackParameters'
    indices = get(refs + '.index')
    selected = list(checked_slice(begin, end, len(indices)))
    if len(selected) != 1:
        raise ValueError('Expected exactly one fitted perigee parameter per trajectory')
    i = selected[0]
    par = prefix + 'TrackParameters'
    index = checked_index(indices[i], get(refs + '.collectionID')[i], ids[par], len(get(par + '.qOverP')))
    if get(par + '.surface')[index] != 0:
        raise ValueError('Expected origin perigee (surface ID zero) from CKFTracking')
    state = [get(par + '.' + field)[index] for field in ('loc.a', 'loc.b', 'phi', 'theta', 'qOverP', 'time')]
    seed_refs = '_' + prefix + 'Trajectories_seed'
    seed_id = get(seed_refs + '.collectionID')[trajectory]
    seed_index = get(seed_refs + '.index')[trajectory]
    checked_index(seed_index, seed_id, seed_id, seed_sizes.get(seed_id, 0))
    return state, get(par + '.covariance.covariance[21]')[index], [seed_id, seed_index]


def collection_alias(name, prefix, suffix):
    """Map internal fitted-collection names to filtered or unfiltered PODIO names."""
    for kind in ('Tracks', 'Trajectories', 'TrackParameters'):
        for leading in ('', '_'):
            stem = leading + prefix + kind
            if name == stem or name.startswith(stem + '.') or name.startswith(stem + '_'):
                return stem + suffix + name[len(stem):]
    return name


def analyze(path, prefix, layers_per_station, threshold_keV, crossing_angle, suffix=''):
    import hashlib
    import json
    import awkward as ak
    import numpy as np
    import uproot
    fields = COMMON_FIELDS + track_fields(prefix)
    fields += ['EventHeader.runNumber', 'EventHeader.eventNumber', 'MCParticles.charge']
    fields += ['MCParticles.' + name + '.' + c for name in ('vertex', 'momentum') for c in 'xyz']
    fields += ['B0TrackerMeasurements.' + f for f in ('surface', 'loc.a', 'loc.b')]
    for relation in ('Tracks_trajectory', 'Trajectories_trackParameters', 'Trajectories_seed'):
        fields += ['_' + prefix + relation + '.' + c for c in ('index', 'collectionID')]
    fields += [prefix + 'Trajectories.trackParameters_' + c for c in ('begin', 'end')]
    fields += [prefix + 'TrackParameters.' + c for c in
               ('surface', 'loc.a', 'loc.b', 'phi', 'theta', 'qOverP', 'time', 'covariance.covariance[21]')]
    fields = sorted(set(fields))
    aliases = {name: collection_alias(name, prefix, suffix) for name in fields}
    rows, track_rows, events, hashes = [], [], [], {}
    with uproot.open(path) as root:
        metadata = ak.to_list(root['podio_metadata']['events___CollectionTypeInfo'].arrays())
        if len(metadata) != 1:
            raise ValueError('Expected one collection metadata entry')
        md = metadata[0]
        ids = dict(parallel(md['events___CollectionTypeInfo.name'], md['events___CollectionTypeInfo.collectionID']))
        if len(ids) != len(md['events___CollectionTypeInfo.name']) or len(ids) != len(set(ids.values())):
            raise ValueError('Duplicate collection identities')
        seed_fields = {cid: name + '.quality' for name, cid, dtype, subset in parallel(
            md['events___CollectionTypeInfo.name'], md['events___CollectionTypeInfo.collectionID'],
            md['events___CollectionTypeInfo.dataType'], md['events___CollectionTypeInfo.isSubset'])
            if dtype == 'edm4eic::TrackSeedCollection' and not subset}
        for kind in ('Tracks', 'Trajectories', 'TrackParameters'):
            canonical = prefix + kind
            ids[canonical] = ids[collection_alias(canonical, prefix, suffix)]
        for batch in root['events'].iterate(list(aliases.values()) + list(seed_fields.values()), step_size=128):
            arrays = {name: ak.to_list(batch[actual]) for name, actual in aliases.items()}
            for ev in range(len(batch)):
                get = lambda name: arrays[name][ev]
                runs, numbers = get('EventHeader.runNumber'), get('EventHeader.eventNumber')
                if len(runs) != 1 or len(numbers) != 1:
                    raise ValueError('Expected one EventHeader')
                key = [runs[0], numbers[0]]
                event_key = str(key)
                if event_key in hashes:
                    raise ValueError('Duplicate event identity')
                events.append(key)
                upstream = {f: get(f) for f in fields if not f.lstrip('_').startswith(prefix)}
                hashes[event_key] = hashlib.sha256(json.dumps(upstream, sort_keys=True).encode()).hexdigest()
                truth = EventTruth(get, ids, layers_per_station, threshold_keV)
                matches, tracks = truth.match(prefix)
                truth_states = {}
                for mc in sorted(truth.primaries):
                    momentum = [get('MCParticles.momentum.' + c)[mc] for c in 'xyz']
                    vertex = [get('MCParticles.vertex.' + c)[mc] for c in 'xyz']
                    state = perigee(vertex, momentum, get('MCParticles.charge')[mc])
                    truth_states[mc] = state
                    ca, sa = math.cos(crossing_angle), math.sin(crossing_angle)
                    x, y, z = momentum
                    rows.append(dict(key=key + [mc], p=float(np.linalg.norm(momentum)),
                                     ion_theta_mrad=1000 * math.atan2(math.hypot(x * ca - z * sa, y), x * sa + z * ca),
                                     stations=len(truth.direct_stations[mc]), truth=state.tolist(),
                                     matched=mc in matches, track=None))
                by_primary = {r['key'][2]: r for r in rows[-len(truth.primaries):]} if truth.primaries else {}
                for track in tracks:
                    seed_sizes = {cid: len(batch[field][ev]) for cid, field in seed_fields.items()}
                    state, cov, seed = fitted_parameters(get, ids, prefix, track['track'], seed_sizes)
                    primary = track['primary'] if track['purity'] >= MIN_PURITY else None
                    audit = covariance_metrics(cov, state, truth_states.get(primary))
                    row = dict(key=key + [track['track']], seed=seed, primary=primary, purity=track['purity'],
                               measurements=track['measurements'], chi2=track['chi2'], **audit)
                    # Keep invalid fit diagnostics JSON-safe without changing acceptance.
                    row['parameters'] = [float(x) if math.isfinite(x) else None for x in state]
                    row['bound_covariance'] = [float(x) if math.isfinite(x) else None for x in cov]
                    track_rows.append(row)
                    if primary in matches and matches[primary]['track'] == track['track']:
                        by_primary[primary]['track'] = row
    return dict(input=str(path), events=events, upstream_hashes=hashes, rows=rows, tracks=track_rows)


def summarize(rows, tracks=None):
    import numpy as np
    accepted = [r for r in rows if r['stations'] >= 3]
    matched = [r['track'] for r in accepted if r['matched']]
    valid = [t for t in matched if t['joint95'] is not None]
    def distribution(values):
        return dict(zip(('p16', 'median', 'p84'), np.percentile(values, [16, 50, 84]).tolist())) if values else None
    out = dict(eligible=len(accepted), matched=len(matched), efficiency=len(matched) / len(accepted) if accepted else None,
               coverage_valid=len(valid), coverage_invalid=len(matched) - len(valid),
               joint95_coverage=sum(t['joint95'] for t in valid) / len(valid) if valid else None,
               pulls=[distribution([t['pulls'][i] for t in valid]) for i in range(5)],
               qop_residual_valid=sum(t['qop_residual_percent'] is not None for t in matched),
               qop_residual_percent=distribution([t['qop_residual_percent'] for t in matched
                                                   if t['qop_residual_percent'] is not None]))
    if tracks is not None:
        counts = defaultdict(int)
        for tr in tracks:
            counts[tr['covariance_status']] += 1
        direct = [tr for tr in tracks if tr['primary'] is not None]
        unique = {(tuple(tr['key'][:2]), tr['primary']) for tr in direct}
        out.update(tracks=len(tracks), direct_primary_unmatched_tracks=len(tracks)-len(direct),
                   direct_primary_fake_fraction=(len(tracks)-len(direct))/len(tracks) if tracks else None,
                   duplicate_tracks=len(direct)-len(unique), covariance_status=dict(counts))
    return out


def compare_candidates(baseline, variant):
    """Compare unique pre-ambiguity candidates, never choose a truth-nearest fit."""
    import numpy as np
    def group(rows):
        result = defaultdict(list)
        for row in rows:
            key = (tuple(row['key'][:2]), tuple(row['seed']), tuple(sorted(row['measurements'])))
            result[key].append(row)
        return result
    left, right = group(baseline), group(variant)
    common = left.keys() & right.keys()
    unique = sorted(k for k in common if len(left[k]) == len(right[k]) == 1)
    ambiguous = sorted(k for k in common if len(left[k]) != 1 or len(right[k]) != 1)
    shifts, ratios = [], []
    invalid = 0
    def unpack(packed):
        return np.array([[packed[min(i, j) + max(i, j) * (max(i, j) + 1) // 2]
                          for j in range(5)] for i in range(5)], dtype=float)
    for key in unique:
        a, b = left[key][0], right[key][0]
        ac, bc = unpack(a['bound_covariance']), unpack(b['bound_covariance'])
        ap, bp = np.array(a['parameters'][:5], dtype=float), np.array(b['parameters'][:5], dtype=float)
        if (not np.isfinite([ac, bc]).all() or not np.isfinite([ap, bp]).all()
                or np.any(np.diag(ac) <= 0) or np.any(np.diag(bc) <= 0)):
            invalid += 1
            continue
        delta = bp - ap
        delta[2] = math.atan2(math.sin(delta[2]), math.cos(delta[2]))
        shifts.append(np.abs(delta) / np.sqrt(np.diag(ac)))
        ratios.append(np.sqrt(np.diag(bc) / np.diag(ac)))
    def distributions(values):
        if not values:
            return None
        percentiles = np.percentile(values, [0, 16, 50, 84, 100], axis=0)
        return {name: dict(zip(('min', 'p16', 'median', 'p84', 'max'), percentiles[:, i].tolist()))
                for i, name in enumerate(('loc0', 'loc1', 'phi', 'theta', 'qOverP'))}
    return dict(baseline_candidates=len(baseline), variant_candidates=len(variant),
                gain_groups=sorted(right.keys() - left.keys()), loss_groups=sorted(left.keys() - right.keys()),
                common_groups=len(common), ambiguous_groups=ambiguous, unique_pairs=len(unique),
                numerical_pairs=len(shifts), invalid_pairs=invalid,
                abs_parameter_shift_in_baseline_sigma=distributions(shifts),
                sigma_ratio_variant_over_baseline=distributions(ratios))


def main():
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('--variant', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--prefix', default='B0TrackerCKF')
    parser.add_argument('--suffix', choices=('', 'Unfiltered'), default='',
                        help='Analyze candidates before ambiguity resolution with Unfiltered')
    parser.add_argument('--layers-per-station', type=int, choices=(1, 2), default=2,
                        help='Current B0_tracker.xml: layer=2*(station-1)+(1=back,2=front)')
    parser.add_argument('--threshold-keV', type=float, required=True)
    parser.add_argument('--crossing-angle', type=float, default=-0.025, help='Outgoing ion-axis angle in lab, radians')
    parser.add_argument('--provenance', type=Path, help='JSON with software, geometry, map, sample/generator and beam settings')
    args = parser.parse_args()
    results = {}
    for name, path in [('baseline', args.baseline), ('variant', args.variant)]:
        if path is None:
            continue
        result = analyze(path, args.prefix, args.layers_per_station, args.threshold_keV,
                         args.crossing_angle, args.suffix)
        result['summary'] = summarize(result['rows'], result['tracks'])
        result['bins'] = {}
        for field, edges in [('p', [0, 8, 15, 25, 35, 45, 1000]),
                             ('ion_theta_mrad', [0, 4, 8, 12, 16, 22, 30, 1000]), ('stations', [0, 1, 2, 3, 4, 5])]:
            result['bins'][field] = [dict(low=lo, high=hi, **summarize([r for r in result['rows'] if lo <= r[field] < hi]))
                                     for lo, hi in zip(edges, edges[1:])]
        results[name] = result
    if args.variant:
        b, v = results['baseline'], results['variant']
        if b['upstream_hashes'] != v['upstream_hashes']:
            raise ValueError('Paired comparison requires identical event identities and upstream input fields')
        br, vr = ({tuple(r['key']): r for r in x['rows']} for x in (b, v))
        if br.keys() != vr.keys():
            raise ValueError('Paired truth identities differ')
        common = [k for k in br if br[k]['matched'] and vr[k]['matched']]
        equal_hits = [k for k in common if br[k]['track']['measurements'] == vr[k]['track']['measurements']]
        results['paired'] = dict(gains=[k for k in br if not br[k]['matched'] and vr[k]['matched']],
                                 losses=[k for k in br if br[k]['matched'] and not vr[k]['matched']],
                                 common=len(common), equal_assigned_hits=len(equal_hits), equal_hit_keys=equal_hits,
                                 baseline_equal_hits=summarize([br[k] for k in equal_hits]),
                                 variant_equal_hits=summarize([vr[k] for k in equal_hits]))
        results['paired_candidates'] = compare_candidates(b['tracks'], v['tracks'])
    results['definitions'] = dict(acceptance=ACCEPTANCE_DEFINITION, matching=MATCHING_DEFINITION,
        covariance='Full 5D origin-perigee pulls; chi-square(5) 95% threshold 11.070497693516351. Float32 normalized eigenvalue guard; invalid retained in counts.',
        fake='Track without >=0.5 direct-primary-proton purity, not generic charged-track fake rate.',
        precision='JSON full floating-point precision; percentiles NumPy linear; units mm,rad,GeV; q/p residual percent.',
        options={k: str(v) if isinstance(v, Path) else v for k,v in vars(args).items()})
    if args.provenance:
        results['provenance'] = json.loads(args.provenance.read_text())
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v['summary'] for k, v in results.items() if k in ('baseline', 'variant')}, indent=2))


if __name__ == '__main__':
    main()
