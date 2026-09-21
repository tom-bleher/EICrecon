"""Adapted from presentation17/scripts/test_b0_truth.py.

Focused synthetic checks for truth ancestry, PODIO relations and matching."""
import unittest
from analyze import EventTruth, SECONDARY_BIT, primary_station_sets, track_fields

P = 'B0TrackerCKF'
IDS = dict(MCParticles=1, B0TrackerHits=2, B0TrackerRawHits=3,
           B0TrackerRecHits=4, B0TrackerMeasurements=5)


def event(sim=((0, 0),), raw=((0,),), measurements=((0,),), tracks=((0,),), pdg=(2212, 2212)):
    cell = 150 | (1 << 8)
    d = {'MCParticles.generatorStatus': [1] * len(pdg), 'MCParticles.PDG': list(pdg),
         'B0TrackerHits.cellID': [cell] * len(sim), 'B0TrackerHits.eDep': [10e-6] * len(sim),
         'B0TrackerHits.quality': [q for _, q in sim],
         '_B0TrackerHits_particle.index': [m for m, _ in sim],
         '_B0TrackerHits_particle.collectionID': [IDS['MCParticles']] * len(sim),
         'B0TrackerRawHits.cellID': [cell] * len(raw),
         'B0TrackerRawHits.charge': [10 * len(set(sims)) for sims in raw],
         'B0TrackerRecHits.cellID': [cell] * len(raw),
         '_B0TrackerRecHits_rawHit.index': list(range(len(raw))),
         '_B0TrackerRecHits_rawHit.collectionID': [IDS['B0TrackerRawHits']] * len(raw),
         '_B0TrackerRawHitAssociations_rawHit.index': [r for r, sims in enumerate(raw) for _ in sims],
         '_B0TrackerRawHitAssociations_rawHit.collectionID': [IDS['B0TrackerRawHits']] * sum(map(len, raw)),
         '_B0TrackerRawHitAssociations_simHit.index': [s for sims in raw for s in sims],
         '_B0TrackerRawHitAssociations_simHit.collectionID': [IDS['B0TrackerHits']] * sum(map(len, raw)),
         P + 'Tracks.chi2': [1.] * len(tracks)}
    for prefix, relation, groups, target in [
            ('B0TrackerMeasurements', 'hits', measurements, 'B0TrackerRecHits'),
            (P + 'Tracks', 'measurements', tracks, 'B0TrackerMeasurements')]:
        start = 0
        d[prefix + '.' + relation + '_begin'] = []
        d[prefix + '.' + relation + '_end'] = []
        for group in groups:
            d[prefix + '.' + relation + '_begin'].append(start)
            start += len(group)
            d[prefix + '.' + relation + '_end'].append(start)
        d['_' + prefix + '_' + relation + '.index'] = [i for g in groups for i in g]
        d['_' + prefix + '_' + relation + '.collectionID'] = [IDS[target]] * start
    return d


def match(d):
    return EventTruth(d.__getitem__, IDS, 1, 10).match(P)


class TruthTests(unittest.TestCase):
    def test_inactive_links_excluded_and_charge_verified(self):
        d = event(sim=((0, 0), (0, SECONDARY_BIT)), raw=((0, 1),))
        d['B0TrackerHits.eDep'][1] = 9e-6
        d['B0TrackerRawHits.charge'][0] = 10
        self.assertEqual(match(d)[0][0]['purity'], 1.)
        d['B0TrackerRawHits.charge'][0] = 19
        with self.assertRaisesRegex(ValueError, 'charge closure'):
            match(d)

    def test_wrong_threshold_and_invalid_energy_fail(self):
        d = event()
        with self.assertRaisesRegex(ValueError, 'charge closure'):
            EventTruth(d.__getitem__, IDS, 1, 11)
        d['B0TrackerHits.eDep'][0] = float('nan')
        with self.assertRaisesRegex(ValueError, 'finite'):
            match(d)

    def test_charge_rounding_matches_positive_cpp_llround(self):
        d = event(); d['B0TrackerHits.eDep'][0] = 10.5e-6
        d['B0TrackerRawHits.charge'][0] = 11
        self.assertEqual(match(d)[0][0]['purity'], 1.)

    def test_direct_and_secondary_shares(self):
        choices, tracks = match(event(sim=((0, 0), (0, SECONDARY_BIT)), raw=((0, 1),)))
        self.assertEqual(choices[0]['purity'], .5)
        self.assertEqual(tracks[0]['fractions'], {0: .5})

    def test_secondary_ancestry_never_matches(self):
        d = event(sim=((0, SECONDARY_BIT),))
        truth = EventTruth(d.__getitem__, IDS, 1, 10)
        self.assertEqual(truth.direct_stations[0], set())
        self.assertEqual(truth.ancestry_stations[0], {1})
        self.assertEqual(truth.match(P)[0], {})

    def test_duplicate_sim_links_do_not_change_share(self):
        d = event(sim=((0, 0), (0, SECONDARY_BIT)), raw=((0, 0, 0, 1),))
        self.assertEqual(match(d)[0][0]['purity'], .5)

    def test_unlinked_measurement_stays_in_denominator(self):
        d = event(raw=((0,), ()), measurements=((0,), (1,)), tracks=((0, 1),))
        self.assertEqual(match(d)[0][0]['purity'], .5)
        self.assertEqual(match(event(raw=((),)))[0], {})

    def test_empty_measurement_stays_in_denominator(self):
        d = event(measurements=((0,), ()), tracks=((0, 1),))
        self.assertEqual(match(d)[0][0]['purity'], .5)

    def test_unknown_share_in_multihit_measurement(self):
        d = event(raw=((0,), ()), measurements=((0, 1),))
        self.assertEqual(match(d)[0][0]['purity'], .5)

    def test_nonproton_contribution_stays_in_denominator(self):
        d = event(sim=((0, 0), (1, 0)), raw=((0, 1),), pdg=(2212, 11))
        self.assertEqual(match(d)[0][0]['purity'], .5)

    def test_primary_tie_is_deterministic(self):
        d = event(sim=((1, 0), (0, 0)), raw=((0, 1),))
        self.assertEqual(list(match(d)[0]), [0])

    def test_exact_half_many_contributors(self):
        d = event(sim=tuple([(0, 0)] * 7 + [(0, SECONDARY_BIT)] * 7), raw=(tuple(range(14)),))
        self.assertEqual(match(d)[0][0]['purity'], .5)

    def test_best_primary_track_prioritizes_purity(self):
        d = event(sim=((0, 0), (0, SECONDARY_BIT)), raw=((0,), (1,)),
                  measurements=((0,), (0,), (1,)), tracks=((0, 1, 2), (0, 1)))
        self.assertEqual(match(d)[0][0]['track'], 1)

    def test_track_ties_prefer_measurements_then_chi2_then_index(self):
        d = event(measurements=((0,), (0,)), tracks=((0,), (0, 1), (0, 1)))
        d[P + 'Tracks.chi2'] = [0., 3., 2.]
        self.assertEqual(match(d)[0][0]['track'], 2)
        d[P + 'Tracks.chi2'][1] = 2.
        self.assertEqual(match(d)[0][0]['track'], 1)

    def test_duplicate_fitted_measurement_fails(self):
        with self.assertRaisesRegex(ValueError, 'repeats'):
            match(event(tracks=((0, 0),)))

    def test_invalid_relation_and_crossed_collection_fail(self):
        for field, bad in [('_B0TrackerHits_particle.index', -1),
                           ('_B0TrackerHits_particle.collectionID', 999),
                           ('_B0TrackerRawHitAssociations_simHit.index', 100),
                           ('_B0TrackerRawHitAssociations_simHit.collectionID', IDS['MCParticles']),
                           ('_B0TrackerRecHits_rawHit.collectionID', 999),
                           ('_B0TrackerMeasurements_hits.index', 100),
                           ('_' + P + 'Tracks_measurements.collectionID', 999)]:
            with self.subTest(field=field):
                d = event(); d[field][0] = bad
                with self.assertRaises(ValueError):
                    match(d)

    def test_mismatched_cells_and_slices_fail(self):
        d = event(); d['B0TrackerHits.cellID'][0] += 1
        with self.assertRaisesRegex(ValueError, 'cell'):
            match(d)
        d = event(); d[P + 'Tracks.measurements_end'][0] = 2
        with self.assertRaisesRegex(ValueError, 'slice'):
            match(d)

    def test_physical_stations_not_cad_front_back_layers(self):
        d = event(sim=((0, 0), (0, 0), (0, SECONDARY_BIT)))
        d['B0TrackerHits.cellID'] = [150 | (l << 8) for l in (1, 2, 3)]
        direct, ancestry, owners = primary_station_sets(d.__getitem__, IDS, 2)
        self.assertEqual(direct[0], {1})
        self.assertEqual(ancestry[0], {1, 2})
        self.assertEqual(owners, [0, 0, None])




class NumericalTests(unittest.TestCase):
    def test_candidate_pairing_uses_seed_and_hits_not_winner_or_track_index(self):
        from analyze import compare_candidates
        packed = [float(i == j) for j in range(6) for i in range(j + 1)]
        def row(track, seed, hits, shift=0., covariance=None):
            return dict(key=[1, 2, track], seed=[7, seed], measurements=hits,
                        parameters=[shift, 0., 0., 0., .1, 0.],
                        bound_covariance=packed if covariance is None else covariance)
        left = [row(0, 3, [2, 1]), row(1, 4, [1, 2]), row(2, 5, [3, 4])]
        right = [row(9, 3, [1, 2], 2., [4 * x for x in packed]),
                 row(8, 4, [1, 2]), row(7, 4, [1, 2]), row(6, 6, [3, 4])]
        result = compare_candidates(left, right)
        self.assertEqual(result['unique_pairs'], 1)
        self.assertEqual(result['numerical_pairs'], 1)
        self.assertEqual(len(result['ambiguous_groups']), 1)
        self.assertEqual(len(result['gain_groups']), 1)
        self.assertEqual(len(result['loss_groups']), 1)
        self.assertEqual(result['abs_parameter_shift_in_baseline_sigma']['loc0']['median'], 2.)
        self.assertEqual(result['sigma_ratio_variant_over_baseline']['qOverP']['median'], 2.)
        right[0]['parameters'][0] = None
        invalid = compare_candidates(left, right)
        self.assertEqual(invalid['invalid_pairs'], 1)
        self.assertIsNone(invalid['abs_parameter_shift_in_baseline_sigma'])

    def test_filtered_and_unfiltered_collection_aliases(self):
        from analyze import collection_alias
        for prefix in ('B0TrackerCKF', 'B0TrackerCKFTruthSeeded'):
            for kind in ('Tracks', 'Trajectories', 'TrackParameters'):
                name = prefix + kind
                for field in ('', '.chi2', '_measurements.index'):
                    self.assertEqual(collection_alias(name + field, prefix, ''), name + field)
                    self.assertEqual(collection_alias(name + field, prefix, 'Unfiltered'),
                                     name + 'Unfiltered' + field)
                self.assertEqual(collection_alias('_' + name + '_trajectory.index', prefix, 'Unfiltered'),
                                 '_' + name + 'Unfiltered_trajectory.index')
            for untouched in ('B0TrackerMeasurements.loc.a', '_B0TrackerMeasurements_hits.index',
                              prefix + 'TrackAssociations.weight'):
                self.assertEqual(collection_alias(untouched, prefix, 'Unfiltered'), untouched)

    def test_small_positive_eigenvalue_is_accepted(self):
        import numpy as np
        from analyze import covariance_metrics
        cov = np.eye(6)
        cov[0, 1] = cov[1, 0] = np.nextafter(np.float32(1), np.float32(0))
        packed = [cov[i, j] for j in range(6) for i in range(j + 1)]
        result = covariance_metrics(packed, [0.] * 6)
        self.assertGreater(result['min_normalized_eigenvalue'], 0)
        self.assertLess(result['min_normalized_eigenvalue'], 1e-6)
        self.assertEqual(result['covariance_status'], 'positive_definite')

    def test_float32_nonpositive_matrix_is_roundoff_limited(self):
        import numpy as np
        from analyze import covariance_metrics
        cov = np.eye(6)
        # An SPD block with correlation 1-1e-8 loses its positive eigenvalue
        # when serialized to float32. Do not clamp it or report valid coverage.
        cov[0, 1] = cov[1, 0] = 1 - 1e-8
        cov = cov.astype(np.float32)
        packed = [cov[i, j] for j in range(6) for i in range(j + 1)]
        result = covariance_metrics(packed, [0.] * 6)
        self.assertLessEqual(result['min_normalized_eigenvalue'], 0)
        self.assertEqual(result['covariance_status'], 'roundoff_limited')
        self.assertGreater(result['rounding_eigenvalue_bound'], 0)
        self.assertIsNone(result['joint95'])

    def test_displaced_perigee_and_signed_charge(self):
        import numpy as np
        from analyze import perigee
        for charge in (-1, 1):
            state = perigee([2., 3., 100.], [1., 0., 2.], charge)
            np.testing.assert_allclose(state, [3., 96., 0., np.arctan2(1., 2.), charge / np.sqrt(5.)])
            # Moving the truth point along its ray leaves the perigee invariant.
            np.testing.assert_allclose(state, perigee([12., 3., 120.], [1., 0., 2.], charge))
        with self.assertRaises(ValueError):
            perigee([0, 0, 0], [0, 0, 1], 1)

    def test_packed_covariance_joint_pulls_and_angle_wrap(self):
        import numpy as np
        from analyze import covariance_metrics
        cov = np.diag([4., 9., .01, .04, .0001, 1.])
        cov[0, 1] = cov[1, 0] = 1.5
        packed = [cov[i, j] for j in range(6) for i in range(j + 1)]
        truth = np.array([0., 0., np.pi - .01, .2, .1])
        fitted = [2., 3., -np.pi + .01, .4, .11, 0.]
        result = covariance_metrics(packed, fitted, truth)
        residual = np.array([2., 3., .02, .2, .01])
        self.assertEqual(result['covariance_status'], 'positive_definite')
        np.testing.assert_allclose(result['pulls'], [1., 1., .2, 1., 1.])
        self.assertAlmostEqual(result['mahalanobis2'], residual @ np.linalg.solve(cov[:5, :5], residual))
        self.assertAlmostEqual(result['qop_residual_percent'], 10.)
        self.assertTrue(result['joint95'])
        packed[0] = -1
        self.assertEqual(covariance_metrics(packed, fitted)['covariance_status'], 'nonpositive')
        packed[0] = float('nan')
        self.assertEqual(covariance_metrics(packed, fitted)['covariance_status'], 'nonfinite')




class MomentumQualityTests(unittest.TestCase):
    def test_charge_tails_and_invalid_are_not_hidden(self):
        from analyze import momentum_quality
        rows = [dict(stations=3, matched=True, truth=[0, 0, 0, 0, .025],
                     track=dict(parameters=[0, 0, 0, 0, qop]))
                for qop in (.025, -.00125, 0, None)]
        rows.append(dict(stations=2, matched=True, truth=[0, 0, 0, 0, .025],
                         track=dict(parameters=[0, 0, 0, 0, -.025])))
        result = momentum_quality(rows)
        self.assertEqual(result['matched_protons'], 4)
        self.assertEqual(result['wrong_charge'], 1)
        self.assertEqual(result['correct_charge'], 1)
        self.assertEqual(result['invalid_or_zero_qop'], 2)
        self.assertEqual(result['momentum_abs_residual_gt100_percent'], 1)
        self.assertAlmostEqual(result['momentum_max_abs_residual_percent'], 1900)
        self.assertAlmostEqual(result['momentum_residual_percent']['median'], 950)

    def test_empty_acceptance(self):
        from analyze import momentum_quality
        self.assertIsNone(momentum_quality([])['momentum_residual_percent'])


if __name__ == '__main__':
    unittest.main()
