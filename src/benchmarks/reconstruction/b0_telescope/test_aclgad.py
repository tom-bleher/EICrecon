# SPDX-License-Identifier: LGPL-3.0-or-later
"""Pure adapters/provenance tests; these are not PODIO or detector smoke tests."""
from pathlib import Path
import json
import math
import tempfile
import unittest
from aclgad_analysis import measurement_votes, seed_votes
from run_aclgad_comparison import command, digest, simulation_inputs, CHAINS


class WeightTests(unittest.TestCase):
    def setUp(self):
        self.a, self.b = (2, 0), (2, 1)
        self.hits = [{'id': (4, i), 'raw': (3, i)} for i in range(3)]
        self.raw = {(3, 0): {self.a: .2}, (3, 1): {self.b: 1.0}}

    def test_noise_fraction_not_promoted(self):
        self.assertEqual(measurement_votes(self.hits[:1], [1], self.raw), {self.a: .2})

    def test_cluster_charge_weights_and_noise(self):
        result = measurement_votes(self.hits, [.5, .25, .25], self.raw)
        self.assertAlmostEqual(result[self.a], .1)
        self.assertAlmostEqual(result[self.b], .25)
        self.assertLess(sum(result.values()), 1)

    def test_legacy_weights_only_normalized_downward(self):
        self.raw[(3, 0)] = {self.a: 2, self.b: 1}
        result = measurement_votes(self.hits[:1], [], self.raw)
        self.assertAlmostEqual(result[self.a], 2 / 3)
        self.assertAlmostEqual(sum(result.values()), 1)

    def test_invalid_weights(self):
        for weights in ([1], [1, -1, 1], [1, math.nan, 1], [0, 0, 0]):
            with self.assertRaises(ValueError): measurement_votes(self.hits, weights, self.raw)
        with self.assertRaises(ValueError): measurement_votes([self.hits[0]] * 2, [], self.raw)
        for bad in (math.nan, -1):
            self.raw[(3, 0)] = {self.a: bad}
            with self.assertRaises(ValueError): measurement_votes(self.hits, [], self.raw)

    def test_seed_weights_measurements_not_channels(self):
        ms = {0: {'hits': self.hits[:2], 'weights': [.5, .5]},
              1: {'hits': self.hits[2:], 'weights': [1]}}
        reverse = {(4, 0): 0, (4, 1): 0, (4, 2): 1}
        self.raw[(3, 0)] = {self.a: 1}; self.raw[(3, 1)] = {self.a: 1}; self.raw[(3, 2)] = {self.b: 1}
        self.assertEqual(seed_votes([(4, 0), (4, 1), (4, 2)], ms, reverse, self.raw), {self.a: .5, self.b: .5})
        with self.assertRaises(ValueError): seed_votes([(4, 0)], ms, reverse, self.raw)
        with self.assertRaises(ValueError): seed_votes([(4, 4)], ms, reverse, self.raw)
        with self.assertRaises(ValueError): seed_votes([(4, 0), (4, 0)], ms, reverse, self.raw)


class CommandTests(unittest.TestCase):
    def make(self, **kwargs):
        args = dict(executable='/bin/eicrecon', sim=Path('/tmp/sim'), compact=Path('/tmp/compact'),
            geometry_calibration=Path('/tmp/calib'), material=Path('/tmp/map'), response=Path('/tmp/my response.cfg'),
            out=Path('/tmp/out'), chains=CHAINS, allow_experimental=False, prior_scale=1, settings=[])
        args.update(kwargs)
        return command(**args)

    def test_optin_collections_and_literal_paths(self):
        result = self.make()
        self.assertIn('-Pplugins=b0_telescope,b0_aclgad', result)
        self.assertIn('-Pb0_telescope:response=aclgad', result)
        self.assertIn('-Pb0_aclgad:B0ACLGADResponse:calibrationFile=/tmp/my response.cfg', result)
        cols = next(x for x in result if x.startswith('-Ppodio:output_include_collections=')).split('=', 1)[1].split(',')
        self.assertIn('B0ACLGADMeasurements', cols)
        self.assertNotIn('B0TrackerMeasurements', cols)
        self.assertIn('-Pb0_aclgad:B0ACLGADResponse:allowExperimental=false', result)

    def test_reserved_settings_and_invalid_controls(self):
        for setting in ('B0TelescopeSeeds:diagnosticsFile=x', 'B0TelescopeCKF:weakPriorScale=0',
                        'response=effective', 'B0TelescopeSeeds:unknown=', 'evil:parameter=1'):
            with self.assertRaises(ValueError): self.make(settings=[setting])
        for chains in ((), ('bad',), (CHAINS[0], CHAINS[0])):
            with self.assertRaises(ValueError): self.make(chains=chains)
        for scale in (0, -1, math.inf):
            with self.assertRaises(ValueError): self.make(prior_scale=scale)


class ReceiptTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        self.root = Path(temp.name); self.profile = self.root / 'profile'; self.profile.mkdir()
        (self.profile / 'compact.xml').write_text('<lccdd/>')
        (self.profile / 'B0.xml').write_text('B0 physical profile fixture')
        self.sim = self.root / 'sim.root'; self.sim.write_bytes(b'not a detector sample')
        self.manifest = self.profile / 'b0_readout_contract.json'
        self.manifest.write_text(json.dumps({'schema_version': 1, 'kind': 'B0ACLGADReadout', 'response_version': 1,
            'compact': 'compact.xml', 'readout_file': 'B0.xml',
            'files': {name: digest(self.profile / name) for name in ('compact.xml', 'B0.xml')}}))
        self.receipt = self.root / 'simulation.json'
        self.data = {'schema_version': 1, 'kind': 'B0ACLGADSimulation', 'status': 'completed',
            'profile': str(self.profile), 'profile_sha256': digest(self.manifest),
            'output': {'path': str(self.sim), 'sha256': digest(self.sim)}}
        self.save()

    def save(self): self.receipt.write_text(json.dumps(self.data))

    def test_matching_simulation(self):
        sim, compact, hashes = simulation_inputs(self.profile, self.receipt)
        self.assertEqual(sim, self.sim)
        self.assertEqual(compact, self.profile / 'compact.xml')
        self.assertIn(str(self.receipt), hashes)

    def test_simulation_tamper(self):
        self.sim.write_bytes(b'changed')
        with self.assertRaises(ValueError): simulation_inputs(self.profile, self.receipt)

    def test_profile_tamper(self):
        (self.profile / 'B0.xml').write_text('changed')
        with self.assertRaises(ValueError): simulation_inputs(self.profile, self.receipt)

    def test_failed_receipt_and_other_profile(self):
        self.data['status'] = 'failed'; self.save()
        with self.assertRaises(ValueError): simulation_inputs(self.profile, self.receipt)
        self.data['status'] = 'completed'; self.data['profile'] = '/other'; self.save()
        with self.assertRaises(ValueError): simulation_inputs(self.profile, self.receipt)

    def test_unrecorded_asset(self):
        (self.profile / 'extra').write_text('extra')
        with self.assertRaises(ValueError): simulation_inputs(self.profile, self.receipt)

if __name__ == '__main__': unittest.main()
