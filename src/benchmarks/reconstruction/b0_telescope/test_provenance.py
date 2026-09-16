# SPDX-License-Identifier: LGPL-3.0-or-later
import copy
import tempfile
import unittest
from pathlib import Path
import numpy as np
from compare_priors import compare
from provenance import describe, verify_file, verify_outputs, xml_dependencies
from run_b0_telescope import command

class ProvenanceTests(unittest.TestCase):
    def test_hash_and_tamper(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'input'; p.write_text('original')
            record = describe(p); verify_file(record)
            manifest = {'schema_version': 1, 'status': 'completed', 'outputs': {'reco': record, 'truth': record}}
            verify_outputs(manifest, p, p)
            p.write_text('changed')
            with self.assertRaises(ValueError): verify_file(record)
            with self.assertRaises(ValueError): verify_outputs(manifest, p, p)

    def test_recursive_xml_and_cycles(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / 'a.xml', Path(tmp) / 'b.xml'
            a.write_text('<lccdd><include ref="b.xml"/></lccdd>')
            b.write_text('<lccdd><include ref="a.xml"/></lccdd>')
            self.assertEqual(len(xml_dependencies(a)), 2)
            b.write_text('<lccdd><include ref="missing.xml"/></lccdd>')
            with self.assertRaises(FileNotFoundError): xml_dependencies(a)

    def test_unresolved_or_network_include(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'x.xml'
            for value in ('https://example.test/a.xml', '${B0_UNDEFINED_VARIABLE}/a.xml'):
                p.write_text(f'<lccdd><include ref="{value}"/></lccdd>')
                with self.assertRaises(ValueError): xml_dependencies(p)

    def test_command_and_reserved_options(self):
        p = Path('/a path with spaces')
        c = command('eicrecon', p, p, p, p, p, ('B0TelescopeCKF',), 10, ['B0TelescopeSeeds:useTiming=true'])
        self.assertEqual(c[-1], str(p))
        self.assertIn('-Pb0_telescope:B0TelescopeCKF:weakPriorScale=10', c)
        with self.assertRaises(ValueError): command('x', p, p, p, p, p, ('B0TelescopeCKF',), 1, ['B0TelescopeSeeds:diagnosticsFile=/elsewhere'])
        with self.assertRaises(ValueError): command('x', p, p, p, p, p, ('unknown',), 1, [])

    def test_prior_comparison(self):
        row = {'chain': 'a', 'event': 1, 'particle_id': [3, 4], 'measurement_ids': [[9, 1], [9, 2]],
            'surface': '123', 'state': [0, 0, 0, .02, .1, 0], 'covariance': np.eye(6).tolist()}
        same = compare([row], [row]); self.assertEqual(same['max_mean_shift'], 0)
        changed = copy.deepcopy(row); changed['state'][2] += 2 * np.pi
        self.assertLess(compare([row], [changed])['max_mean_shift'], 1e-12)
        changed['measurement_ids'][0][1] = 3
        result = compare([row], [changed]); self.assertEqual(result['changed_hits'], 1)
        self.assertIsNone(result['max_mean_shift'])
        with self.assertRaises(ValueError): compare([row, row], [row])
        changed = copy.deepcopy(row); changed['covariance'][0][0] = -1
        with self.assertRaises(np.linalg.LinAlgError): compare([row], [changed])

if __name__ == '__main__': unittest.main()
