# SPDX-License-Identifier: LGPL-3.0-or-later
import json
import math
from pathlib import Path
import tempfile
import unittest
import numpy as np
from analyze_b0_telescope import (Selection, TruthReference, evaluate, match_particle,
                                  rate, residual_diagnostics, strict_json, write_report)


def fixtures():
    truth = TruthReference(events={7}, configuration={"schema_version": 1})
    for identity, pdg in (((3, 0), 2212), ((3, 1), 211), ((4, 0), 2212)):
        key = (7, *identity)
        truth.particles[key] = {"supported": True, "pdg": pdg, "n_stations": 4,
            "p_production_gev": 20, "vertex_mm": [0, 0, 1000], "particle_id": identity}
        truth.states[(key, "18446744073709551600")] = {
            "surface": "18446744073709551600", "state": [0, 0, .2, .03, .05, 20]}
    def track(index, identity, weight=1):
        return {"id": (12, index), "matches": [] if identity is None else [{"particle_id": identity, "weight": weight}],
            "n_measurements": 8, "chi2": 4, "ndf": 11, "measurement_ids": [(10, i) for i in range(8)],
            "parameters": {"surface": "18446744073709551600", "state": [0, 0, .2, .03, .05, 20],
                           "covariance": np.eye(6).tolist()}}
    event = {"event": 7, "n_hits": 16, "chains": {"test": {
        "seeds": [{"matches": [{"particle_id": (3, 0), "weight": 1}]}],
        "tracks": [track(0, (3, 0)), track(1, (3, 0)), track(2, (3, 1)), track(3, None)]}}}
    return truth, event


class AnalysisTests(unittest.TestCase):
    def test_real_non_target_is_not_fake(self):
        truth, event = fixtures()
        report, rows = evaluate([event], truth, Selection(frozenset({2212})), ("test",))
        result = report["chains"]["test"]
        self.assertEqual(report["target_particles"], 2)
        self.assertEqual(result["total_efficiency"]["value"], .5)
        self.assertEqual(result["fake_fraction"]["value"], .25)
        self.assertEqual(result["counts"]["matched_outside_target"], 1)
        self.assertEqual(result["counts"]["duplicate_tracks"], 1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["surface"], "18446744073709551600")

    def test_collection_id_is_part_of_truth_identity(self):
        truth, event = fixtures()
        event["chains"]["test"]["tracks"][1]["matches"][0]["particle_id"] = (4, 0)
        report, _ = evaluate([event], truth, Selection(frozenset({2212})), ("test",))
        self.assertEqual(report["chains"]["test"]["total_efficiency"]["value"], 1)
        self.assertEqual(report["chains"]["test"]["recovered_without_matching_seed"], 1)

    def test_empty_denominator_is_undefined(self):
        self.assertIsNone(rate(0, 0)["value"])
        self.assertEqual(rate(1, 1)["value"], 1)
        self.assertGreater(rate(0, 10)["high95"], 0)
        with self.assertRaises(ValueError): rate(2, 1)

    def test_strict_majority_ties_and_order(self):
        rows = [{"particle_id": [3, 1], "weight": .5}, {"particle_id": [3, 2], "weight": .5}]
        self.assertEqual(match_particle(rows, .5), (None, "ambiguous"))
        self.assertEqual(match_particle(rows[::-1], .5), (None, "ambiguous"))
        self.assertEqual(match_particle([], .5), (None, "unassociated"))
        with self.assertRaises(ValueError): match_particle([{"particle_id": [3, 0], "weight": float("nan")}], .5)

    def test_wrapped_phi_and_spd_covariance(self):
        estimate = {"surface": "1", "state": [0, 0, math.pi - .01, .1, -.1], "covariance": np.eye(5).tolist()}
        truth = {"surface": "1", "state": [0, 0, -math.pi + .01, .1, .1]}
        result = residual_diagnostics(estimate, truth)
        self.assertAlmostEqual(result["residual"][2], -.02)
        self.assertTrue(result["wrong_charge"])
        self.assertAlmostEqual(result["d2"], .0004 + .04)
        estimate["covariance"][0][1] = estimate["covariance"][1][0] = 2
        with self.assertRaises(ValueError): residual_diagnostics(estimate, truth)

    def test_reference_surface_must_match(self):
        estimate = {"surface": "1", "state": [0]*5, "covariance": np.eye(5).tolist()}
        with self.assertRaises(ValueError): residual_diagnostics(estimate, {"surface": "2", "state": [0]*5})
        estimate["covariance"] = [1, 2, 3]
        with self.assertRaises(ValueError): residual_diagnostics(estimate, {"surface": "1", "state": [0]*5})

    def test_missing_reference_does_not_change_efficiency(self):
        truth, event = fixtures(); truth.states.clear()
        report, rows = evaluate([event], truth, Selection(frozenset({2212})), ("test",))
        self.assertEqual(report["chains"]["test"]["total_efficiency"]["value"], .5)
        self.assertEqual(report["chains"]["test"]["counts"]["missing_same_surface_truth"], 1)
        self.assertEqual(rows, [])

    def test_duplicate_events_fail(self):
        truth, event = fixtures()
        with self.assertRaises(ValueError): evaluate([event, event], truth, Selection(), ("test",))

    def test_missing_chain_fails(self):
        truth, event = fixtures()
        with self.assertRaises(ValueError): evaluate([event], truth, Selection(), ("missing",))

    def test_truth_file_validation_and_report_serialization(self):
        truth, event = fixtures()
        report, rows = evaluate([event], truth, Selection(frozenset({2212})), ("test",))
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory); write_report(report, rows, out)
            self.assertEqual(strict_json((out / "summary.json").read_text())["events"], 1)
            path = out / "truth.jsonl"
            path.write_text(json.dumps({"type": "configuration", "schema_version": 1}) + '\n' +
                            json.dumps({"type": "event", "event": 7, "seeds": 0}) + '\n')
            self.assertEqual(TruthReference.read(path).events, {7})
            with path.open("a") as stream: stream.write(json.dumps({"type": "event", "event": 7}) + '\n')
            with self.assertRaises(ValueError): TruthReference.read(path)

    def test_invalid_configuration(self):
        for threshold in (.4, 1, float("nan")):
            with self.assertRaises(ValueError): Selection(match_threshold=threshold)
        with self.assertRaises(ValueError): strict_json('{"x":NaN}')
        with self.assertRaises(ValueError): Selection(min_stations=2)

    def test_zero_curvature_does_not_invent_zero_momentum(self):
        estimate = {"surface": "1", "state": [0, 0, 0, .1, 0], "covariance": np.eye(5).tolist()}
        result = residual_diagnostics(estimate, {"surface": "1", "state": [0, 0, 0, .1, .05]})
        self.assertTrue(result["unresolved_charge"])
        self.assertIsNone(result["relative_momentum_residual"])
        self.assertTrue(math.isfinite(result["d2"]))

if __name__ == "__main__": unittest.main()
