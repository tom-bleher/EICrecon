"""Regression: malformed covariance entries remain visible in JSON denominators."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class CovarianceAuditTest(unittest.TestCase):
    def test_invalid_covariances_are_counted(self):
        covariance = [0.0] * 21
        for index in [0, 2, 5, 9, 14, 20]:
            covariance[index] = 1.0
        base = dict(
            event=0,
            mc=0,
            clean_forward=True,
            nhits=4,
            truth_momentum=[1.0, 0.0, 10.0],
            truth_vertex=[0.0] * 3,
            q_truth_origin=0.1,
            seed_state=[0.0, 0.0, 0.0, 0.1, 0.1],
            covariance=covariance,
        )
        rows = [copy.deepcopy(base) for _ in range(6)]
        for index, row in enumerate(rows):
            row["event"] = index
        rows[1]["covariance"][0] = float("nan")
        rows[2]["covariance"][0] = -1.0
        rows[3]["covariance"][0] = 0.0
        rows[4]["covariance"][1] = float("nan")
        rows[5]["covariance"][1] = 2.0  # Positive diagonals, indefinite matrix.
        for clean in [True, False]:
            with self.subTest(clean_subset=clean):
                for row in rows:
                    row["clean_forward"] = clean
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "rows.json"
                    path.write_text(json.dumps(rows))
                    run = subprocess.run(
                        [
                            sys.executable,
                            str(Path(__file__).with_name("seed_covariance_audit.py")),
                            str(path),
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                report = json.loads(run.stdout)
                result = report["all_associated"]
                if not clean:
                    self.assertEqual(report["first_clean_per_truth"], {"n": 0})
                else:
                    self.assertEqual(report["first_clean_per_truth"]["n"], 6)
                self.assertEqual(result["nonfinite_covariance"], 2)
                self.assertEqual(result["nonpositive_or_nonfinite_diagonal"], 3)
                self.assertEqual(result["correlation_eigenvalue_tested"], 2)
                self.assertEqual(result["nonpositive_correlation_eigenvalues"], 1)
                joint = result["joint_5d"]
                self.assertEqual(
                    (joint["valid"], joint["invalid"], joint["coverage_denominator"]),
                    (1, 5, 6),
                )
                self.assertAlmostEqual(
                    joint["coverage_at_95percent_gaussian_ellipsoid"], 1 / 6
                )
                self.assertEqual(result["pulls"]["loc0"]["invalid"], 3)


if __name__ == "__main__":
    unittest.main()
