"""Run: python -m unittest discover -s src/tests/scripts_test -v"""

import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import awkward as ak
import numpy as np

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "compare_tracking_outputs.py"
spec = importlib.util.spec_from_file_location("compare_tracking_outputs", SCRIPT)
comparator = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = comparator
spec.loader.exec_module(comparator)


class TrackingComparisonTest(unittest.TestCase):
    def snapshot(self, values, identities=None, name="B0Params.qOverP"):
        return comparator.Snapshot(identities or [(1, 10), (1, 11)],
                                   {"B0Params": {"collectionID": 3}},
                                   {name: ak.Array(values)})

    def test_event_reordering_is_allowed(self):
        a = self.snapshot([[1.0, 2.0], [3.0]])
        b = self.snapshot([[3.0], [1.0, 2.0]], [(1, 11), (1, 10)])
        self.assertEqual(comparator.compare(a, b), [])

    def test_small_float_and_signed_zero_differences_are_detected(self):
        for x, y in [(0.125, np.nextafter(0.125, 1.0)), (0.0, -0.0)]:
            with self.subTest(x=x, y=y):
                self.assertEqual(len(comparator.compare(self.snapshot([[x], []]),
                                                        self.snapshot([[y], []]))), 1)

    def test_covariance_relation_and_object_order_are_checked(self):
        for field in ["B0Params.covariance", "_B0Tracks_measurements.index",
                      "_B0Tracks_measurements.collectionID", "B0Params.qOverP"]:
            with self.subTest(field=field):
                a = self.snapshot([[1, 2], []], name=field)
                b = self.snapshot([[2, 1], []], name=field)
                self.assertEqual(comparator.compare(a, b), [(field, 1, (1, 10))])

    def test_missing_events_and_schema_fail(self):
        a = self.snapshot([[1], []])
        with self.assertRaisesRegex(ValueError, "identity sets"):
            comparator.compare(a, self.snapshot([[1], []], [(1, 10), (1, 12)]))
        with self.assertRaisesRegex(ValueError, "schemas"):
            comparator.compare(a, self.snapshot([[1], []], name="B0Params.theta"))
        b = self.snapshot([[1], []])
        b.metadata["B0Params"]["collectionID"] = 4
        with self.assertRaisesRegex(ValueError, "metadata"):
            comparator.compare(a, b)

    def test_identity_must_be_single_and_unique(self):
        for run, event in [([[1], [1]], [[10], [10]]), ([[1]], [[]]),
                           ([[1, 2]], [[10, 11]])]:
            with self.subTest(run=run, event=event), self.assertRaises(ValueError):
                comparator.event_identities(ak.Array(run), ak.Array(event))

    def test_events_after_2000_are_checked(self):
        keys = [(1, i) for i in range(5000)]
        a = self.snapshot([[1]] * 5000, keys)
        b = self.snapshot([[1]] * 4999 + [[2]], keys)
        self.assertEqual(comparator.compare(a, b), [("B0Params.qOverP", 1, (1, 4999))])


class ReaderTest(unittest.TestCase):
    class Branch:
        def __init__(self, values):
            self.values = ak.Array(values)

        def keys(self):
            return []

        def array(self):
            return self.values

        def arrays(self):
            return self.values

    class Root(dict):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    def root(self, target=2):
        group = "events___CollectionTypeInfo"
        metadata = self.Branch([{
            group + ".name": ["EventHeader", "B0Tracks", "Hits"],
            group + ".collectionID": [0, 1, 2],
            group + ".dataType": ["Header", "Track", "Hit"],
        }])
        events = {
            "EventHeader.runNumber": self.Branch([[1]]),
            "EventHeader.eventNumber": self.Branch([[10]]),
            "B0Tracks.chi2": self.Branch([[1.0]]),
            "_B0Tracks_hits.index": self.Branch([[0]]),
            "_B0Tracks_hits.collectionID": self.Branch([[target]]),
            "Hits.position.x": self.Branch([[2.0]]),
        }
        return self.Root(events=events, podio_metadata={group: metadata})

    def test_relation_closure_is_compared(self):
        with patch.object(comparator.uproot, "open", return_value=self.root()):
            a = comparator.read_snapshot("unused")
        root = self.root()
        root["events"]["Hits.position.x"] = self.Branch([[2.1]])
        with patch.object(comparator.uproot, "open", return_value=root):
            b = comparator.read_snapshot("unused")
        self.assertEqual(set(a.metadata), {"EventHeader", "B0Tracks", "Hits"})
        self.assertEqual(comparator.compare(a, b), [("Hits.position.x", 1, (1, 10))])

    def test_dangling_relations_are_rejected(self):
        with patch.object(comparator.uproot, "open", return_value=self.root(target=99)):
            with self.assertRaisesRegex(ValueError, "Dangling relation"):
                comparator.read_snapshot("unused")

    def test_subset_collection_follows_owning_collection(self):
        root = self.root()
        events = root["events"]
        del events["B0Tracks.chi2"]
        events["B0Tracks_objIdx.index"] = events.pop("_B0Tracks_hits.index")
        events["B0Tracks_objIdx.collectionID"] = events.pop("_B0Tracks_hits.collectionID")
        with patch.object(comparator.uproot, "open", return_value=root):
            snapshot = comparator.read_snapshot("unused")
        self.assertEqual(set(snapshot.metadata), {"EventHeader", "B0Tracks", "Hits"})
        self.assertIn("B0Tracks_objIdx.index", snapshot.fields)

    def test_missing_header_is_rejected(self):
        root = self.root()
        del root["events"]["EventHeader.eventNumber"]
        with patch.object(comparator.uproot, "open", return_value=root):
            with self.assertRaisesRegex(ValueError, "Missing EventHeader"):
                comparator.read_snapshot("unused")


if __name__ == "__main__":
    unittest.main()
