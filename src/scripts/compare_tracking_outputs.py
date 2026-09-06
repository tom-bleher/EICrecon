#!/usr/bin/env python3
"""Compare exact persisted tracking data, allowing only event-entry reordering.

Defaults to B0* collections and their recursive PODIO relation closure. Use
--collections to select another explicit comma-separated collection list.
Requires one EventHeader per entry with unique (runNumber, eventNumber).
Older validation files without EventHeader cannot establish this comparison.
All events, primitive fields, object order and relation indices/IDs are checked;
float storage bits (including signed zero and NaNs) must agree. ROOT container
bytes and unrelated collections are outside the comparison scope.

Run with Python, uproot, awkward and numpy in eic-shell.
"""

import argparse
from dataclasses import dataclass
import sys

import awkward as ak
import numpy as np
import uproot


@dataclass
class Snapshot:
    identities: list
    metadata: dict
    fields: dict


def event_identities(run, event):
    if len(run) != len(event):
        raise ValueError("EventHeader run/event entry counts differ")
    identities = []
    for r, e in zip(ak.to_list(run), ak.to_list(event)):
        if len(r) != 1 or len(e) != 1:
            raise ValueError("Expected exactly one EventHeader per event")
        identities.append((int(r[0]), int(e[0])))
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate EventHeader (runNumber, eventNumber)")
    return identities


def collection_fields(tree, name):
    """Return canonical leaf names, retaining vector and relation branches."""
    fields = {}
    for path in tree.keys():
        key = path.rsplit("/", 1)[-1]
        if not (key.startswith(name + ".") or key.startswith(name + "_objIdx.")
                or key.startswith("_" + name + "_")
                or key.startswith("_" + name + ".")):
            continue
        branch = tree[path]
        if branch.keys():
            continue
        if key in fields:
            raise ValueError(f"Ambiguous branch {key}")
        fields[key] = branch.array()
    if not fields:
        raise ValueError(f"No persisted fields for collection {name}")
    return fields


def read_snapshot(path, collections=None):
    with uproot.open(path) as root:
        tree = root["events"]
        try:
            identities = event_identities(tree["EventHeader.runNumber"].array(),
                                          tree["EventHeader.eventNumber"].array())
        except KeyError as exc:
            raise ValueError("Missing EventHeader identity; persist EventHeader and rerun") from exc
        group = "events___CollectionTypeInfo"
        records = ak.to_list(root["podio_metadata"][group].arrays())
        if len(records) != 1:
            raise ValueError("Expected one collection metadata entry")
        record = records[0]
        # TTree has dotted fields; RNTuple has a nested record.
        record = record.get(group, record)
        columns = {k.rsplit(".", 1)[-1]: v for k, v in record.items()}
        names = columns["name"]
        metadata = {name: {k: v[i] for k, v in columns.items()} for i, name in enumerate(names)}
        by_id = {v["collectionID"]: k for k, v in metadata.items()}
        if len(by_id) != len(names) or len(metadata) != len(names):
            raise ValueError("Duplicate collection name or ID in metadata")
        requested = set(collections) if collections else {n for n in names if n.startswith("B0")}
        if not requested:
            raise ValueError("No collections selected")
        pending = requested | {"EventHeader"}
        fields = {}
        selected = {}
        while pending:
            name = pending.pop()
            if name in selected:
                continue
            if name not in metadata:
                raise ValueError(f"Collection {name} is not persisted")
            selected[name] = metadata[name]
            current = collection_fields(tree, name)
            fields.update(current)
            for key, values in current.items():
                if not key.endswith(".collectionID"):
                    continue
                index_key = key.removesuffix(".collectionID") + ".index"
                if index_key not in current:
                    raise ValueError(f"Missing relation indices for {key}")
                ids = np.asarray(ak.flatten(values, axis=None))
                indices = np.asarray(ak.flatten(current[index_key], axis=None))
                if len(ids) != len(indices):
                    raise ValueError(f"Relation ID/index size mismatch for {key}")
                for cid in np.unique(ids[indices >= 0]):
                    if int(cid) not in by_id:
                        raise ValueError(f"Dangling relation {key}: collectionID {cid} is not persisted")
                    pending.add(by_id[int(cid)])
        return Snapshot(identities, selected, fields)


def exact_value(value):
    """Preserve numeric dtype, shape and bits without float conversion."""
    array = np.asarray(value)
    if array.dtype.kind == "O":
        raise ValueError("Unsupported non-primitive event field")
    return array.dtype.str, array.shape, array.tobytes()


def compare(left, right):
    if set(left.identities) != set(right.identities):
        raise ValueError("Event identity sets differ "
                         f"({len(set(left.identities) - set(right.identities))} left-only, "
                         f"{len(set(right.identities) - set(left.identities))} right-only)")
    if left.metadata != right.metadata:
        raise ValueError("Selected collection metadata/IDs differ")
    if left.fields.keys() != right.fields.keys():
        raise ValueError("Selected branch schemas differ")
    order = {key: i for i, key in enumerate(right.identities)}
    differences = []
    for name in sorted(left.fields):
        a, b = left.fields[name], right.fields[name]
        if len(a) != len(left.identities) or len(b) != len(right.identities):
            raise ValueError(f"Event count mismatch in {name}")
        count = 0
        first = None
        for i, key in enumerate(left.identities):
            if exact_value(a[i]) != exact_value(b[order[key]]):
                count += 1
                if first is None:
                    first = key
        if count:
            differences.append((name, count, first))
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left")
    parser.add_argument("right")
    parser.add_argument("--collections", help="Comma-separated collection names (default: all B0*)")
    args = parser.parse_args()
    collections = args.collections.split(",") if args.collections else None
    try:
        left = read_snapshot(args.left, collections)
        right = read_snapshot(args.right, collections)
        differences = compare(left, right)
    except (ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print("Collections (including relation closure): " + ", ".join(sorted(left.metadata)))
    print(f"Compared all {len(left.identities)} events, {len(left.fields)} fields; "
          "exact storage bits and within-event order, matched by EventHeader")
    for name, count, first in differences:
        print(f"DIFF {name}: {count} events; first (run, event)={first}")
    print(f"{'FAIL' if differences else 'PASS'}: {len(differences)} differing fields")
    return 1 if differences else 0


if __name__ == "__main__":
    sys.exit(main())
