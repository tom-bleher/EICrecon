"""Paired comparison on the same truth-associated candidates, with no truth ranking."""

import argparse
import json
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input", help="Numeric estimator input (reads its .json sidecar)")
parser.add_argument("result", help="C++ estimator output")
args = parser.parse_args()
with open(args.input + ".json") as stream:
    rows = json.load(stream)
if not rows:
    raise ValueError("No associated seeds to compare")
x = np.loadtxt(args.result, ndmin=2)
if x.shape != (len(rows), 12) or not np.array_equal(x[:, 0], np.arange(len(rows))):
    raise ValueError("Result IDs or row/column counts do not match extracted seeds")
qfirst = np.array([r["q_truth_first"] for r in rows])
qorig = np.array([r["q_truth_origin"] for r in rows])
qc = np.array([r["q_current"] for r in rows])
dt = np.array([r["direction_truth_first"] for r in rows])


def stats(v):
    good = np.isfinite(v)
    valid = v[good]
    if not len(valid):
        return dict(n=len(v), nonfinite=len(v))
    p = np.quantile(valid, [0.16, 0.5, 0.84])
    return dict(
        n=len(v),
        nonfinite=int((~good).sum()),
        median=float(p[1]),
        half68=float((p[2] - p[0]) / 2),
        rms=float(np.sqrt(np.mean(valid**2))),
        p95_abs=float(np.quantile(abs(valid), 0.95)),
    )


def first(clean=False):
    seen = set()
    mask = []
    for r in rows:
        key = (r["event"], r["mc"])
        take = (not clean or r["clean_forward"]) and key not in seen
        mask.append(take)
        if take:
            seen.add(key)
    return np.array(mask)


clean = first(True)
result = {
    "selection": "First emitted truth-associated seed per primary proton; no closest truth q/p or charge selection. Clean-forward requires every associated sim hit to retain >50% origin momentum and lab pz>0. All q/p comparisons below use truth at the first selected sensor. Direction comparison also at that sensor; baseline propagated with full field RK4, no material."
}
for name, mask in [
    ("all_associated", np.ones(len(x), bool)),
    ("first_per_truth", first()),
    ("first_clean_per_truth", clean),
    ("first_clean_three_hits", clean & np.array([r["nhits"] == 3 for r in rows])),
    ("first_clean_four_hits", clean & np.array([r["nhits"] == 4 for r in rows])),
]:
    if not mask.any():
        continue
    d = {}
    for alg, q in [("current", qc), ("acts", x[:, 1])]:
        d[alg + "_relative_qop_percent"] = stats((100 * (q - qfirst) / qfirst)[mask])
        d[alg + "_wrong_charge"] = int(np.sum(q[mask] < 0))
    for alg, direction in [("acts", x[:, 2:5]), ("current_transported", x[:, 5:8])]:
        dot = np.sum(direction * dt, axis=1)
        angle = np.arctan2(np.linalg.norm(np.cross(direction, dt), axis=1), dot) * 1000
        d[alg + "_angle_mrad"] = stats(angle[mask])
        for i, c in enumerate(["x", "y"]):
            d[alg + "_" + c + "_slope_residual_mrad"] = stats(
                (1000 * (direction[:, i] / direction[:, 2] - dt[:, i] / dt[:, 2]))[mask]
            )
    d["current_transport_position_distance_mm"] = stats(x[mask, 11])
    result[name] = d
result["field_T_ranges"] = [
    np.min(x[:, 8:11], axis=0).tolist(),
    np.max(x[:, 8:11], axis=0).tolist(),
]
# Paired bootstrap, so both estimates use the same resampled particles.
rng = np.random.default_rng(3821)
indices = np.flatnonzero(clean)
changes = []
for _ in range(1000 if len(indices) else 0):
    pick = rng.choice(indices, len(indices))
    w = []
    for q in [qc, x[:, 1]]:
        lo, hi = np.quantile(
            100 * (q[pick] - qfirst[pick]) / qfirst[pick], [0.16, 0.84]
        )
        w.append((hi - lo) / 2)
    changes.append(w[1] - w[0])
result[
    "bootstrap_95percent_interval_acts_minus_current_qop_width_percentage_points"
] = (np.quantile(changes, [0.025, 0.975]).tolist() if changes else None)
print(json.dumps(result, indent=2))
