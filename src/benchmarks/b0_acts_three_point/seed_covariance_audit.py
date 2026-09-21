"""Truth-associated seed audit; no nearest-q/p matching or residual clipping."""

import argparse
import json
import numpy as np
from scipy.stats import chi2

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input", help="JSON sidecar produced by extract.py")
args = parser.parse_args()
with open(args.input) as stream:
    rows = json.load(stream)
if not rows:
    raise ValueError("No associated seeds to audit")


def first_mask(rows, clean=False):
    seen = set()
    mask = []
    for r in rows:
        key = (r["event"], r["mc"])
        take = (not clean or r["clean_forward"]) and key not in seen
        mask.append(take)
        if take:
            seen.add(key)
    return np.array(mask)


def summary(indices):
    if not len(indices):
        return {"n": 0}
    residual = []
    cov = []
    for i in indices:
        r = rows[i]
        p = np.array(r["truth_momentum"])
        p /= np.linalg.norm(p)
        v = np.array(r["truth_vertex"])
        ph = np.arctan2(p[1], p[0])
        th = np.arccos(p[2])
        truth = [
            -v[0] * np.sin(ph) + v[1] * np.cos(ph),
            v[2] - (v[0] * np.cos(ph) + v[1] * np.sin(ph)) / np.tan(th),
            ph,
            th,
            r["q_truth_origin"],
        ]
        d = np.array(r["seed_state"]) - truth
        d[2] = np.remainder(d[2] + np.pi, 2 * np.pi) - np.pi
        flat = r["covariance"]
        c = np.array(
            [
                [flat[max(a, b) * (max(a, b) + 1) // 2 + min(a, b)] for b in range(5)]
                for a in range(5)
            ]
        )
        residual.append(d)
        cov.append(c)
    d = np.array(residual)
    c = np.array(cov)
    count = len(indices)
    diagonal = np.diagonal(c, axis1=1, axis2=2)
    finite_covariance = np.isfinite(c).all(axis=(1, 2))
    positive_diagonal = (np.isfinite(diagonal) & (diagonal > 0)).all(axis=1)
    finite_residual = np.isfinite(d).all(axis=1)
    # Never pass NaN/negative/zero variances to sqrt or the eigensolver.
    marginal_valid = np.isfinite(diagonal) & (diagonal > 0) & np.isfinite(d)
    pull = np.full_like(d, np.nan, dtype=float)
    pull[marginal_valid] = d[marginal_valid] / np.sqrt(diagonal[marginal_valid])
    correlation_valid = finite_covariance & positive_diagonal
    corr = np.zeros_like(c)
    sigma = np.sqrt(diagonal[correlation_valid])
    corr[correlation_valid] = (
        c[correlation_valid] / sigma[:, :, None] / sigma[:, None, :]
    )
    minimum_eigen = np.full(count, np.nan)
    if correlation_valid.any():
        minimum_eigen[correlation_valid] = np.linalg.eigvalsh(corr[correlation_valid])[
            :, 0
        ]
    joint_valid = correlation_valid & finite_residual & (minimum_eigen > 0)
    result = {
        "n": count,
        "nonfinite_covariance": int((~finite_covariance).sum()),
        "nonpositive_or_nonfinite_diagonal": int((~positive_diagonal).sum()),
        "nonfinite_residual": int((~finite_residual).sum()),
        "correlation_eigenvalue_tested": int(correlation_valid.sum()),
        "nonpositive_correlation_eigenvalues": int(
            (minimum_eigen[correlation_valid] <= 0).sum()
        ),
        "minimum_correlation_eigenvalue": (
            float(minimum_eigen[correlation_valid].min())
            if correlation_valid.any()
            else None
        ),
        "pulls": {},
    }
    for j, name in enumerate(["loc0", "loc1", "phi", "theta", "q/p"]):
        valid = marginal_valid[:, j]
        values = pull[valid, j]
        quantiles = np.quantile(values, [0.16, 0.5, 0.84]) if len(values) else None
        result["pulls"][name] = {
            "valid": int(valid.sum()),
            "invalid": int((~valid).sum()),
            "median": float(quantiles[1]) if quantiles is not None else None,
            "half68": (
                float((quantiles[2] - quantiles[0]) / 2)
                if quantiles is not None
                else None
            ),
            "coverage_denominator": count,
            "inside_1sigma": float(np.sum(abs(values) <= 1) / count),
            "abs_pull_gt3": float(np.sum(abs(values) > 3) / count),
        }
    mahal = np.einsum(
        "ni,ni->n",
        pull[joint_valid],
        np.linalg.solve(corr[joint_valid], pull[joint_valid, :, None])[:, :, 0],
    )
    result["joint_5d"] = {
        "valid": int(joint_valid.sum()),
        "invalid": int((~joint_valid).sum()),
        "coverage_denominator": count,
        "median_chi2": float(np.median(mahal)) if len(mahal) else None,
        "coverage_at_68percent_gaussian_ellipsoid": float(
            np.sum(mahal <= chi2.ppf(0.68, 5)) / count
        ),
        "coverage_at_95percent_gaussian_ellipsoid": float(
            np.sum(mahal <= chi2.ppf(0.95, 5)) / count
        ),
    }
    return result


result = {
    "matching": "Each seed hit has exactly one raw->sim relation; all hits point to same generatorStatus1 PDG2212 particle. First emitted seed per truth particle, never closest truth momentum. Clean-forward additionally requires every associated simulated hit to retain >50% origin momentum and positive lab pz."
}
for name, mask in [
    ("all_associated", np.ones(len(rows), bool)),
    ("first_per_truth", first_mask(rows)),
    ("first_clean_per_truth", first_mask(rows, True)),
]:
    result[name] = summary(np.flatnonzero(mask))
base = first_mask(rows, True)
for n in [3, 4]:
    mask = base & np.array([r["nhits"] == n for r in rows])
    if mask.any():
        result[f"first_clean_{n}hits"] = summary(np.flatnonzero(mask))
print(json.dumps(result, indent=2, allow_nan=False))
