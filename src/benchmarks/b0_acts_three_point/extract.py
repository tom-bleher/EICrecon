"""Extract the same persisted candidates for both estimators; never truth-rank seeds."""

import argparse
import json
from pathlib import Path
import awkward as ak
import numpy as np
import uproot

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "input", help="PODIO ROOT output including seed hits and truth relation closure"
)
parser.add_argument(
    "output", help="Numeric estimator input, with .json and .counts.json sidecars"
)
parser.add_argument("--events", type=int, default=100)
parser.add_argument(
    "--crossing-angle", type=float, default=-0.025, help="Ion frame rotation in radians"
)
args = parser.parse_args()
if args.events <= 0 or not np.isfinite(args.crossing_angle):
    parser.error("events must be positive and angle finite")
f = uproot.open(args.input)
a = f["events"].arrays(entry_stop=args.events)
group = "events___CollectionTypeInfo"
record = ak.to_list(f["podio_metadata"][group].arrays())[0]
record = record.get(group, record)
columns = {k.rsplit(".", 1)[-1]: v for k, v in record.items()}
collection_ids = dict(zip(columns["name"], columns["collectionID"]))
relation_targets = {
    "_B0TrackerRecHits_rawHit": "B0TrackerRawHits",
    "_B0TrackerRawHitAssociations_rawHit": "B0TrackerRawHits",
    "_B0TrackerRawHitAssociations_simHit": "B0TrackerHits",
    "_B0TrackerHits_particle": "MCParticles",
    "_B0TrackerSeeds_hits": "B0TrackerRecHits",
    "_B0TrackerSeeds_params": "B0TrackerSeedParameters",
}
for relation, target in relation_targets.items():
    valid = a[relation + ".index"] >= 0
    if not ak.all(a[relation + ".collectionID"][valid] == collection_ids[target]):
        raise ValueError(f"{relation} does not reference expected collection {target}")
identities = list(
    zip(
        ak.to_list(a["EventHeader.runNumber"]), ak.to_list(a["EventHeader.eventNumber"])
    )
)
if any(len(r) != 1 or len(e) != 1 for r, e in identities):
    raise ValueError("Expected one EventHeader per event")
identities = [(r[0], e[0]) for r, e in identities]
if len(set(identities)) != len(identities):
    raise ValueError("Duplicate run/event identities")


def g(c, v, e):
    return ak.to_list(a[c + "." + v][e])


def rel(c, e):
    return g(c, "index", e)


rows = []
counts = {
    "seeds_total": 0,
    "ambiguous_hit_association": 0,
    "mixed_truth": 0,
    "nonprimary_proton": 0,
}
with open(args.output, "w") as out:
    for e in range(len(a)):
        pos = np.array([g("B0TrackerRecHits", "position." + k, e) for k in "xyz"]).T
        simmom = np.array([g("B0TrackerHits", "momentum." + k, e) for k in "xyz"]).T
        truth = np.array([g("MCParticles", "momentum." + k, e) for k in "xyz"]).T
        vertex = np.array([g("MCParticles", "vertex." + k, e) for k in "xyz"]).T
        raw = rel("_B0TrackerRecHits_rawHit", e)
        ar = rel("_B0TrackerRawHitAssociations_rawHit", e)
        ass = rel("_B0TrackerRawHitAssociations_simHit", e)
        particles = rel("_B0TrackerHits_particle", e)
        sh = rel("_B0TrackerSeeds_hits", e)
        pars = rel("_B0TrackerSeeds_params", e)
        # Collection IDs are checked above; validate indices before dereferencing.
        raw_to_sim = {}
        bounds = [
            (raw, len(g("B0TrackerRawHits", "cellID", e))),
            (ar, len(g("B0TrackerRawHits", "cellID", e))),
            (ass, len(simmom)),
            (particles, len(truth)),
            (sh, len(pos)),
            (pars, len(g("B0TrackerSeedParameters", "qOverP", e))),
        ]
        if any(any(i < 0 or i >= size for i in values) for values, size in bounds):
            raise ValueError(f"Invalid relation index in event {identities[e]}")
        for r, si in zip(ar, ass):
            raw_to_sim.setdefault(r, []).append(si)
        for s, (begin, end) in enumerate(
            zip(
                g("B0TrackerSeeds", "hits_begin", e), g("B0TrackerSeeds", "hits_end", e)
            )
        ):
            counts["seeds_total"] += 1
            if not 0 <= begin <= end <= len(sh) or end - begin < 3:
                raise ValueError(
                    f"Invalid seed hit span in event {identities[e]} seed {s}"
                )
            ids = sorted(
                sh[begin:end],
                key=lambda i: pos[i, 0] * np.sin(args.crossing_angle)
                + pos[i, 2] * np.cos(args.crossing_angle),
            )
            sims = []
            for h in ids:
                match = raw_to_sim.get(raw[h], [])
                if len(match) != 1 or match[0] < 0:
                    break
                sims.append(match[0])
            if len(sims) != len(ids):
                counts["ambiguous_hit_association"] += 1
                continue
            mc = [particles[i] for i in sims]
            if len(set(mc)) != 1 or mc[0] < 0:
                counts["mixed_truth"] += 1
                continue
            m = mc[0]
            if (
                g("MCParticles", "generatorStatus", e)[m] != 1
                or g("MCParticles", "PDG", e)[m] != 2212
            ):
                counts["nonprimary_proton"] += 1
                continue
            p = pos[ids]
            zion = p[:, 0] * np.sin(args.crossing_angle) + p[:, 2] * np.cos(
                args.crossing_angle
            )
            mid = min(
                range(1, len(p) - 1),
                key=lambda i: abs(zion[i] - (zion[0] + zion[-1]) / 2),
            )
            j = pars[s]
            q = g("B0TrackerSeedParameters", "qOverP", e)[j]
            ph = g("B0TrackerSeedParameters", "phi", e)[j]
            th = g("B0TrackerSeedParameters", "theta", e)[j]
            loc = g("B0TrackerSeedParameters", "loc.a", e)[j]
            z = g("B0TrackerSeedParameters", "loc.b", e)[j]
            r = [-loc * np.sin(ph), loc * np.cos(ph), z]
            d = [np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
            # Unweighted approximation of the current seeder's fitted midpoint query.
            # Coincides with its weighted fit only for equal per-hit variances.
            ca, sa = np.cos(args.crossing_angle), np.sin(args.crossing_angle)
            xion = p[:, 0] * ca - p[:, 2] * sa
            zref = zion.mean()
            scale = np.max(abs(zion - zref))
            u = (zion - zref) / scale
            um = ((zion[0] + zion[-1]) / 2 - zref) / scale
            ax = np.polynomial.polynomial.polyfit(u, xion, 2)
            ay = np.polynomial.polynomial.polyfit(u, p[:, 1], 1)
            xm = np.polynomial.polynomial.polyval(um, ax)
            ym = np.polynomial.polynomial.polyval(um, ay)
            zm = (zion[0] + zion[-1]) / 2
            fieldpoint = [xm * ca + zm * sa, ym, -xm * sa + zm * ca]
            out.write(
                " ".join(
                    map(
                        str, [len(rows), *p[0], *p[mid], *p[-1], *r, *d, *fieldpoint, q]
                    )
                )
                + "\n"
            )
            p_orig = np.linalg.norm(truth[m])
            hitmom = simmom[sims]
            mag = np.linalg.norm(hitmom, axis=1)
            # Diagnostic subset: exclude heavily degraded and backward/re-entering hits,
            # without selecting by either reconstructed estimate.
            clean = bool(np.all(mag > 0.5 * p_orig) and np.all(hitmom[:, 2] > 0))
            rows.append(
                dict(
                    event=e,
                    run_number=identities[e][0],
                    event_number=g("EventHeader", "eventNumber", e)[0],
                    seed=s,
                    nhits=len(ids),
                    mc=m,
                    clean_forward=clean,
                    triplet=[ids[0], ids[mid], ids[-1]],
                    q_current=q,
                    q_truth_origin=1 / p_orig,
                    q_truth_first=1 / mag[0],
                    direction_truth_first=(hitmom[0] / mag[0]).tolist(),
                    truth_momentum=truth[m].tolist(),
                    truth_vertex=vertex[m].tolist(),
                    seed_state=[loc, z, ph, th, q],
                    covariance=g(
                        "B0TrackerSeedParameters", "covariance.covariance[21]", e
                    )[j],
                )
            )
Path(args.output + ".json").write_text(json.dumps(rows))
counts.update(
    selected=len(rows),
    events=len(a),
    clean_forward=sum(r["clean_forward"] for r in rows),
)
counts.update(input=args.input, crossing_angle_rad=args.crossing_angle)
Path(args.output + ".counts.json").write_text(json.dumps(counts, indent=2))
print(counts)
