#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Offline B0 tracking efficiency, residual, pull, and D^2 analysis.

This script does not run reconstruction. It reads an EICrecon PODIO file and
optional seeder/CKF counter JSON, then writes a text report and a multi-page
PDF. Final-track pulls and D^2 are diagnostics until the independent refit
(PR 3) lands.

Physical B0 stations are clustered from ion-frame hit z with the same 50 mm
gap used by B0TrackerStubSeeder. A generator-stable proton is reconstructible
when those hits occupy at least three distinct stations.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

STATION_Z_GAP_MM = 50.0
MIN_STATIONS = 3
CROSSING_ANGLE = -0.025  # rad, DD4hep CrossingAngle
MATCH_WEIGHT_MIN = 0.5
PDG_PROTON = 2212


def wrap_phi(delta: float) -> float:
    return (delta + math.pi) % (2.0 * math.pi) - math.pi


def ion_z(x: float, z: float, crossing_angle: float = CROSSING_ANGLE) -> float:
    return x * math.sin(crossing_angle) + z * math.cos(crossing_angle)


def cluster_stations(z_values: list[float], gap: float = STATION_Z_GAP_MM) -> list[tuple[float, float]]:
    """Single-linkage clustering of ion-frame z; returns (z_min, z_max) intervals."""
    ordered = sorted(z for z in z_values if math.isfinite(z))
    if not ordered or not (gap > 0.0):
        return []
    stations: list[tuple[float, float]] = []
    z_min = z_max = ordered[0]
    for z in ordered[1:]:
        if z - z_max > gap:
            stations.append((z_min, z_max))
            z_min = z_max = z
        else:
            z_max = z
    stations.append((z_min, z_max))
    return stations


def perigee_from_ray(ref, direction, origin=(0.0, 0.0, 0.0)):
    """Straight-line perigee (loc0, loc1, phi, theta) at the origin, matching the seeder."""
    dx, dy, dz = direction
    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if not (norm > 0.0):
        return 0.0, 0.0, 0.0, 0.0
    dx, dy, dz = dx / norm, dy / norm, dz / norm
    rx = ref[0] - origin[0]
    ry = ref[1] - origin[1]
    denom = dx * dx + dy * dy
    t = 0.0 if denom <= 0.0 else -(rx * dx + ry * dy) / denom
    x_pca = rx + t * dx
    y_pca = ry + t * dy
    z_pca = ref[2] + t * dz - origin[2]
    phi = math.atan2(dy, dx)
    theta = math.acos(max(-1.0, min(1.0, dz)))
    loc0 = -x_pca * math.sin(phi) + y_pca * math.cos(phi)
    return loc0, z_pca, phi, theta


def mahalanobis_d2(delta, cov5) -> float | None:
    """D^2 = d^T C^{-1} d for a 5-vector. None if C is not invertible."""
    try:
        import numpy as np
    except ImportError:
        return None
    cov = np.array(cov5, dtype=float)
    if cov.shape != (5, 5):
        return None
    try:
        return float(np.dot(delta, np.linalg.solve(cov, delta)))
    except np.linalg.LinAlgError:
        return None


def _self_test() -> int:
    official = [5902.0, 6172.0, 6442.0, 6712.0]
    assert len(cluster_stations(official)) == 4
    faces = [5902.0, 5909.0, 6172.0, 6179.0, 6442.0, 6449.0, 6712.0, 6719.0]
    assert len(cluster_stations(faces)) == 4
    assert abs(wrap_phi(0.1) - 0.1) < 1e-12
    assert abs(wrap_phi(math.pi + 0.2) + (math.pi - 0.2)) < 1e-9
    loc0, loc1, phi, theta = perigee_from_ray((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    assert abs(loc0) < 1e-12 and abs(loc1) < 1e-12
    assert abs(theta) < 1e-12
    identity = [[1.0 if i == j else 0.0 for j in range(5)] for i in range(5)]
    d2 = mahalanobis_d2([1.0, 0.0, 0.0, 0.0, 0.0], identity)
    assert d2 is None or abs(d2 - 1.0) < 1e-9
    print("b0 tracking analysis self-test passed")
    return 0


def _particle_id(particle) -> tuple[int, int]:
    oid = particle.getObjectID()
    return int(oid.collectionID), int(oid.index)


def _momentum(particle):
    p = particle.getMomentum()
    return float(p.x), float(p.y), float(p.z)


def _vertex(particle):
    v = particle.getVertex()
    return float(v.x), float(v.y), float(v.z)


def _p_theta_phi(px, py, pz):
    p = math.sqrt(px * px + py * py + pz * pz)
    theta = math.acos(max(-1.0, min(1.0, pz / p))) if p > 0.0 else 0.0
    phi = math.atan2(py, px)
    return p, theta, phi


def _cov5(track_params):
    cov = track_params.getCovariance()
    return [[float(cov(i, j)) for j in range(5)] for i in range(5)]


def _params5(track_params):
    loc = track_params.getLoc()
    return [
        float(loc.a),
        float(loc.b),
        float(track_params.getPhi()),
        float(track_params.getTheta()),
        float(track_params.getQOverP()),
    ]


def _truth_params(particle):
    px, py, pz = _momentum(particle)
    p, theta, phi = _p_theta_phi(px, py, pz)
    charge = float(particle.getCharge())
    q_over_p = charge / p if p > 0.0 else 0.0
    loc0, loc1, pphi, ptheta = perigee_from_ray(_vertex(particle), (px, py, pz))
    return [loc0, loc1, pphi, ptheta, q_over_p], p, theta, phi


@dataclass
class ParticleInfo:
    pid: tuple[int, int]
    p: float
    theta: float
    phi: float
    theta_ion: float
    n_stations: int
    charge: float
    truth5: list[float]
    reconstructible: bool


@dataclass
class Sample:
    particles: list[ParticleInfo] = field(default_factory=list)
    n_events: int = 0
    n_stub_seeds: int = 0
    n_stub_tracks: int = 0
    n_truth_tracks: int = 0
    stub_fakes: int = 0
    truth_fakes: int = 0
    stub_duplicates: int = 0
    truth_duplicates: int = 0
    seed_residuals: list = field(default_factory=list)
    stub_residuals: list = field(default_factory=list)
    truth_residuals: list = field(default_factory=list)


def _get_collection(frame, name):
    try:
        return frame.get(name)
    except Exception:
        return None


def _load_podio(path: Path):
    from podio.reading import get_reader

    reader = get_reader(str(path))
    return reader.get("events")


def _hit_particle_id(hit) -> tuple[int, int] | None:
    try:
        particle = hit.getParticle()
    except Exception:
        return None
    if particle is None:
        return None
    try:
        return _particle_id(particle)
    except Exception:
        return None


def _assoc_weight_map(assocs):
    best: dict[tuple[int, int], tuple[float, tuple[int, int]]] = {}
    if assocs is None:
        return best
    for assoc in assocs:
        rec = assoc.getRec()
        sim = assoc.getSim()
        weight = float(assoc.getWeight())
        rec_id = (int(rec.getObjectID().collectionID), int(rec.getObjectID().index))
        sim_id = _particle_id(sim)
        prev = best.get(rec_id)
        if prev is None or weight > prev[0]:
            best[rec_id] = (weight, sim_id)
    return best


def _seed_mc_id(seed, raw_assocs) -> tuple[int, int] | None:
    votes: dict[tuple[int, int], int] = defaultdict(int)
    try:
        hits = seed.getHits()
    except Exception:
        return None
    for hit in hits:
        try:
            raw = hit.getRawHit()
        except Exception:
            raw = None
        if raw_assocs is not None and raw is not None:
            for assoc in raw_assocs:
                try:
                    if assoc.getRawHit() == raw:
                        votes[_particle_id(assoc.getSimHit().getParticle())] += 1
                except Exception:
                    continue
        else:
            pid = _hit_particle_id(hit)
            if pid is not None:
                votes[pid] += 1
    if not votes:
        return None
    return max(votes.items(), key=lambda kv: kv[1])[0]


def _residual_record(truth: ParticleInfo, reco5, cov5, chi2ndf, n_stations_reco):
    dphi = wrap_phi(reco5[2] - truth.truth5[2])
    delta = [
        reco5[0] - truth.truth5[0],
        reco5[1] - truth.truth5[1],
        dphi,
        reco5[3] - truth.truth5[3],
        reco5[4] - truth.truth5[4],
    ]
    pulls = []
    if cov5 is not None:
        for i, d in enumerate(delta):
            var = cov5[i][i]
            pulls.append(d / math.sqrt(var) if var > 0.0 and math.isfinite(var) else float("nan"))
    else:
        pulls = [float("nan")] * 5
    q_over_p_t = truth.truth5[4]
    p_reco = abs(1.0 / reco5[4]) if reco5[4] != 0.0 else float("nan")
    dp_over_p = (p_reco - truth.p) / truth.p if truth.p > 0.0 else float("nan")
    return {
        "p": truth.p,
        "theta": truth.theta,
        "theta_ion": truth.theta_ion,
        "n_stations": truth.n_stations,
        "n_stations_reco": n_stations_reco,
        "dloc0": delta[0],
        "dloc1": delta[1],
        "dphi": dphi,
        "dtheta": delta[3],
        "dqoverp": delta[4],
        "dp_over_p": dp_over_p,
        "chi2ndf": chi2ndf,
        "pulls": pulls,
        "d2": mahalanobis_d2(delta, cov5) if cov5 is not None else None,
        "charge_ok": (reco5[4] == 0.0) or (math.copysign(1.0, reco5[4]) == math.copysign(1.0, q_over_p_t)),
    }


def analyze_events(events) -> Sample:
    sample = Sample()
    for frame in events:
        sample.n_events += 1
        mc = _get_collection(frame, "MCParticles")
        sim_hits = _get_collection(frame, "B0TrackerHits")
        seeds = _get_collection(frame, "B0TrackerSeeds")
        raw_assocs = _get_collection(frame, "B0TrackerRawHitAssociations")
        stub_tracks = _get_collection(frame, "B0TrackerCKFTracks")
        stub_assocs = _get_collection(frame, "B0TrackerCKFTrackAssociations")
        truth_tracks = _get_collection(frame, "B0TrackerCKFTruthSeededTracks")
        truth_assocs = _get_collection(frame, "B0TrackerCKFTruthSeededTrackAssociations")
        if mc is None:
            continue

        hits_by_particle: dict[tuple[int, int], list[float]] = defaultdict(list)
        if sim_hits is not None:
            for hit in sim_hits:
                pid = _hit_particle_id(hit)
                if pid is None:
                    continue
                pos = hit.getPosition()
                hits_by_particle[pid].append(ion_z(float(pos.x), float(pos.z)))

        event_particles: dict[tuple[int, int], ParticleInfo] = {}
        for particle in mc:
            if int(particle.getGeneratorStatus()) != 1:
                continue
            if abs(int(particle.getPDG())) != PDG_PROTON:
                continue
            pid = _particle_id(particle)
            px, py, pz = _momentum(particle)
            p, theta, phi = _p_theta_phi(px, py, pz)
            stations = cluster_stations(hits_by_particle.get(pid, []))
            truth5, _, _, _ = _truth_params(particle)
            tx = px * math.cos(CROSSING_ANGLE) - pz * math.sin(CROSSING_ANGLE)
            tz = px * math.sin(CROSSING_ANGLE) + pz * math.cos(CROSSING_ANGLE)
            theta_ion = math.atan2(math.hypot(tx, py), tz)
            info = ParticleInfo(
                pid=pid,
                p=p,
                theta=theta,
                phi=phi,
                theta_ion=theta_ion,
                n_stations=len(stations),
                charge=float(particle.getCharge()),
                truth5=truth5,
                reconstructible=len(stations) >= MIN_STATIONS,
            )
            event_particles[pid] = info
            sample.particles.append(info)

        if seeds is not None:
            sample.n_stub_seeds += len(seeds)
            for seed in seeds:
                mc_id = _seed_mc_id(seed, raw_assocs)
                info = event_particles.get(mc_id) if mc_id is not None else None
                if info is None or not info.reconstructible:
                    continue
                info.__dict__["_seeded"] = True
                try:
                    params = seed.getParams()
                    sample.seed_residuals.append(
                        _residual_record(info, _params5(params), _cov5(params), float("nan"), info.n_stations)
                    )
                except Exception:
                    pass

        stub_best = _assoc_weight_map(stub_assocs)
        truth_best = _assoc_weight_map(truth_assocs)

        def consume_tracks(tracks, best, flag, residuals):
            n = 0 if tracks is None else len(tracks)
            if tracks is None:
                return n, 0, 0
            fake = 0
            counts: dict[tuple[int, int], int] = defaultdict(int)
            for track in tracks:
                rec_id = (int(track.getObjectID().collectionID), int(track.getObjectID().index))
                hit = best.get(rec_id)
                matched_pid = hit[1] if hit is not None and hit[0] >= MATCH_WEIGHT_MIN else None
                info = event_particles.get(matched_pid) if matched_pid is not None else None
                if info is not None and info.reconstructible:
                    info.__dict__[flag] = True
                    counts[matched_pid] += 1
                    try:
                        traj = track.getTrajectory()
                        pars = traj.getTrackParameters()[0]
                        ndf = float(track.getNdf())
                        chi2ndf = float(track.getChi2()) / ndf if ndf > 0.0 else float("nan")
                        residuals.append(
                            _residual_record(info, _params5(pars), _cov5(pars), chi2ndf, info.n_stations)
                        )
                    except Exception:
                        pass
                else:
                    fake += 1
            duplicates = sum(max(0, c - 1) for c in counts.values())
            return n, fake, duplicates

        n_stub, f_stub, d_stub = consume_tracks(
            stub_tracks, stub_best, "_stub_matched", sample.stub_residuals
        )
        n_truth, f_truth, d_truth = consume_tracks(
            truth_tracks, truth_best, "_truth_matched", sample.truth_residuals
        )
        sample.n_stub_tracks += n_stub
        sample.n_truth_tracks += n_truth
        sample.stub_fakes += f_stub
        sample.truth_fakes += f_truth
        sample.stub_duplicates += d_stub
        sample.truth_duplicates += d_truth

    return sample


def _summarize(sample: Sample) -> dict:
    reconstructible = [p for p in sample.particles if p.reconstructible]
    seeded = [p for p in reconstructible if p.__dict__.get("_seeded")]
    n_rec = len(reconstructible)
    n_seeded = len(seeded)
    n_stub = sum(1 for p in reconstructible if p.__dict__.get("_stub_matched"))
    n_stub_given_seed = sum(1 for p in seeded if p.__dict__.get("_stub_matched"))
    n_truth = sum(1 for p in reconstructible if p.__dict__.get("_truth_matched"))

    def rate(num, den):
        return None if den == 0 else num / den

    stub_duplicates = sample.stub_duplicates
    truth_duplicates = sample.truth_duplicates

    def residual_stats(rows, key):
        vals = [r[key] for r in rows if r.get(key) is not None and math.isfinite(r[key])]
        if not vals:
            return None
        mean = sum(vals) / len(vals)
        var = sum((v - mean) ** 2 for v in vals) / len(vals)
        return {"n": len(vals), "mean": mean, "rms": math.sqrt(var)}

    def pull_stats(rows, index):
        vals = [r["pulls"][index] for r in rows if r.get("pulls") and math.isfinite(r["pulls"][index])]
        if not vals:
            return None
        mean = sum(vals) / len(vals)
        var = sum((v - mean) ** 2 for v in vals) / len(vals)
        return {"n": len(vals), "mean": mean, "width": math.sqrt(var)}

    charge_ok = [r for r in sample.stub_residuals if r.get("charge_ok")]
    return {
        "n_events": sample.n_events,
        "n_truth_protons": len(sample.particles),
        "n_reconstructible": n_rec,
        "n_seeded": n_seeded,
        "n_stub_seeds": sample.n_stub_seeds,
        "n_stub_tracks": sample.n_stub_tracks,
        "n_truth_tracks": sample.n_truth_tracks,
        "epsilon_seed": rate(n_seeded, n_rec),
        "epsilon_track_given_seed": rate(n_stub_given_seed, n_seeded),
        "epsilon_total_stub": rate(n_stub, n_rec),
        "epsilon_total_truth": rate(n_truth, n_rec),
        "fake_rate_stub": rate(sample.stub_fakes, sample.n_stub_tracks),
        "fake_rate_truth": rate(sample.truth_fakes, sample.n_truth_tracks),
        "duplicate_rate_stub": rate(stub_duplicates, n_rec),
        "duplicate_rate_truth": rate(truth_duplicates, n_rec),
        "charge_correct_stub": rate(len(charge_ok), len(sample.stub_residuals)),
        "seed_dp_over_p": residual_stats(sample.seed_residuals, "dp_over_p"),
        "stub_dp_over_p": residual_stats(sample.stub_residuals, "dp_over_p"),
        "stub_dtheta": residual_stats(sample.stub_residuals, "dtheta"),
        "stub_dphi": residual_stats(sample.stub_residuals, "dphi"),
        "seed_pull_qoverp": pull_stats(sample.seed_residuals, 4),
        "stub_pull_qoverp": pull_stats(sample.stub_residuals, 4),
        "stub_d2": residual_stats(
            [{"d2": r["d2"]} for r in sample.stub_residuals if r.get("d2") is not None], "d2"
        ),
        "by_stations": {
            str(n): {
                "reconstructible": sum(1 for p in reconstructible if p.n_stations == n),
                "seeded": sum(1 for p in seeded if p.n_stations == n),
                "stub_matched": sum(
                    1 for p in reconstructible if p.n_stations == n and p.__dict__.get("_stub_matched")
                ),
            }
            for n in sorted({p.n_stations for p in reconstructible})
        },
    }


def _fmt(value) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, dict) and "mean" in value:
        width = value.get("width", value.get("rms"))
        return f"{value['mean']:.4g} (rms {width:.4g}, n={value['n']})"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def write_report(path: Path, summary: dict, counters: dict | None, provenance: dict | None) -> None:
    lines = ["# B0 tracking benchmark", ""]
    if provenance:
        lines += ["## Provenance", ""]
        for key, value in provenance.items():
            lines.append(f"- **{key}**: `{value}`")
        lines.append("")
        lines.append("Quoted config-comment yields (698/765 stub, 723/765 truth) are **historical**")
        lines.append("until this command reproduces them on the pinned sample.")
        lines.append("")
    lines += ["## Efficiency", ""]
    for key in (
        "n_events",
        "n_truth_protons",
        "n_reconstructible",
        "n_seeded",
        "epsilon_seed",
        "epsilon_track_given_seed",
        "epsilon_total_stub",
        "epsilon_total_truth",
        "fake_rate_stub",
        "duplicate_rate_stub",
        "charge_correct_stub",
    ):
        lines.append(f"- **{key}**: {_fmt(summary.get(key))}")
    lines += ["", "## Residuals and diagnostic pulls", ""]
    for key in (
        "seed_dp_over_p",
        "stub_dp_over_p",
        "stub_dtheta",
        "stub_dphi",
        "seed_pull_qoverp",
        "stub_pull_qoverp",
        "stub_d2",
    ):
        lines.append(f"- **{key}**: {_fmt(summary.get(key))}")
    lines.append("")
    lines.append("Final-track pulls and D^2 are diagnostics, not validated physics results.")
    lines += ["", "## By physical B0 stations", ""]
    for n, row in summary.get("by_stations", {}).items():
        lines.append(
            f"- **{n} stations**: reconstructible={row['reconstructible']}, "
            f"seeded={row['seeded']}, stub matched={row['stub_matched']}"
        )
    if counters:
        lines += ["", "## Failure accounting", ""]
        seeder = counters.get("seeder", {})
        lines.append("### Seeder")
        for key, value in seeder.items():
            lines.append(f"- **{key}**: {value}")
        if "overlapRejected" in counters:
            lines.append(f"- **overlapRejected**: {counters['overlapRejected']}")
        for chain_name, payload in counters.get("chains", {}).items():
            lines.append(f"### CKF / ambiguity ({chain_name})")
            for key, value in payload.get("ckf", {}).items():
                lines.append(f"- **{key}**: {value}")
            for key, value in payload.get("ambiguity", {}).items():
                lines.append(f"- **{key}**: {value}")
    path.write_text("\n".join(lines) + "\n")


def write_plots(path: Path, sample: Sample) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping PDF plots", file=sys.stderr)
        return

    reconstructible = [p for p in sample.particles if p.reconstructible]
    seeded_mask = [bool(p.__dict__.get("_seeded")) for p in reconstructible]
    stub_mask = [bool(p.__dict__.get("_stub_matched")) for p in reconstructible]

    def efficiency_xy(values, matched_mask, bins):
        import numpy as np

        values = np.asarray(values)
        matched_mask = np.asarray(matched_mask, dtype=bool)
        tot, edges = np.histogram(values, bins=bins)
        ok, _ = np.histogram(values[matched_mask], bins=edges)
        with np.errstate(divide="ignore", invalid="ignore"):
            eff = np.true_divide(ok, tot)
            eff[tot == 0] = float("nan")
        centers = 0.5 * (edges[:-1] + edges[1:])
        return centers, eff, tot

    with PdfPages(path) as pdf:
        if reconstructible:
            import numpy as np

            pvals = [p.p for p in reconstructible]
            tvals = [1000.0 * p.theta_ion for p in reconstructible]
            fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
            for ax, vals, bins, xlabel, mask, label in (
                (axes[0], pvals, np.linspace(8.0, 41.0, 12), r"$p$ [GeV]", seeded_mask, r"$\epsilon_{\mathrm{seed}}$"),
                (axes[1], tvals, np.linspace(4.0, 22.0, 10), r"$\theta_{\mathrm{ion}}$ [mrad]", seeded_mask, r"$\epsilon_{\mathrm{seed}}$"),
            ):
                x, y, _ = efficiency_xy(vals, mask, bins)
                ax.plot(x, y, marker="o", label=label)
                x2, y2, _ = efficiency_xy(vals, stub_mask, bins)
                ax.plot(x2, y2, marker="s", label=r"$\epsilon_{\mathrm{total}}$")
                ax.set_xlabel(xlabel)
                ax.set_ylabel("efficiency")
                ax.set_ylim(0.0, 1.05)
                ax.legend()
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        def hist_page(rows, key, xlabel, title):
            vals = [r[key] for r in rows if r.get(key) is not None and math.isfinite(r[key])]
            fig, ax = plt.subplots(figsize=(6.5, 4.5))
            if vals:
                ax.hist(vals, bins=40, histtype="step")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("tracks")
            ax.set_title(title)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        hist_page(sample.stub_residuals, "dp_over_p", r"$\Delta p/p$", r"Stub-seeded $\Delta p/p$")
        hist_page(sample.stub_residuals, "dqoverp", r"$\Delta(q/p)$ [1/GeV]", r"Stub-seeded $\Delta(q/p)$")
        hist_page(sample.stub_residuals, "dtheta", r"$\Delta\theta$ [rad]", r"Stub-seeded $\Delta\theta$")
        hist_page(sample.stub_residuals, "dphi", r"$\Delta\phi$ [rad]", r"Stub-seeded $\Delta\phi$")
        hist_page(sample.stub_residuals, "chi2ndf", r"$\chi^2/\mathrm{ndf}$", r"Stub-seeded $\chi^2/\mathrm{ndf}$")
        hist_page(sample.seed_residuals, "dtheta", r"$\Delta\theta$ [rad]", r"Seed $\Delta\theta$")
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
        for ax, rows, title in (
            (axes[0], sample.seed_residuals, "Seed $q/p$ pull (diagnostic)"),
            (axes[1], sample.stub_residuals, "Final $q/p$ pull (diagnostic)"),
        ):
            vals = [r["pulls"][4] for r in rows if r.get("pulls") and math.isfinite(r["pulls"][4])]
            if vals:
                ax.hist(vals, bins=40, histtype="step")
            ax.set_xlabel(r"$(q/p-\mathrm{truth})/\sigma$")
            ax.set_title(title)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        hist_page(
            [{"d2": r["d2"]} for r in sample.stub_residuals if r.get("d2") is not None],
            "d2",
            r"$D^2$",
            "Final-track $D^2$ (diagnostic)",
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reco", type=Path, help="EICrecon PODIO file (.edm4eic.root)")
    parser.add_argument("--counters", type=Path, help="JSON from B0_TRACKING_COUNTERS_FILE")
    parser.add_argument("--provenance", type=Path, help="JSON provenance written by the run script")
    parser.add_argument("--out", type=Path, required=False, default=Path("b0_tracking_benchmark"))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)

    if args.self_test:
        return _self_test()
    if args.reco is None:
        parser.error("--reco is required unless --self-test is set")

    args.out.mkdir(parents=True, exist_ok=True)
    events = _load_podio(args.reco)
    sample = analyze_events(events)
    summary = _summarize(sample)
    counters = json.loads(args.counters.read_text()) if args.counters and args.counters.exists() else None
    provenance = json.loads(args.provenance.read_text()) if args.provenance and args.provenance.exists() else None
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(args.out / "report.md", summary, counters, provenance)
    write_plots(args.out / "b0_tracking_validation.pdf", sample)
    print(f"wrote {args.out / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
