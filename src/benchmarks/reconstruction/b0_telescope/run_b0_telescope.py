#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Run explicitly selected B0 comparison chains and hash inputs/outputs."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time

from analyze_b0_telescope import DEFAULT_CHAINS
from provenance import describe, git_state, verify_file, write, xml_dependencies


def command(executable: str, sim: Path, compact: Path, calibration: Path, material: Path,
            out: Path, chains: tuple[str, ...], scale: float, settings: list[str]) -> list[str]:
    if not chains or len(set(chains)) != len(chains) or any(c not in DEFAULT_CHAINS for c in chains):
        raise ValueError('Choose unique supported comparison chains')
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError('Weak-prior scale must be finite and positive')
    collections = ['EventHeader', 'MCParticles', 'B0TrackerHits', 'B0TrackerRawHits',
        'B0TrackerRecHits', 'B0TrackerRawHitAssociations', 'B0TrackerMeasurements',
        'B0TelescopeSeeds', 'B0TelescopeSeedParameters',
        'B0TelescopeTruthSeeds', 'B0TelescopeTruthSeedParameters']
    for stem in chains:
        collections += [stem + suffix for suffix in ('Tracks', 'Parameters', 'Trajectories', 'Associations', 'Links')]
    args = [executable, '-Pplugins=b0_telescope', f'-Pdd4hep:xml_files={compact}',
        f'-Pdd4hep:calib={calibration}', f'-Pacts:MaterialMap={material}',
        f'-Ppodio:output_file={out / "reco.edm4eic.root"}',
        '-Ppodio:output_include_collections=' + ','.join(collections),
        f'-Pb0_telescope:B0TelescopeSeeds:stationMapFile={out / "stations.json"}',
        f'-Pb0_telescope:B0TelescopeSeeds:diagnosticsFile={out / "seeds.jsonl"}',
        f'-Pb0_telescope:B0TelescopeTruthSeeds:referenceFile={out / "truth.jsonl"}']
    for stem in chains:
        args += [f'-Pb0_telescope:{stem}:diagnosticsFile={out / (stem + ".jsonl")}',
                 f'-Pb0_telescope:{stem}:weakPriorScale={scale}']
    reserved = ('diagnosticsFile', 'referenceFile', 'stationMapFile', 'weakPriorScale')
    allowed = set(DEFAULT_CHAINS) | {'B0TelescopeSeeds', 'B0TelescopeTruthSeeds'}
    for setting in settings:
        key, separator, value = setting.partition('=')
        parts = key.split(':')
        if not separator or not value or len(parts) != 2 or parts[0] not in allowed or parts[1] in reserved:
            raise ValueError('Settings must be FactoryTag:parameter=value, excluding controlled paths/prior scale')
        args.append('-Pb0_telescope:' + setting)
    return args + [str(sim)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('sim', 'compact', 'calibration', 'material', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--chain', choices=DEFAULT_CHAINS, action='append', default=[])
    parser.add_argument('--weak-prior-scale', type=float, default=1)
    parser.add_argument('--set', dest='settings', action='append', default=[])
    parser.add_argument('--pdg', type=int, action='append', default=[])
    parser.add_argument('--executable', default='eicrecon')
    parser.add_argument('--source', type=Path, action='append', default=[], help='Source checkout to record, not assumed to match binary')
    parser.add_argument('--asset', type=Path, action='append', default=[], help='Additional detector/calibration/plugin input')
    parser.add_argument('--dry-run', action='store_true', help='Record command/provenance without executing reconstruction')
    args = parser.parse_args()
    paths = {name: getattr(args, name).resolve(strict=True) for name in ('sim', 'compact', 'calibration', 'material')}
    executable = shutil.which(args.executable)
    if not executable:
        parser.error(f'Executable not found: {args.executable}')
    out = args.out.resolve()
    chains = tuple(args.chain or DEFAULT_CHAINS)
    invocation = command(executable, **paths, out=out, chains=chains,
                         scale=args.weak_prior_scale, settings=args.settings)
    inputs = {name: describe(path) for name, path in paths.items()}
    dependencies = {r['path']: r for root in (paths['compact'], paths['calibration']) for r in xml_dependencies(root)}
    for path in args.asset:
        record = describe(path); dependencies[record['path']] = record
    manifest = {'schema_version': 1, 'status': 'planned',
        'created_utc': datetime.now(timezone.utc).isoformat(), 'command': invocation,
        'chains': chains, 'inputs': inputs, 'xml_and_explicit_assets': list(dependencies.values()),
        'executable': describe(Path(executable)), 'source_checkouts': [git_state(p) for p in args.source],
        'limitations': ['Source checkout identity does not prove installed-library identity.',
                       'Dynamic plugin and calibration files must be supplied as --asset to be hashed.',
                       'Requesting all chains is a comparison run, not single-chain timing.']}
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'provenance.json', manifest)
    if args.dry_run:
        print(json.dumps(invocation, indent=2)); return 0
    started = time.perf_counter()
    try:
        manifest['status'] = 'running'; write(out / 'provenance.json', manifest)
        with (out / 'reconstruction.log').open('w') as log:
            subprocess.run(invocation, check=True, stdout=log, stderr=subprocess.STDOUT)
        for record in [*inputs.values(), *dependencies.values(), manifest['executable']]:
            verify_file(record)
        manifest['outputs'] = {name: describe(out / file) for name, file in (
            ('reco', 'reco.edm4eic.root'), ('truth', 'truth.jsonl'), ('stations', 'stations.json'), ('seeds', 'seeds.jsonl'))}
        manifest['status'] = 'completed'; manifest['reconstruction_wall_seconds'] = time.perf_counter() - started
        write(out / 'provenance.json', manifest)
    except Exception as error:
        manifest['status'] = 'failed'; manifest['error'] = str(error)
        manifest['reconstruction_wall_seconds'] = time.perf_counter() - started
        write(out / 'provenance.json', manifest)
        raise
    analysis = [sys.executable, str(Path(__file__).with_name('analyze_b0_telescope.py')),
        '--reco', str(out / 'reco.edm4eic.root'), '--truth', str(out / 'truth.jsonl'),
        '--out', str(out / 'analysis'), '--provenance', str(out / 'provenance.json'),
        '--diagnostics', str(out / 'seeds.jsonl')]
    for stem in chains:
        analysis += ['--chain', stem, '--diagnostics', str(out / (stem + '.jsonl'))]
    for pdg in args.pdg:
        analysis += ['--pdg', str(pdg)]
    subprocess.run(analysis, check=True)
    print(out / 'analysis' / 'report.md')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
