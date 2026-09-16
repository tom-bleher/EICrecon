#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Run physical-pitch AC-LGAD through all four existing local B0 comparison chains."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time

CHAINS = ('B0TelescopeCKF', 'B0TelescopeDirect', 'B0TelescopeTruthCKF', 'B0TelescopeTruthDirect')


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def simulation_inputs(profile: Path, receipt_path: Path) -> tuple[Path, Path, dict[str, str]]:
    """Require the completed simulation receipt from epic's profile runner."""
    profile = profile.resolve(strict=True)
    manifest_path = profile / 'b0_readout_contract.json'
    manifest = json.loads(manifest_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    if (manifest.get('schema_version'), manifest.get('kind'), manifest.get('response_version')) != (1, 'B0ACLGADReadout', 1):
        raise ValueError('Unsupported physical B0 profile')
    if (receipt.get('schema_version'), receipt.get('kind'), receipt.get('status')) != (1, 'B0ACLGADSimulation', 'completed'):
        raise ValueError('Require a completed B0 simulation receipt, not an arbitrary existing ROOT file')
    if Path(receipt['profile']).resolve() != profile or digest(manifest_path) != receipt['profile_sha256']:
        raise ValueError('Simulation/profile identity mismatch')
    files = manifest.get('files', {})
    if not files or manifest.get('compact') not in files or manifest.get('readout_file') not in files:
        raise ValueError('Incomplete B0 profile contract')
    inputs = {str(manifest_path): digest(manifest_path), str(receipt_path.resolve()): digest(receipt_path)}
    for name, expected in files.items():
        if Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError('Profile path traversal')
        path = (profile / name).resolve(strict=True)
        if not path.is_relative_to(profile) or not path.is_file() or digest(path) != expected:
            raise ValueError('Profile asset missing, changed, or outside root')
        inputs[str(path)] = expected
    present = {str(p.relative_to(profile)) for p in profile.rglob('*') if p.is_file() and p != manifest_path}
    if present != set(files):
        raise ValueError('Profile file set changed')
    simulation = Path(receipt['output']['path']).resolve(strict=True)
    if digest(simulation) != receipt['output']['sha256']:
        raise ValueError('Simulation output hash mismatch')
    inputs[str(simulation)] = receipt['output']['sha256']
    return simulation, profile / manifest['compact'], inputs


def command(executable: str, sim: Path, compact: Path, geometry_calibration: Path,
            material: Path, response: Path, out: Path, chains: tuple[str, ...],
            allow_experimental: bool, prior_scale: float, settings: list[str]) -> list[str]:
    if not chains or len(set(chains)) != len(chains) or any(c not in CHAINS for c in chains):
        raise ValueError('Choose unique supported comparison chains')
    if not math.isfinite(prior_scale) or prior_scale <= 0:
        raise ValueError('Prior scale must be positive and finite')
    collections = ['EventHeader', 'MCParticles', 'B0TrackerHits', 'B0ACLGADRawHits', 'B0ACLGADRawHitLinks',
        'B0ACLGADRawHitAssociations', 'B0ACLGADChannelHits', 'B0ACLGADMeasurements',
        'B0TelescopeSeeds', 'B0TelescopeSeedParameters', 'B0TelescopeTruthSeeds', 'B0TelescopeTruthSeedParameters']
    for stem in chains:
        collections.extend(stem + suffix for suffix in ('Tracks', 'Parameters', 'Trajectories', 'Associations', 'Links'))
    args = [executable, '-Pplugins=b0_telescope,b0_aclgad', '-Pb0_telescope:response=aclgad',
        f'-Pdd4hep:xml_files={compact}', f'-Pdd4hep:calib={geometry_calibration}', f'-Pacts:MaterialMap={material}',
        f'-Pb0_aclgad:B0ACLGADResponse:calibrationFile={response}',
        f'-Pb0_aclgad:B0ACLGADResponse:allowExperimental={str(allow_experimental).lower()}',
        f'-Ppodio:output_file={out / "reco.edm4eic.root"}',
        '-Ppodio:output_include_collections=' + ','.join(collections),
        f'-Pb0_telescope:B0TelescopeSeeds:stationMapFile={out / "stations.json"}',
        f'-Pb0_telescope:B0TelescopeSeeds:diagnosticsFile={out / "seeds.jsonl"}',
        f'-Pb0_telescope:B0TelescopeTruthSeeds:referenceFile={out / "truth.jsonl"}']
    for stem in chains:
        args.extend((f'-Pb0_telescope:{stem}:diagnosticsFile={out / (stem + ".jsonl")}',
                     f'-Pb0_telescope:{stem}:weakPriorScale={prior_scale}'))
    reserved = {'diagnosticsFile', 'referenceFile', 'stationMapFile', 'weakPriorScale'}
    allowed = set(CHAINS) | {'B0TelescopeSeeds', 'B0TelescopeTruthSeeds'}
    for setting in settings:
        key, separator, value = setting.partition('=')
        parts = key.split(':')
        if not separator or not value or len(parts) != 2 or parts[0] not in allowed or parts[1] in reserved:
            raise ValueError('Settings must be FactoryTag:parameter=value, excluding controlled output/prior settings')
        args.append('-Pb0_telescope:' + setting)
    return args + [str(sim)]


def main() -> None:
    # Reuse the provenance machinery from the preceding telescope benchmark PR.
    from provenance import describe, git_state, write, xml_dependencies
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('profile', 'simulation-receipt', 'geometry-calibration', 'material', 'response', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--chain', choices=CHAINS, action='append', default=[])
    parser.add_argument('--allow-experimental', action='store_true')
    parser.add_argument('--prior-scale', type=float, default=1)
    parser.add_argument('--set', dest='settings', action='append', default=[])
    parser.add_argument('--pdg', type=int, action='append', default=[])
    parser.add_argument('--source', type=Path, action='append', default=[])
    parser.add_argument('--asset', type=Path, action='append', default=[])
    parser.add_argument('--executable', default='eicrecon')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    profile = args.profile.resolve(strict=True)
    sim, compact, hashes = simulation_inputs(profile, args.simulation_receipt)
    executable = shutil.which(args.executable)
    if not executable:
        parser.error('Reconstruction executable not found')
    response = args.response.resolve(strict=True)
    calibration = args.geometry_calibration.resolve(strict=True)
    material = args.material.resolve(strict=True)
    out = args.out.absolute()
    if out.resolve().is_relative_to(profile):
        parser.error('Output must be outside the immutable profile')
    chains = tuple(args.chain or CHAINS)
    invocation = command(executable, sim, compact, calibration, material, response, out,
                         chains, args.allow_experimental, args.prior_scale, args.settings)
    for record in xml_dependencies(calibration):
        hashes[record['path']] = record['sha256']
    for path in [response, material, Path(executable), *args.asset]:
        path = path.resolve(strict=True); hashes[str(path)] = digest(path)
    manifest = {'schema_version': 1, 'status': 'planned', 'response': 'aclgad',
        'created_utc': datetime.now(timezone.utc).isoformat(), 'command': invocation,
        'chains': chains, 'response_calibration': describe(response),
        'inputs_sha256': hashes, 'source_checkouts': [git_state(p) for p in args.source],
        'limitations': ['The configured response is not automatically a calibrated detector model.',
                       'Explicit --asset inputs must include dynamic libraries and externally loaded calibrations.',
                       'Same final-refit settings do not prove covariance coverage or material correctness.']}
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'provenance.json', manifest)
    if args.dry_run:
        print(json.dumps(invocation, indent=2)); return
    started = time.perf_counter()
    try:
        manifest['status'] = 'running'; write(out / 'provenance.json', manifest)
        with (out / 'reconstruction.log').open('w') as log:
            subprocess.run(invocation, check=True, stdout=log, stderr=subprocess.STDOUT)
        for name, expected in hashes.items():
            if digest(Path(name)) != expected:
                raise ValueError(f'Input changed during reconstruction: {name}')
        simulation_inputs(profile, args.simulation_receipt)
        manifest['outputs'] = {name: describe(out / file) for name, file in (
            ('reco', 'reco.edm4eic.root'), ('truth', 'truth.jsonl'), ('stations', 'stations.json'), ('seeds', 'seeds.jsonl'))}
        manifest['status'] = 'completed'; manifest['wall_seconds'] = time.perf_counter() - started
        write(out / 'provenance.json', manifest)
    except Exception as error:
        manifest['status'] = 'failed'; manifest['error'] = str(error)
        write(out / 'provenance.json', manifest)
        raise
    analysis = [sys.executable, str(Path(__file__).with_name('aclgad_analysis.py')),
        '--reco', str(out / 'reco.edm4eic.root'), '--truth', str(out / 'truth.jsonl'),
        '--out', str(out / 'analysis'), '--provenance', str(out / 'provenance.json'),
        '--diagnostics', str(out / 'seeds.jsonl')]
    for stem in chains:
        analysis.extend(('--chain', stem, '--diagnostics', str(out / (stem + '.jsonl'))))
    for pdg in args.pdg:
        analysis.extend(('--pdg', str(pdg)))
    manifest['analysis_status'] = 'running'; write(out / 'provenance.json', manifest)
    try:
        subprocess.run(analysis, check=True)
        manifest['analysis_status'] = 'completed'
        manifest['analysis_summary'] = describe(out / 'analysis' / 'summary.json')
    except Exception as error:
        manifest['analysis_status'] = 'failed'; manifest['analysis_error'] = str(error)
        write(out / 'provenance.json', manifest)
        raise
    write(out / 'provenance.json', manifest)


if __name__ == '__main__':
    main()
