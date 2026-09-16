# SPDX-License-Identifier: LGPL-3.0-or-later
"""Content-addressed run provenance; no network calls or inferred source SHAs."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def describe(path: Path) -> dict:
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError(f'Not a regular file: {path}')
    return {'path': str(path), 'size': path.stat().st_size, 'sha256': digest(path)}


def xml_dependencies(path: Path) -> list[dict]:
    """Follow local include/ref references; fail on unresolved files, URLs or variables.

    This captures compact XML dependencies, not dynamically loaded libraries or
    files opened internally by detector plugins. Those require explicit assets.
    """
    seen: set[Path] = set()
    def visit(candidate: Path) -> None:
        candidate = candidate.resolve(strict=True)
        if candidate in seen:
            return
        seen.add(candidate)
        root = ET.parse(candidate).getroot()
        for node in root.iter():
            tag = node.tag.rsplit('}', 1)[-1].lower()
            if tag not in ('include', 'gdmlfile'):
                continue
            value = node.get('ref') or node.get('href') or node.get('name')
            if not value:
                raise ValueError(f'Include without a reference in {candidate}')
            value = os.path.expandvars(value)
            if '$' in value or '://' in value:
                raise ValueError(f'Unresolved/nonlocal include {value!r}')
            nested = Path(value)
            visit(nested if nested.is_absolute() else candidate.parent / nested)
    visit(path)
    return [describe(item) for item in sorted(seen)]


def git_state(path: Path) -> dict:
    def run(*args: str) -> str:
        return subprocess.run(['git', '-C', str(path), *args], check=True, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip()
    return {'path': str(path.resolve()), 'commit': run('rev-parse', 'HEAD'),
            'dirty': bool(run('status', '--porcelain'))}


def verify_file(record: dict, path: Path | None = None) -> None:
    path = path or Path(record['path'])
    if describe(path)['sha256'] != record['sha256']:
        raise ValueError(f'Content hash mismatch: {path}')


def verify_outputs(manifest: dict, reco: Path, truth: Path) -> None:
    if manifest.get('schema_version') != 1 or manifest.get('status') != 'completed':
        raise ValueError('Run manifest is not a completed supported run')
    verify_file(manifest['outputs']['reco'], reco)
    verify_file(manifest['outputs']['truth'], truth)


def write(path: Path, data: dict) -> None:
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)
