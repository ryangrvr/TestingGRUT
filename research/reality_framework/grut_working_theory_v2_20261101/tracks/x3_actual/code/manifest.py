#!/usr/bin/env python3
"""Seal/check X3 bytes and the exact inherited inputs. No physics is rerun."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRACKS = TRACK.parent
REPO = Path(subprocess.check_output(
    ['git', 'rev-parse', '--show-toplevel'], cwd=HERE, text=True).strip())
BASE = '9f71d34f8a768e4ec7bc60a46a65a833a818c1de'
TARGET = TRACK / 'MANIFEST.json'
INPUTS = [
    'e2e_kinetic/code/ek_common.py',
    'e2e_kinetic/code/ek_modes.c',
    'e2e_kinetic/code/ek_run.py',
    'e2e_kinetic/code/ek_ptable.py',
    'e2e_kinetic/results/ptable_summary.json',
    'mf_hybrid/results/bench_HFa.json',
    'higgs_reheat/code/sm_rge.py',
    'lf_fermion/code/lf1d.py',
    'lf_chain/results/handover_hohLF55_mid.json',
    'x2_strip/results/lower_H_lattice.json',
    'x2_strip/results/lower_H_diag.txt.gz',
]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def snapshot():
    inputs = []
    for name in INPUTS:
        path = TRACKS / name
        rel = str(path.relative_to(REPO))
        original = subprocess.check_output(['git', 'show', BASE + ':' + rel], cwd=REPO)
        assert original == path.read_bytes(), 'Inherited input changed: ' + rel
        blob = subprocess.check_output(
            ['git', 'rev-parse', BASE + ':' + rel], cwd=REPO, text=True).strip()
        inputs.append(dict(path=rel, sha256=sha(original), git_blob=blob))
    files = []
    for path in sorted(TRACK.rglob('*')):
        if not path.is_file() or path == TARGET or '__pycache__' in path.parts:
            continue
        files.append(dict(path=str(path.relative_to(REPO)),
                          bytes=path.stat().st_size, sha256=sha(path.read_bytes())))
    return dict(base_commit=BASE, inherited_inputs=inputs, package_files=files,
                scientific_acceptance=False,
                meaning='Byte identity and source provenance, not physical validity')


def check_engine_sources():
    original = (TRACKS / 'e2e_kinetic/code/ek_modes.c').read_text()
    half = original
    for old, new in [
        ('double lim = 0.02 / H;', 'double lim = 0.5 * 0.02 / H;'),
        ('0.01 * rhoP[P] / Pc[P]', '0.5 * 0.01 * rhoP[P] / Pc[P]'),
        ('0.2 * tk)', '0.5 * 0.2 * tk)'),
    ]:
        assert half.count(old) == 1
        half = half.replace(old, new)
    for path in (TRACK / 'results/kinetic').glob('*.json'):
        run = json.loads(path.read_text())
        source = half if run['settings']['dt_scale'] == 0.5 else original
        assert sha(source.encode()) == run['engine_sha256'], path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true', help='Seal the current completed packet')
    args = ap.parse_args()
    current = snapshot()
    check_engine_sources()
    if args.write:
        TARGET.write_text(json.dumps(current, indent=2) + '\n')
    else:
        assert json.loads(TARGET.read_text()) == current, 'Manifest mismatch'
    print(len(current['inherited_inputs']), 'unchanged inherited inputs;',
          len(current['package_files']), 'package files;',
          'all seven recorded engine source hashes match')


if __name__ == '__main__':
    main()
