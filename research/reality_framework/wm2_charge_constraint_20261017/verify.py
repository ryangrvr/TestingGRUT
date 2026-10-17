"""Integrity and scoped dispositions; optional controls, not a theorem checker."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-checks', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'MANIFEST.json').read_text())
    names = {p.name for p in HERE.iterdir() if p.is_file() and p.name != 'MANIFEST.json'}
    assert names == set(manifest['files'])
    for name, digest in manifest['files'].items():
        assert sha(HERE/name) == digest, name
    assert sha(HERE/'CHARTER.md') == manifest['frozen_charter_sha256']
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    for name, digest in pins['sha256'].items():
        assert sha(ROOT/name) == digest, name
    executed = json.loads((HERE/'EXECUTED_SOURCES.json').read_text())
    for name, digest in executed.items():
        assert sha(HERE/name) == digest, name
    checks = json.loads((HERE/'CHECKS.json').read_text())
    model_path = ROOT/'research/reality_framework/wm2_prepared_disorder_20261016/model.py'
    assert sha(model_path) == checks['frozen_S4_implementation_sha256']
    findings = json.loads((HERE/'FINDINGS.json').read_text())
    assert findings['full_universal_locality_gate'] == 'FAIL_ALLOWED_INTERNAL_CHANNEL'
    assert findings['constraint_selection'] == 'SUPPLIED_NOT_DERIVED'
    assert findings['external_review'] == 'PENDING'
    assert findings['new_formation_trajectories'] == 0
    assert findings['baseline_WM2_modified'] is False
    assert checks['status'] == 'PASS_DIAGNOSTICS_FULL_GATE_NEGATIVE'
    assert all(row['formation_component_leakage'] == 0 for row in checks['kernel_isometry'])
    assert all(row['Gauss_commutator_frobenius'] > 0 or row['flux_support_commutator_frobenius'] > 0 for row in checks['witnesses'])
    print(f"PASS: {len(names)} packet hashes, frozen charter, upstream/executed sources and negative full-gate disposition")
    if args.run_checks:
        from checks import run
        fresh = run()
        assert fresh['status'] == checks['status']
        assert [r['complete_kernel_dimension'] for r in fresh['kernel_isometry']] == [141,1776,24780]
        assert sha(HERE/'CHECKS.json') == manifest['files']['CHECKS.json']
        print('REPLAY PASS: witness, complete-kernel, local/weighted response and cut-count diagnostics; archive preserved')
    print('CHARGED_POSITION: AUTHOR_DERIVED_IN_CLASS; UNIVERSAL_LOCALITY: FAIL; INDEPENDENT_REVIEW: PENDING')


if __name__ == '__main__':
    main()
