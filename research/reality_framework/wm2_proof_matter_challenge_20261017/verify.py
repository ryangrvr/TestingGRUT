"""Archive checks by default; optional bounded diagnostics without rewriting it."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-checks', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    actual = {str(p.relative_to(HERE)) for p in HERE.iterdir() if p.is_file() and p.name != 'MANIFEST.json'}
    assert actual == set(manifest['files']), 'packet file inventory changed'
    for name, sha in manifest['files'].items():
        assert digest(HERE / name) == sha, name
    assert digest(HERE / 'CHARTER.md') == manifest['frozen_charter_sha256']
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    for name, sha in pins['sha256'].items():
        assert digest(ROOT / name) == sha, name
    executed = json.loads((HERE / 'EXECUTED_SOURCES.json').read_text())
    for name, sha in executed.items():
        assert digest(HERE / name) == sha, name
    checks = json.loads((HERE / 'CHECKS.json').read_text())
    findings = json.loads((HERE / 'FINDINGS.json').read_text())
    keys = ['charge_blocks', 'interaction_words', 'distinct_matter', 'receiver_dependence']
    assert all(checks[key]['status'] == 'PASS' for key in keys)
    assert findings['external_review'] == 'PENDING'
    assert findings['matter_universality_in_C'] == 'NOT_FORCED_COUNTEREXAMPLE'
    assert findings['new_formation_trajectories'] == 0
    assert findings['benchmark_law_changed'] is False
    assert checks['base_commit'] == pins['base_commit'] == findings['base_commit']
    print(f"PASS: {len(manifest['files'])} packet hashes, {len(pins['sha256'])} upstream pins, frozen charter and scoped statuses")
    if args.run_checks:
        from charge_controls import run as charge
        from dyson_controls import run as dyson
        from matter_controls import run as matter
        from resource_controls import run as resource
        for key, work in zip(keys, [charge, dyson, matter, resource]):
            fresh = work()
            assert fresh['status'] == 'PASS'
            print('REPLAY PASS:', key)
        assert digest(HERE / 'CHECKS.json') == manifest['files']['CHECKS.json']
        print('Archived results preserved; no formation solver invoked.')
    print('S9: AUTHOR_DERIVED_EXTERNAL_REVIEW_PENDING; C: COUNTEREXAMPLE; not a theorem checker')


if __name__ == '__main__':
    main()
