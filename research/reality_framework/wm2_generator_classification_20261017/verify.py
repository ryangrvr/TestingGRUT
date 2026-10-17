"""Verify packet integrity and optionally replay finite diagnostics.

No archive mutation and no formation workers. This is not proof certification.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-checks', action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    repo = here.parents[2]
    manifest = json.loads((here / 'MANIFEST.json').read_text())
    for name, expected in manifest['sha256'].items():
        assert digest(here / name) == expected, f'Packet changed: {name}'
    pins = json.loads((here / 'SOURCE_PINS.json').read_text())
    assert digest(here / 'CHARTER.md') == pins['charter_sha256']
    for name, expected in pins['files'].items():
        assert digest(repo / name) == expected, f'Upstream changed: {name}'
    executed = json.loads((here / 'EXECUTED_SOURCES.json').read_text())
    for name, expected in executed['sha256'].items():
        assert digest(here / name) == expected, f'Executed source changed: {name}'
    archived = json.loads((here / 'CHECKS.json').read_text())
    assert archived['status'] == 'PASS_FINITE_DIAGNOSTICS'
    assert archived['derived_universal_locality'] == 'NOT_ESTABLISHED'
    assert archived['new_formation_trajectories'] == 0
    findings = json.loads((here / 'FINDINGS.json').read_text())
    assert findings['independent_proof_review'] == 'PENDING'
    assert findings['new_locality_rule'] == 'SUPPLIED_NOT_DERIVED'
    print(f"PASS: {len(manifest['sha256'])} packet hashes, frozen charter, "
          f"{len(pins['files'])} upstream pins and executed source hashes")
    if args.run_checks:
        archive_before = digest(here / 'CHECKS.json')
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        result = subprocess.run([sys.executable, str(here / 'checks.py')],
                                cwd=here, env=env, capture_output=True,
                                text=True, timeout=60, check=True)
        replay = json.loads(result.stdout)
        assert replay['status'] == 'PASS_FINITE_DIAGNOSTICS'
        assert digest(here / 'CHECKS.json') == archive_before
        print('REPLAY PASS: obstruction, diagonal witness, full physical kernel, '
              'common-space comparison, prefix and forbidden-product controls; archive preserved')
    print('COMMUTANT_ONLY_UNIVERSALITY: FAIL_IN_CLASS; '
          'RELATIONAL_GENERATOR_RULE: SUPPLIED; INDEPENDENT_REVIEW: PENDING')


if __name__ == '__main__':
    main()
