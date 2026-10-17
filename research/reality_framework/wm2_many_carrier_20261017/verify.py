"""Integrity and optional replay; no formation solver or archive replacement."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-checks', action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    repo = here.parents[2]
    manifest = json.loads((here / 'MANIFEST.json').read_text())
    for path, expected in manifest['sha256'].items():
        assert sha(here / path) == expected, f'Packet changed: {path}'
    pins = json.loads((here / 'SOURCE_PINS.json').read_text())
    assert sha(here / 'CHARTER.md') == pins['charter_sha256']
    for path, expected in pins['files'].items():
        assert sha(repo / path) == expected, f'Upstream changed: {path}'
    executed = json.loads((here / 'EXECUTED_SOURCES.json').read_text())
    for path, expected in executed['sha256'].items():
        assert sha(here / path) == expected, f'Executed code changed: {path}'
    checks = json.loads((here / 'CHECKS.json').read_text())
    predictions = json.loads((here / 'PREDICTIONS.json').read_text())
    findings = json.loads((here / 'FINDINGS.json').read_text())
    assert checks['status'] == 'PASS_FINITE_DIAGNOSTICS'
    assert checks['formation_trajectories'] == 0
    assert predictions['experimental_observations'] == 0
    assert predictions['target_time_tau'] == 32
    assert predictions['future_test_rule']['Hamiltonian_error_coefficient_with_time_calibration'] == 66
    assert findings['independent_mathematical_review'] == 'PENDING'
    assert findings['locality_status'] == 'SUPPLIED_PRIMITIVE'
    assert findings['neutral_seed_populated_formation'] == 'EXACTLY_ZERO_BY_NUMBER_CONSERVATION'
    print(f"PASS: {len(manifest['sha256'])} packet hashes, eight upstream pins, "
          "charter and executed-code hashes; pending-review/zero-data boundaries retained")
    if args.run_checks:
        before = {name: sha(here / name) for name in ['CHECKS.json', 'PREDICTIONS.json']}
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        for script in ['checks.py', 'predictions.py']:
            process = subprocess.run([sys.executable, str(here / script)], cwd=here,
                                     env=env, capture_output=True, text=True,
                                     timeout=60, check=True)
            replay = json.loads(process.stdout)
            if script == 'checks.py':
                assert replay['status'] == 'PASS_FINITE_DIAGNOSTICS'
            else:
                assert replay['status'] == 'CALCULATED_NOT_MEASURED'
                for species in ['A', 'B']:
                    for new, old in zip(replay['species'][species]['probabilities'],
                                        predictions['species'][species]['probabilities']):
                        for key in ['probability_lower', 'probability_upper', 'rounded_probability']:
                            assert new[key] == old[key]
            print(f'REPLAY PASS: {script}; archive preserved')
        assert all(sha(here / name) == expected for name, expected in before.items())
    print('FINITE_DENSITY_USE: AUTHOR_DERIVED_CONDITIONAL; '
          'AUTONOMOUS_POPULATED_PHASE: NOT_ESTABLISHED; EMPIRICAL_CONFIRMATION: NONE')


if __name__ == '__main__':
    main()
