"""Reproduce exact controls and fail closed on stale lane evidence."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
import unittest
from controls import run_controls

ROOT = Path(__file__).resolve().parent
HASHED = ('controls.py', 'test_controls.py', 'run.py', 'verify_fixtures.py',
          'PROOFS.md', 'AUDIT.md', 'EXPERIMENTS.md', 'REPORT.md', 'README.md',
          'SOURCES.json', 'INPUTS.json')


def canonical(value):
    if isinstance(value, Fraction):
        return {'rational': [value.numerator, value.denominator]}
    if isinstance(value, float):
        return round(value, 12)
    if isinstance(value, dict):
        return {key: canonical(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [canonical(v) for v in value]
    return value


def table(evidence):
    c = evidence['controls']
    lines = ['# Generated physical-comparator control results', '',
             'Standard physics only. DRAFT, same-author validation; no candidate or Card.', '',
             '| Control | Result | Boundary |', '|---|---|---|',
             f"| Feedback record comparisons / cells | {c['C01_fixed_CPTP_lift']['comparisons']} / {c['C01_fixed_CPTP_lift']['record_cells']} | All depth-three deterministic binary feedback policies and eight initial states |",
             f"| Exact maximum record difference | {c['C01_fixed_CPTP_lift']['maximum_difference']} | Classical equations versus fixed CP instrument |",
             '| Choi rank, each setting | 16 | Minimum environment for the chosen dephasing quantum extension only |',
             '| Register response | I4 | Perfect quantum retained dimension at least 4 |',
             '| Decode success bound, d=1/2/3/4 | 1/4, 1/2, 3/4, 1 | All retained preparation information counted |',
             f"| Reversible control inputs / distinct outputs | {c['C03_reversible_finite_control']['domain_size']} / {c['C03_reversible_finite_control']['distinct_images']} | Inputs/noise/scratch retained |",
             '| N=3 qubit reset environment dimension | at least 8 | Closed unitary, independent inputs, exact reset |',
             '| Finite faithful bath example ranks | 6 initial versus at most 3 after pure reset | Unitary rank obstruction |',
             f"| Finite SWAP beta Q / delta S / D | {c['C05_finite_bath_Landauer']['beta_Q']:.12f} / {c['C05_finite_bath_Landauer']['delta_S']:.12f} / {c['C05_finite_bath_Landauer']['relative_entropy']:.12f} | Numerical log illustration; external work priced |",
             '| Qutrit oriented stationary pair current | 1/6 | DB lift excluded; energy-conserving thermal stroke passes |',
             '| NOT determinant | -1 | Finite-time CTMC excluded; Hamiltonian qubit passes |',
             '| Partial transpose singlet expectation | -1/2 | Positivity alone fails CP |',
             '| PR / quantum CHSH squared | 16 / bound 8 | No-signalling larger than quantum set |',
             '| Three-spin commutator norm | 2 sin(sqrt(2) J t)^2 | Local Hamiltonian has early tails |',
             '| Fixed-pi conductance freedom, N=2/3/4 | 1 / 3 / 6 dimensions | Standard reversible kinetic parameters |',
             '| Thermal GKLS population rates, gamma=1/2 | 4/3 / 8/3 | Same Gibbs ratio 1/3 |',
             '| Thermal GKLS coherence rates, phi=0/1 | 2/3 / 5/3 | Standard dephasing freedom |',
             '| Charge anomaly sums | (1,1) fails; (1,-1),(2,-2) pass | No charge normalization selector |', '',
             'T1–T10 analytical proofs are in PROOFS.md. Finite fixtures do not prove universal claims.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.discover(str(ROOT), pattern='test_controls.py'))
    if not tests.wasSuccessful():
        return 1
    evidence = {'schema': 'PHYSICAL_ADMISSIBILITY_STANDARD_CONTROLS_V1',
                'scope': 'ISSUE_7_NO_K_NO_CANDIDATE_NO_CARD',
                'review': 'DRAFT_PENDING_EXTERNAL_REVIEW_SAME_AUTHOR_CHECKS',
                'python': platform.python_version(), 'tests': tests.testsRun,
                'numerical_serialization': '12 decimal places; tests retain tolerance on raw calculations',
                'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in HASHED},
                'controls': run_controls()}
    encoded = json.dumps(canonical(evidence), indent=2, sort_keys=True, allow_nan=False)+'\n'
    generated = table(evidence)
    results, tables = ROOT/'RESULTS.json', ROOT/'TABLES.md'
    if args.check:
        if (not results.is_file() or not tables.is_file() or
                results.read_text() != encoded or tables.read_text() != generated):
            print('FAIL: missing/stale control evidence, source/code hash or table', file=sys.stderr)
            return 1
        print(f'{tests.testsRun} tests PASS; 12 standard-control groups; evidence matches. DRAFT; Card 1 CLOSED.')
    else:
        results.write_text(encoded)
        tables.write_text(generated)
        print('Generated this lane evidence only; no scientific approval.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
