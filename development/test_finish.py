"""The applied repair cannot suppress a new red, bank delta or altered patch."""
from copy import deepcopy
import hashlib
from pathlib import Path
import subprocess
import sys
import unittest

from finish_checks import (REPAIR_PATHS, PASS_CHANGES, ROOT, repaired_expectation,
                           validate_repair_diff)
from integrity import adjudication_problems, baseline_deltas, locked_manifest


class AppliedRepairTests(unittest.TestCase):
    def setUp(self):
        self.base = locked_manifest(ROOT / 'development/expected_red_manifest.json')
        self.patch = (ROOT / 'development/review_cycle_02/repair.patch').read_bytes()
        self.pin_patch = (ROOT / 'development/finish/pin-fixture.patch').read_bytes()

    def test_only_five_anchored_outcomes_change(self):
        repaired = repaired_expectation(self.base, self.base['protected_inputs'])
        changed = {n for n in self.base['pytest_cases'] if repaired['pytest_cases'][n] != self.base['pytest_cases'][n]}
        self.assertEqual(changed, set(PASS_CHANGES))
        for field in ('declarations', 'open_passes', 'bank_inventory', 'tool_versions'):
            self.assertEqual(repaired[field], self.base[field])

    def test_altered_or_unanchored_failure_is_rejected(self):
        for outcome in ('PASS', 'SKIP', None):
            base = deepcopy(self.base)
            base['pytest_cases'][PASS_CHANGES[0]] = outcome
            with self.subTest(outcome=outcome), self.assertRaises(ValueError):
                repaired_expectation(base, base['protected_inputs'])

    def test_exact_patches_and_only_their_six_protected_paths(self):
        validate_repair_diff(list(REPAIR_PATHS), self.patch, [], self.pin_patch)
        for paths, patch, unknown, pin_patch in (
            (list(REPAIR_PATHS) + ['claims.json'], self.patch, [], self.pin_patch),
            (list(REPAIR_PATHS[:-1]), self.patch, [], self.pin_patch),
            (list(REPAIR_PATHS), self.patch + b'\n', [], self.pin_patch),
            (list(REPAIR_PATHS), self.patch, ['new-scientific-input.txt'], self.pin_patch),
            (list(REPAIR_PATHS), self.patch, [], self.pin_patch + b'\n'),
        ):
            with self.subTest(paths=paths, unknown=unknown), self.assertRaises(ValueError):
                validate_repair_diff(paths, patch, unknown, pin_patch)

    def test_pin_probe_fires_without_writing_to_the_live_registry(self):
        paths = [ROOT / 'provenance/claims.json', ROOT / 'provenance/doc_register_pins.json']
        before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        backup = ROOT / 'provenance/claims.json.pintest'
        self.assertFalse(backup.exists())
        result = subprocess.run([sys.executable, '-m', 'pytest', '-q', 'test_doc_register_pins.py'],
            cwd=ROOT / 'provenance', capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual({str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}, before)
        self.assertFalse(backup.exists())

    def test_new_failure_or_silent_declared_red_pass_stays_a_delta(self):
        m = repaired_expectation(self.base, self.base['protected_inputs'])
        state = {k: m[k] for k in ('declarations', 'open_passes', 'bank_inventory', 'protected_inputs', 'tool_versions')}
        passing = next(n for n, v in m['pytest_cases'].items() if v == 'PASS' and n not in PASS_CHANGES)
        declared = next(iter(m['declarations']))
        for node, outcome in ((passing, 'FAIL'), (declared, 'PASS'), (PASS_CHANGES[0], 'FAIL')):
            cases = dict(m['pytest_cases']); cases[node] = outcome
            with self.subTest(node=node):
                self.assertIn({'kind': 'PYTEST_CASE', 'node': node,
                    'baseline': m['pytest_cases'][node], 'observed': outcome},
                    baseline_deltas(m, cases, state))
                kind = 'STALE_DECLARATION' if node == declared else 'UNDECLARED_FAILING_TEST'
                self.assertIn({'kind': kind, 'node': node}, adjudication_problems(cases, state))

    def test_bank_or_adjudication_refresh_is_rejected(self):
        m = repaired_expectation(self.base, self.base['protected_inputs'])
        state = {k: deepcopy(m[k]) for k in ('declarations', 'open_passes', 'bank_inventory', 'protected_inputs', 'tool_versions')}
        state['bank_inventory']['report']['deletions'] = []
        self.assertIn('BANK_INVENTORY', {d['kind'] for d in baseline_deltas(m, m['pytest_cases'], state)})
        state['open_passes']['P2-TERMINATION-EVENTLOG']['status'] = 'CLOSED'
        self.assertIn('OPEN_PASSES', {d['kind'] for d in baseline_deltas(m, m['pytest_cases'], state)})


if __name__ == '__main__':
    unittest.main()
