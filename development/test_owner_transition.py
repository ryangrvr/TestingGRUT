"""Owner scoping cannot authorize any other failure, pass, outcome or bank delta."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from integrity import adjudication_problems, baseline_deltas, locked_manifest
from owner_transition import reconciled_manifest, ADDED, POINTER, P1A

ROOT = Path(__file__).resolve().parents[1]
TRANSITION = ROOT / 'development/final_precard/P6_TRANSITION.json'


class OwnerTransitionTests(unittest.TestCase):
    def setUp(self):
        self.base = locked_manifest(ROOT / 'development/expected_red_manifest.json')

    def test_only_exact_p2_cases_p1a_and_p6_dispositions_change(self):
        result = reconciled_manifest(self.base, TRANSITION, ROOT)
        expected = deepcopy(self.base)
        expected['declarations'][POINTER]['cases'].update(ADDED)
        expected['open_passes'][P1A]['status'] = 'CLOSED'
        expected['open_passes']['P6-STALE-NETS-IN-STANDING-DOCS']['status']='CLOSED'
        from owner_transition import P6_OUTCOMES
        expected['pytest_cases'].update(P6_OUTCOMES)
        expected['protected_inputs'] = result['protected_inputs']
        self.assertEqual(result, expected)
        self.assertEqual(result['open_passes']['P6-STALE-NETS-IN-STANDING-DOCS']['status'], 'CLOSED')

    def test_even_relocked_unapproved_dispositions_are_rejected(self):
        original = json.loads(TRANSITION.read_text())
        for attack in ('reopen-p6', 'close-p2', 'hide-kappa', 'extra-case', 'outcome-refresh', 'bank-refresh'):
            value = deepcopy(original)
            if attack=='reopen-p6' or attack=='close-p2':
                pid = 'P6-STALE-NETS-IN-STANDING-DOCS' if attack == 'reopen-p6' else 'P2-TERMINATION-EVENTLOG'
                value['open_passes'][pid]['status'] = 'OPEN' if attack=='reopen-p6' else 'CLOSED'
            elif attack == 'hide-kappa':
                key = next(k for k in value['declarations'][POINTER]['cases'] if 'KAPPA' in k)
                value['declarations'][POINTER]['cases'][key] = 'P2-TERMINATION-EVENTLOG'
            elif attack == 'extra-case':
                value['declarations'][POINTER]['cases']['extra'] = 'P2-TERMINATION-EVENTLOG'
            else:
                value['pytest_cases' if attack == 'outcome-refresh' else 'bank_inventory'] = {}
            with self.subTest(attack=attack), tempfile.TemporaryDirectory() as tmp:
                p = Path(tmp) / 'transition.json'
                p.write_text(json.dumps(value))
                p.with_suffix('.json.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest())
                with self.assertRaises(ValueError):
                    reconciled_manifest(self.base, p, ROOT)

    def test_closed_p6_cannot_hide_new_tier_red_or_bank_delta(self):
        result = reconciled_manifest(self.base, TRANSITION, ROOT)
        state = {k: result[k] for k in ('declarations', 'open_passes', 'bank_inventory', 'protected_inputs', 'tool_versions')}
        cases = dict(result['pytest_cases'])
        cases['test_resident.py::TestConsistency::test_tier_contradiction'] = 'FAIL'
        problems = adjudication_problems(cases, state)
        self.assertNotIn({'kind': 'ORPHANED_OPEN_PASS', 'pass': 'P6-STALE-NETS-IN-STANDING-DOCS'}, problems)
        self.assertTrue(any(p['kind'] == 'UNDECLARED_FAILING_TEST' and 'test_tier_contradiction' in p['node'] for p in problems))
        state = deepcopy(state)
        state['bank_inventory']['report']['deletions'] = []
        self.assertIn('BANK_INVENTORY', {d['kind'] for d in baseline_deltas(result, result['pytest_cases'], state)})


if __name__ == '__main__':
    unittest.main()
