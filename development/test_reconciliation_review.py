"""Preconditions and contamination controls for review-only patches."""
from pathlib import Path
import unittest
from reconciliation_review import isolated_fixture_source, replace_once


class ReviewPatchTests(unittest.TestCase):
    def test_source_drift_or_ambiguous_replacement_cannot_be_patched(self):
        for body in ('absent', 'old old'):
            with self.assertRaises(ValueError):
                replace_once(body, 'old', 'new')

    def test_the_synthetic_fixture_restores_live_globals(self):
        root = Path(__file__).resolve().parents[1]
        src = (root / 'provenance/test_expected_red.py').read_text()
        patched = isolated_fixture_source(src)
        import sys
        from copy import deepcopy
        from unittest.mock import patch
        sys.path.insert(0, str(root / 'provenance'))
        import expected_red as live
        namespace = {'__name__': 'review_fixture'}
        exec(compile(patched, 'isolated_fixture.py', 'exec'), namespace)
        declared = live.DECLARED
        passes = live.open_passes
        inventory = deepcopy({n: e['cases'] for n, e in declared.items()})
        hostile_passes = lambda: {'P9-CONTAMINATION': {'status': 'OPEN', 'symptomless': False}}
        with patch.object(live, 'open_passes', hostile_passes):
            case = namespace['TestExpectedRed']('test_the_declared_state_is_accepted')
            case.setUp()
            try:
                rc, _ = case._run()
                self.assertEqual(rc, 0)
                self.assertNotIn('P9-CONTAMINATION', live.open_passes())
                self.assertIsNot(live.DECLARED, declared)
            finally:
                case.tearDown()
            self.assertIs(live.open_passes, hostile_passes)
        self.assertIs(live.DECLARED, declared)
        self.assertIs(live.open_passes, passes)
        self.assertEqual({n: e['cases'] for n, e in declared.items()}, inventory)


if __name__ == '__main__':
    unittest.main()
