"""Abstract arithmetic controls; no private bundle needed in ordinary CI."""
from pathlib import Path
import tempfile
import unittest
from cr5_review import oracle as O
from cr5_review.check import unique_json, verified_bytes


class CR5ReviewTests(unittest.TestCase):
    def test_reference_errors_add_probe_annulus_and_tail_without_fitting(self):
        result = O.reference({'x': (3, '1/16')}, {'x': O.F(1, 8) + O.F(1, 4)})
        self.assertEqual(result, {'x': (O.F(3), O.F(7, 16))})

    def test_the_whole_chart_is_needed_for_exclusion(self):
        current = {'x': (0, 0), 'y': (2, 0)}
        self.assertEqual(O.decide(current, {'x': (0, 0), 'y': (0, 0)}, {'x': 1, 'y': 1}), 'INCLUDE')
        with self.assertRaises(ValueError):
            O.decide(current, {'x': (0, 0)}, {'x': 1})

    def test_invalid_exact_bounds_and_measures_are_rejected(self):
        for val in (True, 0.1, None):
            with self.assertRaises(ValueError):
                O.exact(val)
        for weights in ([-1, 2], [0, 0]):
            with self.assertRaises(ValueError):
                O.measure(weights, [False, False], [False, True])
        with self.assertRaises(ValueError):
            O.reference({'x': (0, 0)}, {'x': -1})

    def test_distinct_fallbacks_and_undefined_ablation_denominator(self):
        self.assertEqual(O.controls(), 7)
        self.assertIsNone(O.measure([1, 1], [True, True], [True, True]))
        self.assertEqual(O.nr17(['EXCLUDE', 'EXCLUDE'], [True, True]), [False, False])
        self.assertEqual(O.measure([1, 1], [False, False], [False, True]), O.F(1, 2))

    def test_a_missing_chain_never_borrows_the_other_certificate(self):
        self.assertEqual(O.q3('EXCLUDE', 'UNRESOLVED'), 'UNRESOLVED')
        self.assertEqual(O.q3('INCLUDE', 'UNRESOLVED'), 'UNRESOLVED')

    def test_transport_uses_absolute_matrix_for_error(self):
        self.assertEqual(O.transport([[1, -1]], {'x': (2, '1/4'), 'y': (-1, '1/8')}),
                         [(O.F(3), O.F(3, 8))])

    def test_tampered_bundle_is_refused_before_any_code_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'tampered.zip'
            path.write_bytes(b'not the reviewed bundle')
            with self.assertRaisesRegex(ValueError, 'no code imported'):
                verified_bytes(path)

    def test_duplicate_manifest_keys_cannot_hide_a_changed_hash(self):
        with self.assertRaises(ValueError):
            unique_json('{"file":"old","file":"new"}')


if __name__ == '__main__':
    unittest.main()
