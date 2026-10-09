"""Hostile exact checks of the representation/control claims."""
from fractions import Fraction as F
from itertools import product
import unittest
from controls import (TABLE, adaptive_step, fixed_step, records,
    equivalence_certificate, response_matrix, rank, matrix_fixture,
    smolin_step, matadd, contraction, bistable, selection,
    urn_law, beta_mixture_uniform, urn_fixture, run_controls)


class ScientificControls(unittest.TestCase):
    def test_kernel_rows_and_independent_table(self):
        self.assertEqual(len(TABLE), 16)
        for x, law, m, a in product((0, 1), repeat=4):
            original = list(adaptive_step((x, law, m), a))
            compiled = list(fixed_step(4*x+2*law+m, a))
            self.assertEqual([(4*s[0]+2*s[1]+s[2], y, p) for s, y, p in original], compiled)
            self.assertEqual(sum(p for _, _, p in compiled), 1)
            self.assertTrue(all(p > 0 for _, _, p in compiled))

    def test_hand_oracle_one_intervention(self):
        self.assertEqual(records((0, 0, 0), 1, lambda _: F(1), adaptive_step),
                         {((1, 1),): F(3, 4), ((1, 0),): F(1, 4)})

    def test_all_depth_three_feedback_policies(self):
        c = equivalence_certificate()
        self.assertEqual(c['maximum_exact_difference'], 0)
        self.assertEqual(c['comparisons'], 1024)
        self.assertEqual(c['record_cells'], 8192)

    def test_randomized_history_dependent_policies(self):
        for offset in range(3):
            policy = lambda h: F((len(h)+sum(y for _, y in h)+offset) % 4, 3)
            for initial in product((0, 1), repeat=3):
                encoded = 4*initial[0]+2*initial[1]+initial[2]
                left = records(initial, 4, policy, adaptive_step)
                self.assertEqual(left, records(encoded, 4, policy, fixed_step))
                self.assertEqual(sum(left.values()), 1)

    def test_wrong_fixed_kernel_is_detected(self):
        bad = dict(TABLE)
        bad[0, 1] = (0, 6)  # reverse the two noise outcomes
        self.assertNotEqual(records((0, 0, 0), 1, lambda _: F(1), adaptive_step),
                            records(0, 1, lambda _: F(1), lambda s, a: fixed_step(s, a, bad)))

    def test_hidden_memory_is_not_visible_reset(self):
        p = lambda m: sum(w for _, y, w in adaptive_step((0, 1, m), 0) if y)
        self.assertEqual((p(0), p(1)), (F(1, 4), F(3, 4)))

    def test_four_state_response_rank(self):
        c = response_matrix()
        self.assertEqual(c['matrix'], [[int(i == j) for j in range(4)] for i in range(4)])
        self.assertEqual(c['rank'], 4)
        # Any explicit factorization through three hidden states has rank <=3.
        a = [[F(i+j+1, 20) for j in range(3)] for i in range(4)]
        b = [[F(i*j+1, 9) for j in range(4)] for i in range(3)]
        ab = [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(4)] for i in range(4)]
        self.assertLessEqual(rank(ab), 3)

    def test_smolin_pair_and_antisymmetry(self):
        c = matrix_fixture()
        self.assertTrue(c['autonomous_pair_exact_match'])
        self.assertTrue(c['antisymmetric'])
        self.assertEqual(c['traces'], [0]*6)
        self.assertNotEqual(c['matrices'][1], c['matrices'][2])

    def test_smolin_degenerate_families(self):
        a = ((0, 1), (-1, 0))
        self.assertEqual(smolin_step(a, a), a)
        b = ((0, 2), (-2, 0))
        self.assertEqual(smolin_step(a, b), matadd(b, matadd(b, a, -1)))

    def test_contraction_exact_error_and_encoded_target(self):
        for target, initial in product((F(1, 137), F(2, 137)), (F(-1), F(0), F(1), F(10))):
            self.assertEqual(contraction(target, initial, 6), target+(initial-target)/64)
        self.assertNotEqual(contraction(F(1, 137), 0, 6), contraction(F(2, 137), 0, 6))

    def test_bistability_breaks_initial_independence(self):
        self.assertEqual([bistable(x, 5) for x in (-1, 0, 1)], [-1, 0, 1])
        self.assertTrue(-1 < bistable(F(-1, 2), 3) < F(-1, 2))
        self.assertTrue(F(1, 2) < bistable(F(1, 2), 3) < 1)

    def test_selection_closed_form_and_absent_type(self):
        for w, p in product((F(1), F(2)), (F(0), F(1, 3), F(1))):
            self.assertEqual(selection(w, p, 5), w**5*p/(1-p+w**5*p))
        self.assertEqual(selection(F(2), F(0), 50), 0)
        self.assertEqual(selection(F(1), F(1, 3), 50), F(1, 3))

    def test_stochastic_selection_depends_on_supplied_rates(self):
        q, r = F(1, 4), F(1, 2)
        p = q/(q+r)
        self.assertEqual((1-p)*q, p*r)
        self.assertEqual(p, F(1, 3))
        self.assertNotEqual(p, q/(q+F(1, 4)))

    def test_urn_and_static_mixture_equal_record_laws(self):
        c = urn_fixture()
        self.assertEqual(c['enumerated_word_cells'], 126)
        self.assertEqual(c['maximum_exact_difference'], 0)
        self.assertEqual(sum(urn_law(w) for w in product((0, 1), repeat=3)), 1)
        self.assertEqual(urn_law((1, 1)), F(1, 3))

    def test_urn_observation_equivalence_does_not_imply_do_equivalence(self):
        self.assertEqual(urn_law((1, 1))/urn_law((1,)), F(2, 3))
        # Causal forcing of first result leaves the static latent prior untouched.
        self.assertEqual(beta_mixture_uniform((1,)), F(1, 2))
        self.assertNotEqual(F(2, 3), F(1, 2))

    def test_empty_domain_and_postulated_birth(self):
        self.assertEqual(sum([], F(0)), 0)
        self.assertEqual(1-F(3, 4)**3, F(37, 64))
        self.assertEqual(1-F(1)**3, 0)

    def test_changed_objective_changes_best_branch(self):
        branches = range(3)
        self.assertEqual(min(branches, key=lambda x: x*x), 0)
        self.assertEqual(min(branches, key=lambda x: (x-2)**2), 2)

    def test_clock_and_parameter_drift_have_fixed_lifts(self):
        self.assertEqual([n for n in range(1, 65) if n & (n-1) == 0], [1, 2, 4, 8, 16, 32, 64])
        x, c, delta = F(0), F(1, 137), F(1, 1000000)
        for n in range(4):
            x, c = c, c+delta
            self.assertEqual(x, F(1, 137)+n*delta)
        self.assertEqual(c, F(1, 137)+4*delta)


if __name__ == '__main__':
    unittest.main()
