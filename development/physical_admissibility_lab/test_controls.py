"""Hostile controls of stated embeddings and exclusions, not law tests."""
from fractions import Fraction as F
from itertools import product
import unittest
import controls as c


class PhysicalControls(unittest.TestCase):
    def test_original_equations_against_frozen_table(self):
        for z, a, b in product(range(8), range(2), range(2)):
            self.assertEqual(c.original_step(z, a, b), c.FIXED[z, a][b])

    def test_cptp_trace_preservation_and_choi_rank(self):
        result = c.cptp_control()
        for action in result['actions']:
            self.assertEqual(action['trace_preserving_diagonal'], [1]*8)
            self.assertEqual(action['choi_rank'], 16)
            self.assertEqual(action['choi_trace'], 8)
            self.assertGreater(min(action['nonzero_choi_eigenvalues']), 0)

    def test_all_feedback_records(self):
        r = c.cptp_control()
        self.assertEqual((r['comparisons'], r['record_cells'], r['maximum_difference']), (1024, 8192, 0))

    def test_instrument_reports_correct_bit_not_relabelled(self):
        for a in (0, 1):
            for _, zp, y, _ in c.channel_edges(a):
                self.assertEqual(y, zp//4)

    def test_register_probe_is_identity(self):
        r = c.memory_control()
        self.assertEqual(r['response_matrix'], c.eye(4))
        self.assertEqual(r['perfect_quantum_dimension_lower_bound'], 4)

    def test_generic_rank_bound_does_not_imply_perfect_discrimination(self):
        # Four qubit states have linearly independent operator coordinates,
        # while a common four-outcome POVM cannot decode equiprobably above 1/2.
        coordinates = [[1, 0, 0, 1], [1, 0, 0, -1], [1, 1, 0, 0], [1, 0, 1, 0]]
        self.assertEqual(c.rank(coordinates), 4)
        self.assertEqual(c.memory_control()['max_equiprobable_success_d1_d2_d3_d4'][1], F(1, 2))

    def test_reversible_controller_counts_garbage(self):
        r = c.reversible_control()
        self.assertEqual(r['domain_size'], r['distinct_images'])

    def test_reset_gram_rank_resource(self):
        for r in c.reset_control()['closed_unitary_bound']:
            self.assertEqual(r['orthogonal_input_gram_rank'], r['pure_environment_lower_bound'])

    def test_faithful_finite_bath_rank_obstruction(self):
        r = c.reset_control()['faithful_bath_example']
        self.assertGreater(r['initial_rank'], r['pure_reset_output_rank_at_most'])

    def test_landauer_identity_and_cost(self):
        r = c.landauer_control()
        self.assertLess(abs(r['identity_residual']), 1e-14)
        self.assertGreater(r['beta_Q'], r['delta_S'])
        self.assertGreater(r['relative_entropy'], 0)
        self.assertFalse(r['exact_pure_reset'])

    def test_Gibbs_preserving_does_not_imply_DB(self):
        r = c.balance_control()
        self.assertEqual(r['stationary_distribution'], r['stationary_output'])
        self.assertEqual(r['pair_current_matrix'][0][1], F(1, 6))

    def test_energy_conserving_thermal_operation_same_visible_kernel(self):
        r = c.balance_control()
        self.assertTrue(r['thermal_operation_permutation_unitary'])
        self.assertTrue(r['thermal_operation_exact_same_kernel'])
        self.assertFalse(r['autonomous_detailed_balance_semigroup'])

    def test_DB_hidden_lift_pair_symmetry(self):
        # Three-state nonuniform reversible Q, a stochastic discrete-time DB P.
        # P=I+Q/10 is DB and stochastic; arbitrary deterministic coarse readout.
        r = c.kinetic_control()
        p = c.add(c.eye(3), c.scale(r['Q'], F(1, 10)))
        labels = [0, 1, 0]
        pair = [[sum(r['pi'][i]*p[i][j] for i in range(3) for j in range(3)
                     if labels[i] == a and labels[j] == b) for b in (0, 1)] for a in (0, 1)]
        self.assertEqual(pair[0][1], pair[1][0])

    def test_NOT_excluded_from_CTMC_but_unitary(self):
        r = c.ctmc_control()
        self.assertEqual(r['NOT_determinant'], -1)
        self.assertEqual(c.mm(c.X, c.X), c.eye(2))
        self.assertGreater(r['minimum_error_numerical'], 0)

    def test_positive_transpose_not_completely_positive(self):
        r = c.transpose_control()
        self.assertEqual(r['singlet_expectation'], F(-1, 2))
        self.assertEqual(sum(r['eigenvalues']), 1)

    def test_PR_normalized_no_signalling(self):
        r = c.bell_control()
        for p in r['PR_distributions'].values():
            self.assertEqual(sum(map(sum, p)), 1)
            self.assertEqual(list(map(sum, p)), [F(1, 2)]*2)
            self.assertEqual([sum(row[j] for row in p) for j in (0, 1)], [F(1, 2)]*2)
        self.assertGreater(r['PR_CHSH_squared'], r['quantum_CHSH_bound_squared'])
        self.assertEqual(r['local_CHSH_max'], 2)

    def test_local_quantum_Hamiltonian_has_tails(self):
        self.assertTrue(all(c.locality_control()['exact_Pauli_identities'].values()))

    def test_classical_shift_strict_cone_and_arrival(self):
        r = c.locality_control()['classical_local_shift_intervention']
        self.assertEqual([row['site_5_difference'] for row in r[:5]], [0]*5)
        self.assertEqual(r[5]['site_5_difference'], 1)

    def test_reversible_conductance_constraints(self):
        r = c.kinetic_control()
        self.assertEqual(r['stationary_generator_product'], [0]*3)
        for i, j in product(range(3), repeat=2):
            self.assertEqual(r['pi'][i]*r['Q'][i][j], r['pi'][j]*r['Q'][j][i])
        self.assertEqual(list(map(sum, r['Q'])), [0]*3)

    def test_GKLS_rates_and_Gibbs_ratio(self):
        r = c.kinetic_control()
        self.assertEqual(r['GKLS_Gibbs_population'][1]/r['GKLS_Gibbs_population'][0], F(1, 3))
        self.assertEqual(r['GKLS_pop_rates_gamma1_gamma2'][1], 2*r['GKLS_pop_rates_gamma1_gamma2'][0])
        self.assertEqual(r['GKLS_stationary_density_derivative'], [[0, 0], [0, 0]])
        self.assertEqual(r['GKLS_coherence_rates_gamma1_phi0_phi1'], [F(2, 3), F(5, 3)])

    def test_anomaly_cancellation_does_not_select_charge_scale(self):
        r = c.anomaly_control()
        self.assertNotEqual(r['1,1']['sum_q3'], 0)
        self.assertEqual(r['1,-1'], r['2,-2'])


if __name__ == '__main__':
    unittest.main()
