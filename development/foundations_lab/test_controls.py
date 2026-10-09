"""Checks of exact invariants and hostile witnesses; no scientific bank gate."""
import math
import unittest
from fractions import Fraction as F
from itertools import product
from controls import (run_controls, eye, mm, add, scale, tr, pr_probability)


class StandardControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = run_controls()

    def test_same_Z_readout_does_not_fix_state(self):
        c = self.r["C01_qubit"]
        self.assertEqual(c["P_Z_plus"], c["P_Z_minus"])
        self.assertEqual((c["P_X_plus_plus"], c["P_X_plus_minus"]), (1, 0))
        self.assertEqual(c["QFI_theta"], 4*c["FS_theta"])

    def test_angular_and_SI_witnesses(self):
        c = self.r["C02_units"]
        self.assertGreater(c["early_cos_formula"], 1)
        for a, b in [("v_formula_units_LT", "velocity_units_LT"),
                     ("Phi_formula_units_LT", "potential_units_LT"),
                     ("theta_second_derivative_units_LT", "acceleration_units_LT")]:
            self.assertNotEqual(c[a], c[b])

    def test_chiral_and_determinant_positive(self):
        c = self.r["C03_spectral"]
        self.assertEqual(c["anticommutator"], scale(0, eye(4)))
        self.assertEqual(c["trace_H"], 0)
        for row in c["determinant_samples"]:
            self.assertGreater(row["paper_ratio"], 0)
            self.assertEqual(row["ordinary_det_ratio"], row["paper_ratio"]**2)
        # General positivity proof is in PROOFS.md; this is only a fixture.

    def test_self_adjoint_does_not_imply_finite_Gibbs_trace(self):
        c = self.r["C04_gibbs"]
        for row in c["cutoffs"]:
            self.assertGreaterEqual(row["Z_N"], row["lower_bound"])
        self.assertGreater(c["cutoffs"][-1]["Z_N"], c["cutoffs"][0]["Z_N"])

    def test_encoding_dimension_obstruction(self):
        c = self.r["C05_encoding"]
        self.assertEqual(c["U_adjoint_U"], [[1, 0, 0], [0, 1, 0], [0, 0, 0]])
        self.assertEqual(c["identity_error_squared"], 1)

    def test_metric_nonuniqueness_witness(self):
        c = self.r["C06_metrics"]
        self.assertGreater(c["BKM_metric_float"], c["SLD_metric"])

    def test_factorization_witness_and_unitarity(self):
        c = self.r["C07_factorization"]
        self.assertEqual(mm(c["unitary"], c["unitary"]), eye(4))
        self.assertEqual((c["purity_before"], c["purity_after"]), (F(1, 2), 1))
        self.assertEqual(tr(c["rho"]), 1)
        self.assertEqual(tr(c["mapped_rho"]), 1)

    def test_PR_normalization_positivity_and_no_signalling(self):
        for x, y in product(range(2), repeat=2):
            self.assertEqual(sum(pr_probability(a, b, x, y) for a, b in product(range(2), repeat=2)), 1)
            for a in range(2):
                self.assertEqual(sum(pr_probability(a, b, x, y) for b in range(2)), F(1, 2))
            for b in range(2):
                self.assertEqual(sum(pr_probability(a, b, x, y) for a in range(2)), F(1, 2))
        self.assertTrue(all(row["p"] >= 0 for row in self.r["C08_PR_box"]["cells"]))
        self.assertEqual(self.r["C08_PR_box"]["PR_CHSH"], 4)
        self.assertEqual(self.r["C08_PR_box"]["local_max"], 2)

    def test_real_density_positivity_and_tomographic_failure(self):
        c = self.r["C09_local_tomography"]
        self.assertEqual(mm(c["YY"], c["YY"]), eye(4))
        for key in ("rho_plus", "rho_minus"):
            rho = c[key]
            self.assertEqual(tr(rho), 1)
            self.assertEqual(mm(rho, rho), scale(F(1, 2), rho))
            # Symmetric rho with eigenvalues 0 or 1/2 is positive semidefinite.
        self.assertTrue(all(row["plus"] == row["minus"] for row in c["local_expectations"]))
        self.assertEqual((c["global_plus"], c["global_minus"]), (1, -1))

    def test_reversible_kinetic_comparator(self):
        c = self.r["C10_kinetics"]
        self.assertAlmostEqual(sum(c["stationary"]), 1)
        for row in c["models"]:
            self.assertLess(abs(row["detailed_balance_residual"]), 1e-14)
            self.assertTrue(0 <= row["P_01_at_t1"] <= 1)
        self.assertGreater(c["models"][1]["P_01_at_t1"], c["models"][0]["P_01_at_t1"])

    def test_anomalies_cancel_with_free_scale(self):
        c = self.r["C11_anomalies"]
        for row in c["families"]:
            for key in ("SU3_squared_U1", "SU2_squared_U1", "gravity_squared_U1", "U1_cubed"):
                self.assertEqual(row[key], 0)
        self.assertEqual(c["families"][1]["u"], F(2, 3))
        self.assertNotEqual(c["families"][1]["q"], c["families"][2]["q"])

    def test_prime_arithmetic_and_noise_reversal(self):
        c = self.r["C12_prime_noise"]
        self.assertEqual(c["first_four_primes_product"], 210)
        self.assertEqual(c["units_mod_30"], [1, 7, 11, 13, 17, 19, 23, 29])
        self.assertEqual(c["nonunit_fraction"], F(11, 15))
        self.assertEqual(c["gcd_49_30"], 1)
        self.assertEqual((c["purity_plus_after"], c["purity_zero_after"]), (1, F(1, 2)))

    def test_May_velocity_definitions_inconsistent(self):
        c = self.r["C13_velocity"]
        self.assertEqual(c["vector_norm_squared_exact"], 1)
        self.assertNotEqual(c["vector_norm_squared_exact"], c["claimed_v_squared"])
        self.assertTrue(math.isclose(c["vector_norm_squared_float"], 1, abs_tol=1e-14))


if __name__ == "__main__":
    unittest.main()
