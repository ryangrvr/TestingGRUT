#!/usr/bin/env python3
"""verify.py -- integrity and headline-number checks for packet grut_closed_whole_theory_20261019.

Usage:  python3 code/verify.py            # hash check against MANIFEST.json + headline checks on results/*.json
        python3 code/verify.py --rerun    # additionally re-run both calculations into a temp dir and compare numerics
Exit status 0 = all checks pass.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

PKT = Path(__file__).resolve().parent.parent
fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# 1. manifest hashes
man = json.loads((PKT / "MANIFEST.json").read_text())
for rel, h in man["files"].items():
    p = PKT / rel
    check(p.exists() and sha(p) == h, f"hash {rel}")

# 2. headline numbers
d = json.loads((PKT / "results" / "closed_whole_demo.json").read_text())
m = json.loads((PKT / "results" / "marginal_compatibility.json").read_text())
S1, S2, S3, S4 = d["S1_preparation_evolution_measurement"], d["S2_nonmonotone_relaxation"], d["S3_derived_reduced_dynamics"], d["S4_initial_correlations"]
S5, S6, S7, S8 = d["S5_nonunique_local_law"], d["S6_records_and_interventions"], d["S7_time_reversal"], d["S8_what_closure_forces"]
inv = S1["invariants"]
check(inv["max_abs_dE_tot"] < 1e-12 and inv["max_abs_dS_tot"] < 1e-12 and inv["max_abs_dpurity_tot"] < 1e-12, "S1 global invariants conserved to 1e-12")
check(S1["born_global_minus_local_max"] < 1e-12, "S1 local Born = global Born")
check(S1["reduced_equation_check_t3"]["max_abs_lhs_minus_rhs_full"] < 1e-8 and S1["reduced_equation_check_t3"]["max_abs_Tr_Abar_commutator_H_Abar"] < 1e-12, "S1 exact reduced equation")
check(S2["number_of_sign_changes_of_dS_A_dt_on_grid"] >= 10 and S2["BLP"]["backflow_sum_of_positive_increments"] > 0.5, "S2 non-monotone S_A and information backflow")
check(S3["cptp_check"]["min_choi_eigenvalue_over_samples"] > -1e-9 and S3["cp_divisibility"]["most_negative"] < -0.1, "S3 CPTP map, not CP-divisible")
check(max(S3["tcl"]["generator_bath_state_dependence_fro_norm_K_thermal-K_ground"].values()) > 0.05, "S3 generator depends on bath state")
q = S4["inf_T_bath_quantum_correlation_Y0X1"]
check(q["max_marginal_deviation_between_rho1_and_rho2"] < 1e-12 and q["min_eigenvalue_rho2"] >= -1e-13 and q["max_trace_distance_rhoA_over_grid"] > 0.05, "S4 identical marginals, different futures")
check(abs(q["norm_of_inhomogeneous_term_at_t0"] - math.sqrt(2) / 4) < 1e-9, "S4 inhomogeneous term = (1/4) Z_0 (norm sqrt(2)/4)")
check(S5["model_b_prime_exact_recurrence"]["trace_distance_rho_global_T_vs_0"] < 1e-10 and S5["model_b_prime_exact_recurrence"]["S_A_at_period"] < 1e-10, "S5 exact recurrence at T = 4 pi")
check(S5["trace_distance_rhoA"]["5.0"]["a_vs_b"] > 0.3, "S5 different complements, different local dynamics")
check(S6["abs_difference"] < 1e-12 and abs(S6["pZ_plus_A_at_t2_with_record_written"] - S6["pZ_plus_A_at_t2_with_measure_and_forget_dephasing"]) < 1e-12, "S6 record = Born probability; unitary record = instrument")
check(abs(S6["pZ_plus_A_at_t2_with_record_written"] - S6["pZ_plus_A_at_t2_without_intervention"]) > 0.05, "S6 intervention back-action")
rp = S6["record_persistence"]["10.0"]
check(rp["idle_register_trace_distance_vs_written"] < 1e-12 and rp["kappa_X_R_optimal_readout_success"] > 1 - 1e-9 and rp["kappa_X_R_X_7_optimal_readout_success"] < 1 - 1e-3, "S6 record persistence depends on H (idle exact; X_R rotation lossless; bath-coupled register loses retrievability)")
check(abs(S6["pZ_plus_A_at_t2_composed_one_time_maps_bath_reprepared"] - S6["pZ_plus_A_at_t2_with_record_written"]) > 0.02, "S6 composed one-time maps differ from the exact closed-whole result")
check(S7["trace_distance_after_forward_then_backward_evolution"] < 1e-12 and S7["entropy_of_A_along_the_conjugated_reversed_preparation"]["40.0"] < 1e-9 and S7["max_abs_rho_rev(t)_minus_rho(T-t)_conj"] < 1e-12, "S7 reversal, rho_rev(t) = rho(T-t)^*, entropy returns to 0")
mx = S8["mixed_S1_trajectory"]
check(S8["pure_global_state_S_A_equals_S_Abar_max_dev"] < 1e-12 and mx["subadditivity_violations"] == 0 and mx["araki_lieb_violations"] == 0 and mx["araki_lieb_slack_at_t3_t20"][0] > 0.1, "S8 entropy theorems (pure: S_A = S_Abar; mixed: subadditivity and Araki-Lieb strict)")
check(S2["first_local_max"]["t"] < 5 and S2["first_local_min_after_it"]["S_A"] < S2["first_local_max"]["S_A"] - 0.05, "S2 early local maximum followed by a drop")
check(m["C2_bell_chain"]["sdp_status_common_global_state"] == "infeasible", "C2 Bell chain infeasible")
check(abs(m["C3_quantitative_obstruction"]["max_F_AB_plus_F_BC_over_all_3qubit_states"] - 1.5) < 1e-5, "C3 max F_AB+F_BC = 3/2")
check(abs(m["C4_werner_thresholds"]["chain_AB_BC_max_p"] - 2 / 3) < 1e-5 and abs(m["C4_werner_thresholds"]["triangle_AB_BC_CA_max_p"] - 1 / 3) < 1e-5, "C4 Werner thresholds 2/3 and 1/3")
check(m["C5_classical_contrast"]["chain_extension_max_abs_marginal_error"] < 1e-12 and m["C5_classical_contrast"]["odd_cycle_anticorrelated_status"] == "infeasible", "C5 classical chain extends, odd cycle does not")
check(m["C6_inconsistent_overlap"]["sdp_status"] == "infeasible" and m["C6_inconsistent_overlap"]["sdp_status_second_solver_SCS"] == "infeasible", "C6 inconsistent overlap infeasible (both solvers)")
check(m["C2_bell_chain"]["sdp_status_second_solver_SCS"] == "infeasible" and m["C2_bell_chain"]["minimum_achievable_marginal_violation_frobenius"] > 0.2, "C2 second solver and positive minimum violation")
check(m["C4_werner_thresholds"]["chain_min_violation_at_p_plus_1e-3"] > 1e-5 and m["C4_werner_thresholds"]["triangle_min_violation_at_p_plus_1e-3"] > 1e-5, "C4 violation certificates above thresholds")
check(m["C5_classical_contrast"]["even_cycle_four_anticorrelated_status"] == "optimal", "C5 even anticorrelated cycle extends")
check(d["provenance"]["script_sha256"] == sha(PKT / "code" / "closed_whole_demo.py"), "provenance hash closed_whole_demo.py matches results")
check(m["provenance"]["script_sha256"] == sha(PKT / "code" / "marginal_compatibility.py"), "provenance hash marginal_compatibility.py matches results")

# 3. optional re-run
if "--rerun" in sys.argv:
    def leaves(o, pre=""):
        if isinstance(o, dict):
            for k, v in o.items():
                yield from leaves(v, pre + "/" + k)
        elif isinstance(o, list):
            for i, v in enumerate(o):
                yield from leaves(v, pre + f"[{i}]")
        else:
            yield pre, o
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "code").mkdir(); (td / "results").mkdir()
        for s in ("closed_whole_demo.py", "marginal_compatibility.py"):
            (td / "code" / s).write_bytes((PKT / "code" / s).read_bytes())
            r = subprocess.run([sys.executable, str(td / "code" / s)], capture_output=True, text=True)
            check(r.returncode == 0, f"re-run {s}")
        for name, ref in (("closed_whole_demo.json", d), ("marginal_compatibility.json", m)):
            new = json.loads((td / "results" / name).read_text())
            worst = 0.0
            for (k1, v1), (k2, v2) in zip(leaves(ref), leaves(new)):
                if k1.startswith("/provenance") or k1 != k2:
                    continue
                if isinstance(v1, (int, float)) and isinstance(v2, (int, float)) and not isinstance(v1, bool):
                    worst = max(worst, abs(v1 - v2) / max(1e-12, abs(v1)) if abs(v1) > 1e-9 else abs(v1 - v2))
                elif v1 != v2:
                    worst = max(worst, 1.0)
            check(worst < 1e-6, f"re-run {name} agrees with committed results (worst rel dev {worst:.2e})")

print("\n%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
