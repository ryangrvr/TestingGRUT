"""Run the relic_yield module end to end and write results/headline.json with all assertions.

Usage:  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py      (about 4 minutes, deterministic)
"""
import hashlib
import json
import math
import os
import runpy
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(HERE)
from common import RESULTS, OMEGA_DM_KNOWN, M_H, M_X, M_D, T_UNIVERSE_S, dump  # noqa: E402

t0 = time.time()
timing = {}
for script in ("r1_charge.py", "r2_spectrum.py", "r3_rates.py", "r4_boltzmann.py", "r5_charge_history.py"):
    t = time.time()
    runpy.run_path(os.path.join(HERE, script), run_name="__main__")
    timing[script] = time.time() - t


def load(name):
    return json.load(open(os.path.join(RESULTS, name)))


C = load("charge_conservation.json")
S = load("spectrum_stability.json")
R = load("rates.json")
B = load("boltzmann.json")
H = load("charge_history.json")

asserts = []


def check(name, ok, value):
    asserts.append(dict(name=name, passed=bool(ok), value=value))


rc = R["checks"]
check("vertex tensors analytic vs finite difference (abs < 1e-9)", max(rc["vertex_tensor_maxdiff_cubic"], rc["vertex_tensor_maxdiff_quartic"]) < 1e-9,
      [rc["vertex_tensor_maxdiff_cubic"], rc["vertex_tensor_maxdiff_quartic"]])
check("hh->XXbar recursion vs hand-written tree amplitude (rel < 1e-12)", rc["hh_to_XXb_vs_analytic_maxrel"] < 1e-12, rc["hh_to_XXb_vs_analytic_maxrel"])
check("charge-violating amplitude vanishes", rc["charge_violating_amplitude_XX_to_hh"] == 0.0, rc["charge_violating_amplitude_XX_to_hh"])
check("Z2_D-violating amplitude vanishes", rc["Z2D_violating_amplitude_hD_to_XXb"] == 0.0, rc["Z2D_violating_amplitude_hD_to_XXb"])
check("Lorentz invariance of 3->2 amplitude (rel < 1e-12)", rc["lorentz_invariance_rel"] < 1e-12, rc["lorentz_invariance_rel"])
check("identical-particle symmetry of 3->2 amplitude (rel < 1e-12)", rc["identical_particle_swap_rel"] < 1e-12, rc["identical_particle_swap_rel"])
check("MC thermal <sigma v^2> -> threshold value at x=3000 (within 2%)", abs(rc["mc_over_threshold_XXXb_to_Xh_x3000"] - 1) < 0.02, rc["mc_over_threshold_XXXb_to_Xh_x3000"])
check("Gondolo-Gelmini vs direct MC for hh->XXbar at x=5 (within 2%)", abs(rc["GG_over_MC_hh_to_XXb_x5"] - 1) < 0.02, rc["GG_over_MC_hh_to_XXb_x5"])
check("pure-h 3->2 threshold amplitude nullified (|M|^2 < 1e-30, known broken-phi^4 property)", rc["hhh_to_hh_threshold_M2"] < 1e-30, rc["hhh_to_hh_threshold_M2"])
check("Noether identity q1 dU/dq2 - q2 dU/dq1 = 0 (rel < 1e-12)", C["noether_identity_max_rel"] < 1e-12, C["noether_identity_max_rel"])
L = C["lattice_FRW_check"]
check("lattice FRW comoving charge conserved (rel < 1e-12)", L["max_rel_change_comoving_Q"] < 1e-12, L["max_rel_change_comoving_Q"])
check("lattice FRW physical charge density ~ a^-3 (rel < 1e-10)", abs(L["physical_density_ratio_final_over_initial"] / L["expected_a_minus3"] - 1) < 1e-10,
      L["physical_density_ratio_final_over_initial"] / L["expected_a_minus3"])
check("lattice fields finite", L["all_finite"], L["max_abs_field_final_GeV"])
E = B["entropy_check_collisionless"]
check("collisionless Boltzmann run conserves hidden entropy and numbers (rel < 1e-7)", max([E["rel_change_entropy"]] + E["rel_change_numbers"]) < 1e-7,
      [E["rel_change_entropy"]] + E["rel_change_numbers"])
low = B["scan_all_species"][0]
check("xi=1e-5 run reproduces frozen initial numbers (rel < 2e-3)", abs(low["Omega_h2"] / low["Omega_frozen_initial_numbers"] - 1) < 2e-3,
      low["Omega_h2"] / low["Omega_frozen_initial_numbers"])
xm = B["xi_max"]
check("Omega_sym(xi_max) = 0.12 (rel < 1e-3)", abs(xm["run_at_xi_max"]["Omega_h2"] / OMEGA_DM_KNOWN - 1) < 1e-3, xm["run_at_xi_max"]["Omega_h2"])
check("3->2 out of equilibrium at xi_max (max Gamma32/H < 1)", xm["run_at_xi_max"]["max_G32_over_H"] < 1, xm["run_at_xi_max"]["max_G32_over_H"])
check("Boltzmann solver status ok for all scan points", all(r["status"] in (0, 1) for r in B["scan_all_species"]), [r["status"] for r in B["scan_all_species"]])
check("h -> X Xbar and h -> D D kinematically closed", M_H < 2 * M_X and M_H < 2 * M_D, [M_H, 2 * M_X, 2 * M_D])
taus = [S["h_decay"]["naive_dimensional"]["tau_s"], S["h_decay"]["loop_h_to_graviton_pair"]["band_tau_s"][0]] + \
       [r["tau_lower_s"] for r in S["h_decay"]["xi_trace_coupling"]["rows"]]
check("all in-action h lifetime estimates exceed 1e3 x age of Universe", min(taus) > 1e3 * T_UNIVERSE_S, min(taus))
tc = H["thermal_charge_capacity"]
check("thermal charge capacity n_Q,max/T^3 vs M/(3T) (rel < 1e-3)", abs(tc["nQ_max_over_T3"] / tc["leading_estimate_M_over_3T"] - 1) < 1e-3,
      tc["nQ_max_over_T3"] / tc["leading_estimate_M_over_3T"])
check("restoring power per quantum maximal at mu_N = 0 (F(z) monotonic)", H["restoration_number_bound"]["F_of_z"]["monotonic_increasing"], H["restoration_number_bound"]["F_of_z"]["F"])
check("kinetic-equilibrium window xi_kin < xi_max exists", B["restoration"]["window"]["kinetic_window_exists"], [B["restoration"]["window"]["xi_kinetic_min_phi"], xm["all_species"]])
check("number-changing equilibrium incompatible with relic bound (xi_chem,low > xi_max) [EXPECTED FALSE-for-compatibility]",
      not B["restoration"]["window"]["chemical_equilibrium_compatible"], B["restoration"]["window"])

rest = {(r["Tc_label"], r["xi_label"]): r["T_SM_min_GeV"] for r in B["restoration"]["rows"]}
headline = dict(
    module="relic_yield",
    conservation_verdict="Comoving net U(1)_X charge is exactly conserved by the unchanged action (Noether; no anomaly, "
                         "gauge-singlet scalar; checked to 4e-16 on an FRW lattice). Candidate 1 cannot generate Y_Q; it is an "
                         "outside-action initial datum supplied at a free time T_Q after the last ~14-20 inflationary e-folds.",
    N_max_efolds_r1=[r["N_max_efolds"] for r in C["inflation_dilution"]["rows"] if r["r"] == 1.0],
    h_lifetime_s=dict(naive_m5_over_MP4=S["h_decay"]["naive_dimensional"]["tau_s"],
                      loop_graviton_pair=S["h_decay"]["loop_h_to_graviton_pair"]["tau_s"],
                      xi_trace_xi_1em4=S["h_decay"]["xi_trace_coupling"]["rows"][1]["tau_s"]),
    Omega_sym_h2_vs_xi=[[r["xi"], r["Omega_h2"]] for r in B["scan_all_species"]],
    Omega_hot_relic_coef=B["hot_relic"]["coef_all_species_GeV"],
    xi_max=dict(all_species=xm["all_species"], no_D=xm["no_D_population"], hot_relic_band=[B["hot_relic"]["xi_max_all"], B["hot_relic"]["xi_max_no_D"]]),
    xi_32_equilibrium=B["xi_32_equilibrium"]["value"],
    Y_N_max=H["restoration_number_bound"]["Y_N_max"],
    T_SM_min_restoring_GeV=dict(leading_Tc=rest[("leading", "all_species")], resummed=[rest[("resummed_low", "all_species")], rest[("resummed_high", "all_species")]],
                                no_D=rest[("resummed_low", "no_D")], number_bound=[r["T_SM_min_thermal_number_GeV"] for r in H["restoration_number_bound"]["rows"]]),
    xi_BEC_formula=H["charge_in_thermal_plasma"]["xi_BEC_formula"],
    floor=dict(coef_GeV_per_omega_sym=H["floor_inequality"]["coefficient_GeV_per_unit_omega_sym"],
               r_max_primary_with_sym012=H["floor_inequality"]["r_max_primary_sym_012"],
               retrospective_DM_r_max=H["floor_inequality"]["retrospective_DM_bound_r_max"]),
    T_Q_max_cold_GeV={r["label"]: r["T_Q_max_nc_placeholder_GeV"] for r in H["charge_supply_cold_sector"]["rows"]},
    prior_art=[(p["ref"], p["Z2_handling"]) for p in H["prior_art"]],
    assertions_passed=sum(a["passed"] for a in asserts), assertions_total=len(asserts),
    timing_s=timing, total_runtime_s=time.time() - t0)
hashes = {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest() for f in sorted(os.listdir(HERE)) if f.endswith(".py")}
dump("headline.json", dict(headline=headline, assertions=asserts, code_sha256=hashes))
print(json.dumps(headline, indent=1, default=float))
for a in asserts:
    print("PASS" if a["passed"] else "FAIL", a["name"])
