#!/usr/bin/env python3
"""Collect the verifier's key numbers into vp_results.json (reads only this directory's JSON files)."""
import json, math, os
here = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(here, f)))
th, pr, sp, sr = L("vp_thermal.json"), L("vp_producer.json"), L("vp_spectator.json"), L("vp_sr_noise.json")
MP = 2.435323e18; sq = math.sqrt(math.pi**2 * 106.75 / 90)
rows = {r["z"]: r for r in th["Tc_of_z"]}
res = {
 "key": "v_production",
 "(1) Tc(z)": {
  "Pi_massless_common_mu": "Pi = (3 lphi + 2 gX^2 + hD^2) T_h^2 Li2(z)/(2 pi^2)  [derived]",
  "Tc1_leading_over_v": th["Tc1_leading_over_v"],
  "Tc1_hartree_over_v": rows[1.0]["Tc_hartree_over_v"],
  "hartree_ratio_to_z1": {z: rows[z]["ratio_hartree_to_hartree_z1"] for z in (0.9, 0.5, 0.1)},
  "leading_factor_sqrt_zeta2_over_Li2": {z: rows[z]["factor_leading"] for z in (0.9, 0.5, 0.1)},
  "hartree_over_leading_max": max(r["hartree_over_leading"] for r in th["Tc_of_z"]),
 },
 "(2) MB exit": {"K": th["MB"]["K"], "coef_GeV": th["MB"]["Tend_coefficient_GeV"],
                 "exact_over_MB_at_mu_-2": th["MB"]["exact_vs_MB_leading"][3]["Tend_exact_over_MB"],
                 "exact_over_MB_at_mu_0": th["MB"]["exact_vs_MB_leading"][0]["Tend_exact_over_MB"],
                 "zero_mode_response_Tend_max_GeV": math.sqrt(0.01 * 1e6 / (2 * math.sqrt(0.01) * 1000) * MP / sq)},
 "(3) T_SM,min": {"hartree_massless_n": th["mu0"]["TSMmin_hartree"], "resummed_range": th["mu0"]["TSMmin_resummed_range"],
                  "G1": th["mu0"]["G1_GeV"], "F1": th["mu0"]["F1"], "G_monotone": th["G_monotone_increasing_as_z_decreases"]},
 "(4) producer": {"C_He_large": pr["late_time"][4]["C"], "C_He_30Gamma": pr["late_time"][0]["C"],
                  "Y_coef_1p5C": pr["late_time"][4]["coef_1p5C"], "eps_check": [pr["late_time_B0.1"]["eps"], pr["late_time_B0.1"]["eps_expected"]],
                  "E_invE": pr["late_time"][4]["E_times_invE"], "A_RD_over_A_MD_He_1e12Gamma": pr["late_time"][4]["A_RD_over_A_MD"],
                  "Bmax_coef": pr["Bmax_coef"], "Bmin_coef_kappa1": pr["Bmin_coef_MD"],
                  "Bmin_cells_over_x2": [c["Bmin_over_x2"] for c in pr["Bmin_cells_Hubble_damped_kappa1"]],
                  "Bmin_kappa_9over4_over_x2": [c["Bmin_over_x2"] for c in pr["Bmin_cells_kappa_9over4"]],
                  "window_closure_x": pr["window_closes_at_mI_over_TRH_MD"], "c5_check": pr["c5_check"]["Bmin"],
                  "mu_over_T_lateRD_at_Bmax": {m["mI_over_TRH"]: m["mu_over_T_at_Bmax"] for m in pr["mu_at_Bmax"]}},
 "(5) spectator": {"moments": sp["equilibrium_moments_H_units"], "c_eps_no_velocity": sp["c_eps_minimal_no_velocity"],
                   "c_eps_slow_roll_exact": sp["c_eps_with_slow_roll_velocity"], "c_eps_slow_roll_MC2e6": sr["c_eps_SR_mean"],
                   "SEM_c_eps_1500_samples": sr["SEM_for_1500"], "frozen_over_potential": sp["frozen_over_potential_term"],
                   "ratio_to_critic": sp["ratio_to_critic"], "cosmic_time_check_ratios": [c["ratio"] for c in sp["cosmic_time_check"]],
                   "J": [sp["dyn"]["J_mean"], sp["dyn"]["J_sem"]], "shares": sp["dyn"]["shares"],
                   "mean_quantum_energy_over_H": sp["dyn"]["mean_quantum_energy_over_H"], "yields": sp["yields"]},
}
json.dump(res, open(os.path.join(here, "vp_results.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
