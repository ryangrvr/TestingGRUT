"""v_relic: collect headline numbers from vr1-vr4 into results/vr_summary.json (no new physics)."""
import json
r1 = json.load(open("results/vr1_spectrum_amp.json"))
r2 = json.load(open("results/vr2_relic.json"))
r3 = json.load(open("results/vr3_bounds.json"))
r4 = json.load(open("results/vr4_threshold_scaling.json"))
b = r2["boltzmann_conv_only"]
S = {
 "item1_stability": {"mX": r1["masses"]["mX"], "mh": r1["masses"]["mh"], "mD": r1["masses"]["mD"],
   "h->XXbar_open": r1["closure"]["h->XXbar_open"], "h->DD_open": r1["closure"]["h->DD_open"],
   "tau_naive_s": r1["h_lifetime"]["tau_naive_s"], "tau_loop_hRiem2_s": r1["h_lifetime"]["tau_loop_s"],
   "tau_trace_xi1e-3_s": r1["h_lifetime"]["trace_xi"]["0.001"]["tau_s_x1"],
   "tau_trace_xi3.7e-5_s": r1["h_lifetime"]["trace_xi"]["3.7e-05"]["tau_s_x1"],
   "age_s": r1["h_lifetime"]["age_s"]},
 "item2_relic": {"Omega_coef_massless_frozen": r2["closed_form"]["Omega_coef_all_frozen_massless"],
   "Omega_frozen_conv_only_xi1e-3": b["conv_0.001"]["Omega"], "module_xi1e-3": 0.272,
   "Omega_frozen_conv_only_xi7.581e-4": b["conv_0.0007581"]["Omega"], "module_xi7.581e-4": 0.120,
   "Omega_frozen_composition_xi1e-3": r2["closed_form"]["Omega_frozen_ref_xi1e-3"],
   "Omega_full_h_to_X_xi1e-3": r2["closed_form"]["Omega_fullconv_ref_xi1e-3"],
   "xi_max_conv_only": b["xi_max_conv_only"], "module_xi_max": 7.581e-4,
   "xi_max_conv_only_noD": b["xi_max_conv_only_noD"], "module_xi_max_noD": 8.364e-4,
   "xi_max_bracket_frozen_vs_fullconv": [r2["closed_form"]["xi_max_frozen_ref"], r2["closed_form"]["xi_max_fullconv_ref"]],
   "Ytot_at_xi_max": b["at_xi_max"]["Ytot"], "mean_mass_at_xi_max": b["at_xi_max"]["mean_mass"],
   "collisionless_Y_drift": b["collisionless_1e-3"]["Y_rel_drift"]},
 "item3_threshold": {"M_thr_numeric": r1["threshold_amplitude"]["M_at_threshold_numeric_various_cos"][0],
   "diagrams": r1["threshold_amplitude"]["diagrams"], "scan": r1["threshold_amplitude"]["scan_a=g2/lam"],
   "sigv_vs_v_log_slope": r4["sigv_vs_vrel"][0]["log_slope_sigv_vs_v"],
   "sigv_hh->XX_GeV-2": {str(x["x=mX/T"]): x["sigv_GeV-2"] for x in r1["thermal_sigv"]["hh->XXbar"]}},
 "item4_TSMmin": {"routeA": r3["routeA_TSMmin=Tc/xi_max"], "routeB": r3["routeB_TSMmin_number"],
   "xi_c_over_xi_ref": r3["xi_c_over_xi_ref_number_matching"], "Pi_monotonic": r3["Pi_fixed_n_vs_z"]["monotonic"]},
 "item5_floor": r3["floor"],
}
json.dump(S, open("results/vr_summary.json", "w"), indent=1, default=float)
print(json.dumps(S, indent=1, default=float)[:3000])
