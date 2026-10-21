"""Run every calculation of module hidden_heating and write results/*.json (+ results/headline.json)."""
import os
import time
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from common import dump, tc, V_BENCH, M_P, c_H
import c1_couplings, c2_graviton_rates, c3_freezein_history, c4_inflation_sources, c5_branching_nonthermal

t0 = time.time()
r1 = c1_couplings.run()
r2 = c2_graviton_rates.run()
r3 = c3_freezein_history.run()
r4 = c4_inflation_sources.run()
r5 = c5_branching_nonthermal.run()
st = r4["stochastic"]
cond = r4["condensate"]["bounds"]
inst = r3["instantaneous"]
head = {
    "coupling_verdict": r1["coupling_verdict"],
    "beta_E_4dof_T9_over_MP4": r2["beta_E_4dof"],
    "beta_E_provenance": "computed here (explicit T_munu contraction, Bose/Fermi statistics, 4 Higgs + 45 Weyl + 12 vectors); "
                         "angular shapes agree with Clery et al. 2112.15214 eqs.(83)-(85); their eq.(86) beta_0 is "
                         "reproduced exactly only if the 45-Weyl channel is doubled (ratio ours/theirs = %.4f)" % r2["ratio_our_pair_rate_to_clery_beta0"],
    "pair_rate_per_real_dof_T8_over_MP4": r2["per_hidden_real_dof"]["pair_rate_coeff_T8_over_MP4"],
    "T_RH_max_instantaneous_GeV": inst["Tc_nominal"]["T_RH_max"],
    "T_RH_max_instantaneous_band_GeV": [inst["Tc_minus30"]["T_RH_max"], inst["Tc_plus30"]["T_RH_max"]],
    "T_RH_max_scaling": "∝ v^(4/7)",
    "T_max_max_noninstantaneous_GeV_by_H_e": {f"{row['H_e']:.0e}": row.get("T_max_at_boundary") for row in r3["noninstantaneous_bath_only"]},
    "T_RH_max_noninstantaneous_GeV_by_H_e": {f"{row['H_e']:.0e}": row.get("T_RH_max") for row in r3["noninstantaneous_bath_only"]},
    "selfthermalisation_Gamma22_over_H_at_T_RH_max": inst["thermalisation_at_T_RH_max"]["Gamma_22_over_H"],
    "nonthermal_restoration_T_RH_threshold_freezein_GeV": [inst["nonthermal_Pi_T_RH_threshold_noHubble"],
                                                            inst["nonthermal_Pi_T_RH_threshold_withHubble"]],
    "H_inf_dS_restoration_equilibrium_GeV": st["H_restoration_equilibrium"],
    "H_inf_origin_locally_stable_GeV": st["H_origin_locally_stable"],
    "central_charge_check": r2["central_charge_check"],
    "mode_function_check_ratio": r4["condensate"]["mode_function_check"]["ratio"],
    "condensate_numeric_vs_analytic_ratio": r3["condensate_numeric_check"]["ratio"],
    "T_max_max_noninst_v_slope": r3["noninstantaneous_v_scaling_H_e_1e13"]["slope_dlnTmax_dlnv"],
    "T_max_max_noninst_Tc_band": r3["noninstantaneous_Tc_band"],
    "H_inf_thresholds_bimodal_GeV": st["thresholds"],
    "H_inf_restoration_free_diffusion_over_v": st["H_restoration_over_v_after_N_efolds_free_diffusion"],
    "H_inf_formal_thermal_dS_energy_GeV": st["H_max_formal_thermal_from_dS_energy"],
    "H_e_bounds_condensate_scattering": cond,
    "B_h_thermal_bound_formula": r5["thermal_branching"]["formula"],
    "B_h_thermal_bound_rows": r5["thermal_branching"]["rows"],
    "B_h_nonthermal_rows": r5["nonthermal_branching"]["rows"],
    "runtime_s": time.time() - t0,
}
# H_e^(1/9) fit of the non-instantaneous T_max boundary
he = np.array([row["H_e"] for row in r3["noninstantaneous_bath_only"] if row.get("binding") is True])
tm = np.array([row["T_max_at_boundary"] for row in r3["noninstantaneous_bath_only"] if row.get("binding") is True])
sl, ic = np.polyfit(np.log10(he), np.log10(tm), 1)
head["T_max_max_fit"] = {"slope_dlog_dlogHe": sl, "prefactor_GeV": 10**ic,
                         "form": "T_max,max ≈ prefactor (H_e/GeV)^slope (v/TeV)^(4/9)"}
dump("headline.json", head)
print("done in %.1f s" % (time.time() - t0))
