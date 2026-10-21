"""One-command reproduction of the thermal_potential module.
    cd modules/thermal_potential && OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 code/run_all.py
Runs every step (about 40 s on one core) and writes results/summary.json with the headline numbers and assertions."""
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
CODE = Path(__file__).resolve().parent
OUT = CODE.parent / "results"
STEPS = ["jb_checks.py", "transition.py", "xd_positivity.py", "kz_thermalization.py", "percolation.py"]


def main():
    t0 = time.time()
    timing = {}
    for s in STEPS:
        t = time.time()
        subprocess.run([sys.executable, str(CODE / s)], check=True, cwd=str(CODE))
        timing[s] = round(time.time() - t, 2)
    import numpy as np
    import scipy
    L = lambda n: json.loads((OUT / n).read_text())
    jb, tr, xd, kz, th, pc = (L("jb_checks.json"), L("transition.json"), L("xd_positivity.json"), L("kibble_zurek.json"),
                              L("thermalization.json"), L("percolation.json"))
    a = tr["adopted"]
    runs = tr["runs"]
    table = {}
    for r in runs:
        table.setdefault(r["rg"], {}).setdefault(r["scale_label"], {})[r["scheme"]] = {
            "Tc_over_v": r["Tc_over_v"], "order": r["order_perturbative"].split(" ")[0], "phic_over_Tc": r["phic_over_Tc"],
            "T0_over_v": r["T0_lower_spinodal_over_v"], "T1_over_v": r["T1_upper_spinodal_over_v"],
            "latent_heat_over_rho_hidden": r["latent_heat_over_rho_hidden"]}
    ol_orders = {r["order_perturbative"].split(" ")[0] for r in runs if r["scheme"] == "oneloop"}
    checks = {
        "JB_series_vs_quadrature_<1e-11": jb["checks_pass"]["series_matches_quadrature_1e-11"],
        "JB_bessel_vs_quadrature_<1e-11": jb["checks_pass"]["bessel_matches_quadrature_1e-11"],
        "JB_table_spline_vs_series_<1e-11": jb["checks_pass"]["spline_matches_series_1e-11"],
        "JB_cubic_coefficient_fit_-pi/6": jb["checks_pass"]["cubic_coefficient_within_1e-4"],
        "JB_log_coefficient_fit_-1/32": jb["checks_pass"]["log_coefficient_within_1e-3"],
        "Pi_phi/T^2_=_0.004175": abs(tr["leading_high_T"]["Pi_phi_over_T2"] - 0.004175) < 1e-15,
        "leading_T/v_=_1.547646465": abs(tr["leading_high_T"]["T_curvature_zero_over_v"] - 1.547646465) < 1e-9,
        "zero_T_phi_min_within_1e-6_of_archived_strict": abs(tr["zero_T_check"]["phi_min_over_v_numerical_V0+V1"] - 1.00012193740) < 1e-6,
        "Tc_independent_of_JB_evaluator_<1e-9": all(abs(v["diff"]) < 1e-9 for v in tr["evaluator_crosscheck_mu0"].values()),
        "X_and_D_curvatures_positive_along_phi_min_all_T": all(v["all_positive"] for v in xd["along_phi_minimum"].values()),
        "global_minimum_has_X=D=0_(|R|,|D|<1e-5)": all(g["max_|R|_or_|D|_at_refined_minima"] < 1e-5 for g in xd["global_search"]),
        "V_eff_even_in_phi_(4-field,<1e-15)": all(v["max_abs_V(phi)-V(-phi)_full_4field_random"] < 1e-15 for v in xd["parity_check"].values()),
        "DR_effective_quartic_positive_at_Tc": all(x > 0 for x in a["lambda_eff_over_lambda_phi_range"]["DR"]),
        "resummed_schemes_Tc_spread_<1%": (a["Tc_over_v_range_resummed_running"][1] / a["Tc_over_v_range_resummed_running"][0] - 1) < 0.01,
        "causal_bound_xi_KZ<=1/H_never_saturated": not kz["causal_bound_ever_active"],
        "both_vacua_percolate_p=1/2_(L=32,64,128)": all(v["all_span"] for v in pc["summary"].values()),
        "ROBUSTNESS_unresummed_order_scale_stable": len(ol_orders) == 1,
    }
    dec = {("%.0e" % r["xi"]): {"T_SM_GeV": r["T_SM_GeV"], "H_GeV": r["H_GeV"], "tau_Q_GeV^-1": r["tau_Q_GeV^-1"],
                                "log10_N_H_z1_matched": r["z1_ising_ginzburg_matched"]["log10_N_H"],
                                "log10_N_H_z1_naive": r["z1_ising_naive"]["log10_N_H"],
                                "log10_N_H_z2_matched": r["z2_ising_ginzburg_matched"]["log10_N_H"],
                                "log10_N_H_z2_naive": r["z2_ising_naive"]["log10_N_H"],
                                "log10_N_H_kibble_ginzburg": r["kibble_ginzburg"]["log10_N_H"],
                                "xi_KZ_z1_matched_GeV^-1": r["z1_ising_ginzburg_matched"]["xi_KZ_GeV^-1"],
                                "xi_KZ_z2_matched_GeV^-1": r["z2_ising_ginzburg_matched"]["xi_KZ_GeV^-1"],
                                "Hubble_length_GeV^-1": r["Hubble_length_GeV^-1"],
                                "eps_freeze_z1_matched": r["z1_ising_ginzburg_matched"]["eps_freeze"],
                                "phi_elastic_Gamma_over_H": [g for g in th["grid"] if abs(g["xi"] - r["xi"]) < 1e-12 * r["xi"]][0]["phi_elastic_Gamma_over_H"]}
           for r in kz["decades"]}
    summary = {
        "module": "thermal_potential",
        "status": "conditional calculation for the zero-charge-chemical-potential, internally thermal hidden bath; "
                  "no statement that Candidate 1 is excluded or validated",
        "Tc_table_over_v": table,
        "Tc_adopted_range_over_v": a["Tc_over_v_range_resummed_running"],
        "Tc_adopted_range_GeV": a["Tc_GeV_range_resummed_running"],
        "Tc_leading_highT_over_v": tr["leading_high_T"]["T_curvature_zero_over_v"],
        "shift_from_leading_percent": a["shift_from_leading_estimate_percent"],
        "two_loop_soft_Tc_shift_estimate_fraction": a["two_loop_soft_sunset_Tc_shift_estimate_fraction"],
        "order": {"resummed_AE_Parwani_DR": a["order_resummed_schemes"], "unresummed_oneloop": a["order_oneloop_unresummed"],
                  "nonperturbative": "continuous, 3D Ising universality (argument, not a lattice computation)"},
        "phic_over_Tc": {"resummed": 0.0, "oneloop_range": a["phic_over_Tc_oneloop_range"]},
        "latent_heat_over_rho_hidden_oneloop_range": a["latent_heat_over_rho_hidden_oneloop_range"],
        "expansion_parameters_at_Tc": {"eps_X": [r for r in runs if r["scheme"] == "DR"][0]["eps_X"],
                                       "eps_D": [r for r in runs if r["scheme"] == "DR"][0]["eps_D"],
                                       "eps_phi_soft_AE_Parwani": a["eps_phi_soft_at_Tc_range_AE_Parwani"],
                                       "lambda_eff_over_lambda_phi": a["lambda_eff_over_lambda_phi_range"]},
        "ginzburg_t_G": a["ginzburg"]["t_G"], "ginzburg_band": a["ginzburg"]["t_G_band_O1_criteria"],
        "kz_by_decade_of_xi": dec,
        "kz_log10_N_H_range_xi_1e-6_to_1": kz["log10_N_H_range_over_xi_grid"],
        "thermalization_xi_min": th["xi_min_for_Gamma>H"],
        "percolation_wall_area_fraction_between_spanning_clusters_L128": pc["summary"]["128"]["mean_wall_area_fraction_between_largest"],
        "checks": checks, "n_checks": len(checks), "n_pass": sum(bool(v) for v in checks.values()),
        "provenance": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                       "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(CODE.glob("*.py"))},
                       "step_seconds": timing, "total_seconds": round(time.time() - t0, 2)},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    cache = OUT / "jb_table.npz"          # binary quadrature cache used by transition.py; rebuilt on demand
    if cache.exists():
        cache.unlink()
    print("checks: %d/%d pass" % (summary["n_pass"], summary["n_checks"]))
    for k, v in checks.items():
        print("  %-55s %s" % (k, "PASS" if v else "FAIL"))
    print("total %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
