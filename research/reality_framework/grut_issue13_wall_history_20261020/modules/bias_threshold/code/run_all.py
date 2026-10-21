"""run_all.py -- module bias_threshold: lattice response runs (2 single-threaded workers, ~4 min wall, ~8 CPU-min),
then analysis, K10 fits, thresholds and the headline.  Existing results/lattice_runs/*.json are reused.

    OMP_NUM_THREADS=1 python3 code/run_all.py
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")


def main():
    t0 = time.time()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    procs = [subprocess.Popen([sys.executable, os.path.join(HERE, "lattice.py"), str(w), "2"], env=env) for w in (0, 1)]
    for p in procs:
        p.wait()
    t_lat = time.time() - t0
    for script in ("analyze_lattice.py", "fit_K10.py", "thresholds.py"):
        subprocess.run([sys.executable, os.path.join(HERE, script)], env=env, check=True, stdout=subprocess.DEVNULL)
    lat = json.load(open(os.path.join(RES, "lattice_response.json")))
    fk = json.load(open(os.path.join(RES, "fit_K10.json")))
    th = json.load(open(os.path.join(RES, "thresholds.json")))
    f3 = fk["prs3d_bias_N128_L128"]
    hs = th["histories"]["H_star"]
    hq = th["histories"]["H_star_Q"]
    head = {
        "labels": "[C] computed here; [D] derived; [Q] quoted (source); [P] preparation input outside the action; "
                  "[E] estimate.  m = 2p - 1 (systematic population bias); thresholds are CONDITIONAL on the scaling "
                  "theory holding over 11-17 decades in L/xi0 beyond the runs.",
        "K10 = c m0^-alpha, 3D PRS 128^3, m0 = initial (quench-start) bias [C, network_sim data]": {
            "alpha": f3["fit_all"]["alpha"], "sig_alpha": f3["fit_all"]["sig_alpha"], "c": f3["fit_all"]["c"],
            "sig_lnc": f3["fit_all"]["sig_lnc"], "first_output_definition": {"alpha": f3["fit_all_first_output_def"]["alpha"],
                                                                            "c": f3["fit_all_first_output_def"]["c"]},
            "jackknife_alpha": [j["alpha"] for j in f3["jackknife"]],
            "implied_gamma = 1/alpha": 1 / f3["fit_all"]["alpha"]},
        "K10 = c m0^-alpha, 2D PRS (m0 = 2(1/2 - p)) [C, network_sim data]": {
            "alpha": fk["prs2d_bias_a2_N1024_L1024"]["fit_all"]["alpha"], "c": fk["prs2d_bias_a2_N1024_L1024"]["fit_all"]["c"],
            "damped_late_start_alpha": fk["prs2d_biasdamp_a2_N1024_L1024"]["fit_all"]["alpha"]},
        "gamma = d - lambda_bar, 3D relativistic (PRS) [C]": {
            "collapse_of_biased_runs (formation bias, erf-linearized)": [f3["formation_fit"]["gamma"], f3["formation_fit"]["sig_gamma"]],
            "x10": f3["formation_fit"]["x10"],
            "unbiased_autocorrelation_lambda_bar_safe/ext": [lat["prs3d_N160"]["lambda_bar_safe"]["mean"],
                                                              lat["prs3d_N160"]["lambda_bar_ext"]["mean"]],
            "direct_small_bias_response_gamma_safe (K <= 2.6)": [lat["prs3d_N160"]["gamma_Ry_safe_m0=0.002"]["mean"],
                                                                 lat["prs3d_N160"]["gamma_Ry_safe_m0=0.002"]["std"]],
            "pre_formation_linear_stage_gamma (theory d/2 = 1.5)": lat["prs3d_N160"]["gamma_pre_formation = ln R_y(form)/ln(xi_form/xi_i)_m0=0.002"]["mean"],
            "model_A_literature_gamma": "1.33 (Liu-Mazenko lambda_bar 1.673, Bray Table 1) - 1.40 (lambda_bar 1.6, brief)",
            "band_used": th["calibration"]["gamma_band"]},
        "gamma 2D [C]": {"relativistic collapse (early/damped)": [fk["prs2d_bias_a2_N1024_L1024"]["formation_fit"]["gamma"],
                                                                 fk["prs2d_biasdamp_a2_N1024_L1024"]["formation_fit"]["gamma"]],
                         "model A, this code, K <= 7: 2 - lambda_bar": lat["modelA2d_N512"]["gamma_from_lambda_safe"]["mean"],
                         "model A literature asymptotic": 0.75},
        "H_star [P: xi = 5e-4, T_f = 3.1e6 GeV]": {
            "log10_K_par": [hs["central"]["log10_K_par"], hs["band"]["log10_K_par"]],
            "log10_K_safe": [hs["central"]["log10_K_safe"], hs["band"]["log10_K_safe"]],
            "T_safe_MeV(f=1..0.01)": hs["band"]["T_safe_MeV"],
            "eps_persist": [hs["central"]["eps_persist"], hs["band"]["log10_eps_persist"]],
            "eps_safe": [hs["central"]["eps_safe"], hs["band"]["log10_eps_safe"]]},
        "H_star_Q": {"eps_persist": [hq["central"]["eps_persist"], hq["band"]["log10_eps_persist"]],
                     "eps_safe": [hq["central"]["eps_safe"], hq["band"]["log10_eps_safe"]],
                     "log10_K_par": hq["band"]["log10_K_par"]},
        "inflationary randomization (H_inf = 2.4e4 - 1e13 GeV [P])": {
            k: {"T_f_GeV": v["central"]["T_f_GeV"], "eps_persist_log10_band": v["band"]["log10_eps_persist"],
                "eps_safe_log10_band": v["band"]["log10_eps_safe"]}
            for k, v in th["histories"].items() if k.startswith("infl")},
        "extrapolation_decades_in_L/xi0": th["extrapolation_decades"],
        "statistical_vs_systematic": th["statistical_vs_systematic"],
        "candidate_residual_biases_H_star": th["candidates_H_star"],
        "damping_e-folds_needed_for_1e-3_remnant": th["damping_required_H_star"],
        "comparison": th["comparison"],
        "runtime_s": {"lattice_wall (2 workers; 0 if reused)": t_lat, "total": time.time() - t0},
    }
    with open(os.path.join(RES, "headline.json"), "w") as f:
        json.dump(head, f, indent=1)
    print(json.dumps(head["runtime_s"]))


if __name__ == "__main__":
    main()
