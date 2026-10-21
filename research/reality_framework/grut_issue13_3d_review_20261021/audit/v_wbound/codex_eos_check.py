#!/usr/bin/env python3
"""Verifier v_wbound: re-read Codex a69bd86 raw diagnostics (read-only extract passed as argv[1]) and recompute the
equation of state independently of wbound/codex_eos.py: w is rebuilt from the K_phys, G_phys, U_phys columns
(rho = K + G + U, p = K - G/3 - U) and compared with Codex's own rho, p columns; rho a^4 ratios, monotonicity and the
integrated continuity check ln(R_end/R_0) vs int (1 - 3w) dln a are reported.  Writes results/codex_eos_check.json.
"""
import json
import os
import sys

import numpy as np

RUNS = ["main_s101", "main_s202", "main_s303", "cold_s101", "kappa4_s101", "dt_half_s101", "grid256_s101",
        "grid384_s101", "large_s101", "gaussian_s101", "gaussian_s202", "gaussian_large_s101"]


def main(raw):
    out = {}
    for r in RUNS:
        arr = np.genfromtxt(os.path.join(raw, r + ".csv"), delimiter=",", names=True)
        a, K, G, U = arr["a"], arr["K_phys"], arr["G_phys"], arr["U_phys"]
        rho = K + G + U
        p = K - G / 3 - U
        w = p / rho
        R = rho * a**4
        lna = np.log(a)
        integ = np.trapezoid(1 - 3 * w, lna)
        out[r] = dict(n=int(len(a)), a_end=float(a[-1]),
                      max_rel_rho_vs_codex=float(np.max(np.abs(rho / arr["rho"] - 1))),
                      max_rel_p_vs_codex=float(np.max(np.abs(p - arr["p"]) / arr["rho"])),
                      w_min=float(w.min()), w_max=float(w.max()), n_w_gt_third=int(np.sum(w > 1 / 3)),
                      R_ratio=float(R[-1] / R[0]), R_monotone=bool(np.all(np.diff(R) > 0)),
                      ln_R_ratio=float(np.log(R[-1] / R[0])), int_1m3w_dlna=float(integ),
                      wbar_rho_dlna=float(np.trapezoid(p, lna) / np.trapezoid(rho, lna)))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "results", "codex_eos_check.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for k, v in out.items():
        print(f"{k:20s} n={v['n']} w[{v['w_min']:+.3f},{v['w_max']:+.3f}] R x{v['R_ratio']:.3f} mono={v['R_monotone']} "
              f"lnR={v['ln_R_ratio']:.3f} int(1-3w)={v['int_1m3w_dlna']:.3f} drho={v['max_rel_rho_vs_codex']:.1e} dp={v['max_rel_p_vs_codex']:.1e}")


if __name__ == "__main__":
    main(sys.argv[1])
