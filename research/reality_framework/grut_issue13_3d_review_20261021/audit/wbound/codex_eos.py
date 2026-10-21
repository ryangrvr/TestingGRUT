#!/usr/bin/env python3
"""Track wbound: equation of state in Codex's own 3D charged-exit diagnostics (a69bd86, raw_diagnostics.tar.gz).
Read-only; pass the extracted raw directory as argv[1].  Columns used: t (conformal), a, rho, p (all fields,
physical, units lambda_phi v^4), K_phys, G_phys, U_phys, U_phi.
For each production run: instantaneous w range, fraction of samples with w > 1/3, ln a-weighted mean w, and
rho a^4 (end/start).  d ln(rho a^4)/d ln a = 1 - 3w for a closed sector, so rho a^4 non-decreasing <=> wbar <= 1/3.
Samples are 1.2 conformal-time units apart (41 per run), so averages are coarse.  Writes results/codex_eos.json.
"""
import csv
import json
import os
import sys

import numpy as np

RUNS = ["main_s101", "main_s202", "main_s303", "cold_s101", "kappa4_s101", "dt_half_s101", "grid256_s101",
        "grid384_s101", "large_s101", "gaussian_s101", "gaussian_s202", "gaussian_large_s101"]


def load(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def main(rawdir):
    out = {}
    for name in RUNS:
        d = load(os.path.join(rawdir, name + ".csv"))
        a, rho, p = d["a"], d["rho"], d["p"]
        w = p / rho
        lna = np.log(a)
        wbar_lna = float(np.trapezoid(p, lna) / np.trapezoid(rho, lna))
        R = rho * a**4
        # continuity check: d ln(rho a^4)/d ln a vs 1 - 3w (finite differences, coarse)
        dlnR = np.diff(np.log(R)) / np.diff(lna)
        wmid = 0.5 * (w[1:] + w[:-1])
        out[name] = dict(w_first=float(w[0]), w_last=float(w[-1]), w_min=float(w.min()), w_max=float(w.max()),
                         n_samples=int(len(w)), frac_w_gt_third=float(np.mean(w > 1 / 3)),
                         wbar_lna_weighted=wbar_lna, rho_a4_end_over_start=float(R[-1] / R[0]),
                         rho_a4_nondecreasing_all_samples=bool(np.all(np.diff(R) >= 0)),
                         max_abs_continuity_mismatch=float(np.max(np.abs(dlnR - (1 - 3 * wmid)))),
                         U_phi_over_rho_last=float(d["U_phi"][-1] / rho[-1]),
                         K_over_rho_last=float(d["K_phys"][-1] / rho[-1]),
                         G_over_rho_last=float(d["G_phys"][-1] / rho[-1]),
                         a_end=float(a[-1]))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "results", "codex_eos.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for k, r in out.items():
        print(f"{k:22s} w[{r['w_min']:+.3f},{r['w_max']:+.3f}] last {r['w_last']:+.3f} wbar {r['wbar_lna_weighted']:+.3f} "
              f"rho_a4 x{r['rho_a4_end_over_start']:.3f} mono {r['rho_a4_nondecreasing_all_samples']} "
              f"cont {r['max_abs_continuity_mismatch']:.3f} >1/3 {r['frac_w_gt_third']:.2f}")


if __name__ == "__main__":
    main(sys.argv[1])
