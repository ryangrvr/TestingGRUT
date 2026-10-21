"""Item 2 and the parity check of item 3.
(a) Along the phi minimum phi_min(T) of each scheme, compute the full one-loop curvatures in the X (q1) and D
    directions by finite differences of the four-field effective potential (4x4 Hessian eigenvalues, exact J_B).
(b) Global grid search over (phi, |X|, D) at several T (AE and Parwani), with local refinement from off-axis starts.
(c) Parity: V_full(phi, q1, D) - V_full(-phi, q1, D) at random points (the off-diagonal 2 g_X^2 phi q1 flips sign).
Writes results/xd_positivity.json."""
import json
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from common import BENCH, MU0, Veff, JBSeries, V_full, debye_coeffs, to_builtin

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent.parent / "results"
p = BENCH
jb = JBSeries()


def phi_min_1d(V, T):
    phi = np.concatenate([[0.0], np.geomspace(1e-6, 1.6, 4000)])
    w = V.W(phi**2, T)
    i = int(np.argmin(w))
    if i == 0:
        return 0.0
    lo, hi = phi[max(i - 1, 0)], phi[min(i + 1, len(phi) - 1)]
    r = minimize(lambda x: V.W(x[0] ** 2, T)[0], [phi[i]], bounds=[(lo, hi)], method="L-BFGS-B",
                 options={"ftol": 1e-16, "gtol": 1e-14})
    return float(r.x[0])


def curv(scheme, T, phi0, which, h=1e-3):
    q1 = np.array([0.0, h]) if which == "X" else np.array([0.0, 0.0])
    D = np.array([0.0, h]) if which == "D" else np.array([0.0, 0.0])
    v = V_full(p, MU0, scheme, jb, np.array([phi0, phi0]), q1, D, T)
    return 2 * (v[1] - v[0]) / h**2


def main():
    t0 = time.time()
    res = {"scale": "mu = mu0, couplings at mu0"}
    Tc_tr = json.loads((OUT / "transition.json").read_text())
    tc = {r["scheme"]: r["Tc_over_v"] for r in Tc_tr["runs"] if r["rg"] == "running" and r["scale_label"] == "mu0"}
    pic = debye_coeffs(p)
    along = {}
    for sch in ("oneloop", "AE", "Parwani"):
        V = Veff(p, MU0, sch, jb)
        Ts = sorted(set(np.round(np.linspace(0.0, 3.0, 31), 6).tolist() + [tc[sch] * (1 - 1e-3), tc[sch], tc[sch] * (1 + 1e-3)]))
        rows = []
        for T in Ts:
            ph = phi_min_1d(V, T)
            mX = curv(sch, T, ph, "X")
            mD = curv(sch, T, ph, "D")
            lead_X = p["c"] * ph**2 + pic[1] * T**2
            lead_D = p["md"] + p["d"] * ph**2 + pic[2] * T**2
            rows.append({"T_over_v": T, "phi_min_over_v": ph, "m2_X_eff_over_v2": mX, "m2_D_eff_over_v2": mD,
                         "leading_m2_X": lead_X, "leading_m2_D": lead_D})
        along[sch] = {"rows": rows, "min_m2_X": min(r["m2_X_eff_over_v2"] for r in rows),
                      "min_m2_D": min(r["m2_D_eff_over_v2"] for r in rows),
                      "all_positive": all(r["m2_X_eff_over_v2"] > 0 and r["m2_D_eff_over_v2"] > 0 for r in rows)}
        print(sch, "min m2_X %.4e  min m2_D %.4e  all positive %s" % (along[sch]["min_m2_X"], along[sch]["min_m2_D"],
                                                                       along[sch]["all_positive"]))
    res["along_phi_minimum"] = along

    # (b) global grid + refinement
    g_phi = np.linspace(0, 1.6, 49)
    g_R = np.linspace(0, 1.6, 49)
    g_D = np.linspace(0, 1.6, 33)
    P, R, Dd = np.meshgrid(g_phi, g_R, g_D, indexing="ij")
    glob = []
    for sch in ("AE", "Parwani"):
        for T in [0.0, 0.5, 1.0, 1.4, tc[sch], 1.7, 2.5]:
            Vg = V_full(p, MU0, sch, jb, P, R, Dd, T)
            i = np.unravel_index(np.argmin(Vg), Vg.shape)
            off = (R > 0) | (Dd > 0)
            gap = Vg[off].min() - Vg.min()
            f = lambda x: V_full(p, MU0, sch, jb, np.array([x[0]]), np.array([x[1]]), np.array([x[2]]), T)[0]
            starts = [[g_phi[i[0]], 0.3, 0.3], [0.2, 0.8, 0.5], [0.05, 0.2, 0.0], [0.05, 0.0, 0.2], [0.9, 0.5, 0.0]]
            ends = []
            for st in starts:
                r = minimize(f, st, method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-16, "maxiter": 20000, "maxfev": 40000})
                ends.append([abs(r.x[0]), abs(r.x[1]), abs(r.x[2]), r.fun])
            glob.append({"scheme": sch, "T_over_v": T, "grid_argmin_(phi,R,D)": [g_phi[i[0]], g_R[i[1]], g_D[i[2]]],
                         "grid_min_off_axis_minus_global": gap, "refined_minima_(|phi|,|R|,|D|,V)": ends,
                         "max_|R|_or_|D|_at_refined_minima": max(max(e[1], e[2]) for e in ends)})
            print(sch, "T=%.4f grid argmin" % T, [g_phi[i[0]], g_R[i[1]], g_D[i[2]]], "max |R|,|D| refined %.2e" % glob[-1]["max_|R|_or_|D|_at_refined_minima"])
    res["global_search"] = glob

    # (c) parity check on the full four-field potential
    rng = np.random.default_rng(13)
    n = 2000
    ph = rng.uniform(-1.5, 1.5, n)
    q = rng.uniform(-1.5, 1.5, n)
    D = rng.uniform(-1.5, 1.5, n)
    par = {}
    for sch in ("oneloop", "AE", "Parwani"):
        worst = 0.0
        scale = 0.0
        for T in [0.0, 0.7, 1.4, 1.56, 2.2]:
            a = V_full(p, MU0, sch, jb, ph, q, D, T)
            b = V_full(p, MU0, sch, jb, -ph, q, D, T)
            worst = max(worst, float(np.max(np.abs(a - b))))
            scale = max(scale, float(np.max(np.abs(a))))
        V = Veff(p, MU0, sch, jb)
        x = np.geomspace(1e-6, 2, 500)
        par[sch] = {"max_abs_V(phi)-V(-phi)_full_4field_random": worst, "max_abs_V_on_sample": scale,
                    "max_abs_1d_background": float(np.max(np.abs(V.V_phi(x, 1.5) - V.V_phi(-x, 1.5))))}
    res["parity_check"] = par
    res["symmetry_statement"] = (
        "U and every scheme's V_eff depend on phi only through phi^2 (and through phi in the 4x4 Hessian only via the "
        "off-diagonal entries 2 g_X^2 phi q1 and 2 h_D^2 phi D, whose sign flip is an orthogonal similarity that leaves the "
        "eigenvalues unchanged).  For a Z2-invariant initial density matrix (e.g. any thermal state above Tc, rho ~ exp(-H/T)) "
        "every Z2-odd expectation vanishes, so the probability that a correlation volume selects +v equals that of -v: "
        "p(+) = p(-) = 1/2 exactly, independent of the order of the transition.  This is a symmetry statement about the "
        "ensemble, not a dynamical result; it fails only if the state entering the transition is not Z2 invariant "
        "(preparation data outside the action).")
    res["seconds"] = time.time() - t0
    (OUT / "xd_positivity.json").write_text(json.dumps(to_builtin(res), indent=1))
    print(json.dumps(to_builtin(par), indent=1), "time %.1fs" % res["seconds"])


if __name__ == "__main__":
    main()
