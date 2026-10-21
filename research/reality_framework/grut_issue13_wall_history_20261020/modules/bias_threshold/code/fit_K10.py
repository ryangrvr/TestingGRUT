"""fit_K10.py -- re-analysis of modules/network_sim biased-decay runs (read-only) and fits K10 = c m0^-alpha.

m0 = 2p - 1 is the initial magnetization of the sign field (|m0| = 2(1/2 - p)); network_sim fitted against
(1/2 - p).  K10 = xi_ref(tau_10)/xi_ref(tau_form): growth factor of the UNBIASED reference network's length
xi_ref = V/A_c = tau/(2 A_ref) between formation and the first time the biased network's area falls below 10% of
the reference.  network_sim took tau_10 as the first output below 10%; here tau_10 is also interpolated
(log A/A_ref linear in log tau between the bracketing outputs) and the output-grid bracket is kept as an error.

Outputs results/fit_K10.json.  [C] = computed here from network_sim's stored time series.
"""
from __future__ import annotations

import glob
import json
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NS = os.path.join(HERE, "..", "..", "network_sim", "results")
OUT = os.path.join(HERE, "..", "results")

GROUPS = {   # group prefix -> (d, tau_form as used by network_sim summary.json)
    "prs3d_bias_N128_L128": 3,
    "prs2d_bias_a2_N1024_L1024": 2,
    "prs2d_biasdamp_a2_N1024_L1024": 2,
}


def load_group(prefix):
    runs = {}
    for fn in glob.glob(os.path.join(NS, "runs", prefix + "_p*_s*.json")):
        m = re.search(r"_p([0-9.]+)_s(\d+)\.json$", fn)
        r = json.load(open(fn))
        t = np.array([e["tau"] for e in r["series"]])
        A = np.array([e["A"] for e in r["series"]])
        fp = np.array([e["f_plus"] for e in r["series"]])
        nl = np.array([e["nlinks"] for e in r["series"]])
        runs[(float(m.group(1)), int(m.group(2)))] = (t, A, fp, nl)
    return runs


def analyse_group(prefix, d, summary):
    runs = load_group(prefix)
    tau_form = summary["biased"][prefix]["tau_formation_ref"]
    ref = [v for (p, s), v in runs.items() if abs(p - 0.5) < 1e-9]
    tmin = max(v[0][0] for v in ref)
    tmax = min(v[0][-1] for v in ref)
    tref = np.geomspace(tmin, tmax, 256)
    Aref = np.mean([np.interp(np.log(tref), np.log(v[0]), v[1]) for v in ref], axis=0)

    def A_ref(t):
        return np.interp(np.log(np.clip(t, tref[0], tref[-1])), np.log(tref), Aref)

    def xi_ref(t):
        return t / (2 * A_ref(t))

    def K_of(t):
        return xi_ref(t) / xi_ref(tau_form)

    rows = []
    for (p, s), (t, A, fp, nl) in sorted(runs.items()):
        if abs(p - 0.5) < 1e-9:
            continue
        m0 = 2 * p - 1
        r = A / A_ref(t)
        ok = (t >= tau_form - 1e-9) & (t <= tref[-1])
        hit = np.nonzero(ok & (r < 0.1))[0]
        row = {"p": p, "seed": s, "m0": m0}
        if hit.size == 0:
            row.update({"K10": None, "censored": "not reached by box end"})
        else:
            i = int(hit[0])
            if t[i] <= tau_form + 1e-9 or i == 0 or not ok[i - 1]:
                row.update({"K10": 1.0, "K10_lo": 1.0, "K10_hi": 1.0, "censored": "below 10% at formation"})
            else:
                # interpolate in log r vs log t between outputs i-1 and i (r may be 0 at i: use first-output)
                r0, r1 = r[i - 1], r[i]
                if r1 > 0 and r0 > 0.1:
                    x = (math.log(0.1) - math.log(r0)) / (math.log(r1) - math.log(r0))
                    lt = math.log(t[i - 1]) + x * (math.log(t[i]) - math.log(t[i - 1]))
                else:
                    lt = math.log(t[i])
                t10 = math.exp(lt)
                row.update({"tau10_first_output": float(t[i]), "tau10_interp": t10,
                            "K10": float(K_of(t10)), "K10_first_output": float(K_of(t[i])),
                            "K10_lo": float(K_of(t[i - 1])), "K10_hi": float(K_of(t[i])), "censored": None})
        # magnetization at formation (interpolated in log tau): the population bias of the formed domains
        row["abs_m_form"] = float(abs(np.interp(math.log(tau_form), np.log(t), 2 * fp - 1)))
        # magnetization trajectory for the scaling-collapse test
        sel = (t >= tau_form) & (t <= tref[-1])
        row["traj"] = {"K": [float(k) for k in K_of(t[sel])], "abs_m": [float(abs(2 * f - 1)) for f in fp[sel]]}
        rows.append(row)
    return {"tau_form": tau_form, "d": d, "rows": rows,
            "xi_ref_form": float(xi_ref(tau_form)), "K_box_end": float(K_of(tref[-1]))}


def wls_powerlaw(m0, K, sig_lnK):
    """ln K = ln c - alpha ln m0, weighted least squares; returns c, alpha and their 1-sigma errors and cov."""
    x = np.log(m0)
    y = np.log(K)
    w = 1.0 / np.asarray(sig_lnK) ** 2
    X = np.vstack([np.ones_like(x), -x]).T
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    beta = cov @ (X.T @ (w * y))
    res = y - X @ beta
    chi2 = float(np.sum(w * res ** 2))
    dof = len(x) - 2
    scale = max(1.0, chi2 / dof) if dof > 0 else 1.0
    cov_s = cov * scale
    return {"c": float(math.exp(beta[0])), "alpha": float(beta[1]),
            "sig_lnc": float(math.sqrt(cov_s[0, 0])), "sig_alpha": float(math.sqrt(cov_s[1, 1])),
            "cov_lnc_alpha": float(cov_s[0, 1]), "chi2": chi2, "dof": dof}


def fit_points(g, m0_max=None):
    """Seed-mean K10 per m0 with an error combining the output-grid bracket and the seed scatter."""
    byp = {}
    for r in g["rows"]:
        if r["K10"] is None or r["censored"]:
            continue
        byp.setdefault(round(abs(r["m0"]), 6), []).append(r)
    pts = []
    for m0, rs in sorted(byp.items()):
        if m0_max is not None and m0 > m0_max + 1e-9:
            continue
        lnK = np.array([math.log(r["K10"]) for r in rs])
        half = np.array([0.5 * (math.log(r["K10_hi"]) - math.log(r["K10_lo"])) for r in rs])
        sig = math.sqrt(np.mean(half ** 2) / 3.0 + (np.var(lnK, ddof=1) / len(rs) if len(rs) > 1 else 0.0))
        sig = max(sig, 0.03)
        pts.append({"m0": m0, "K10_mean": float(math.exp(lnK.mean())), "sig_lnK": sig, "n": len(rs),
                    "K10_first_output_mean": float(np.mean([r["K10_first_output"] for r in rs]))})
    return pts


def formation_fit(g):
    """Fit ln K10 = (ln x10 - ln y_f)/gamma, y_f = erfinv(|m(tau_form)|) (erf-linearized formation bias).
    Scaling form m(K) = erf(y_f K^gamma) (Bray §5.2 eq. (175) with t^(d/4) -> L^(d - lambda_bar)); decay (A < 10%
    of reference) when y_f K^gamma = x10.  Uses every uncensored run (all seeds) with K10 > 1."""
    from scipy.special import erfinv
    pts = [(float(erfinv(min(r["abs_m_form"], 0.999999))), r["K10"], r["K10_lo"], r["K10_hi"], r["m0"], r["seed"])
           for r in g["rows"] if r["K10"] is not None and not r["censored"] and r["K10"] > 1.0]
    if len(pts) < 3:
        return None
    y = np.array([q[0] for q in pts])
    K = np.array([q[1] for q in pts])
    sig = np.array([max(0.5 * math.log(q[3] / q[2]) / math.sqrt(3), 0.03) for q in pts])
    X = np.vstack([np.ones_like(y), np.log(y)]).T
    w = 1 / sig ** 2
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    b = cov @ (X.T @ (w * np.log(K)))
    res = np.log(K) - X @ b
    chi2 = float(np.sum(w * res ** 2))
    dof = len(y) - 2
    cov = cov * max(1.0, chi2 / dof)
    gamma = -1.0 / b[1]
    sig_gamma = math.sqrt(cov[1, 1]) / b[1] ** 2
    lnx10 = b[0] * gamma
    # x10 uncertainty by propagation (lnx10 = -b0/b1)
    J = np.array([-1 / b[1], b[0] / b[1] ** 2])
    sig_lnx10 = float(math.sqrt(J @ cov @ J))
    # x10 at fixed gamma (median over runs)
    x10_fixed = {f"{gm:.3f}": float(np.exp(np.median(np.log(y) + gm * np.log(K)))) for gm in (1.25, 1.327, 1.39, 1.4, 1.5, 0.75, 0.85, 1.0)}
    return {"gamma": float(gamma), "sig_gamma": float(sig_gamma), "x10": float(math.exp(lnx10)),
            "sig_lnx10": sig_lnx10, "cov_lnx10_gamma_note": "propagated from the (b0, b1) covariance",
            "chi2": chi2, "dof": dof, "n_runs": len(pts),
            "points": [{"m0": q[4], "seed": q[5], "y_f": q[0], "K10": q[1]} for q in pts],
            "x10_at_fixed_gamma": x10_fixed}


def collapse_gamma(g, gammas):
    """Data collapse of |m|(K) curves of different m0 onto one function of x = |m0| K^gamma (nonlinear regime)."""
    curves = []
    for r in g["rows"]:
        K = np.array(r["traj"]["K"])
        am = np.array(r["traj"]["abs_m"])
        if K.size < 3:
            continue
        curves.append((abs(r["m0"]), K, am))
    best = None
    scores = []
    for gm in gammas:
        xs = [np.log(m0 * K ** gm) for m0, K, _ in curves]
        lo = max(x.min() for x in xs)
        hi = min(x.max() for x in xs)
        # use the overlap window and |m| < 0.9 (before the network has gone)
        grid = np.linspace(lo, hi, 25) if hi > lo else None
        if grid is None:
            scores.append(None)
            continue
        Y = []
        for (m0, K, am), x in zip(curves, xs):
            order = np.argsort(x)
            Y.append(np.interp(grid, x[order], am[order]))
        Y = np.array(Y)
        keep = Y.mean(axis=0) < 0.9
        sc = float(np.mean(np.var(Y[:, keep], axis=0))) if keep.any() else None
        scores.append(sc)
        if sc is not None and (best is None or sc < best[1]):
            best = (gm, sc)
    return {"gammas": list(gammas), "spread": scores, "best_gamma": best[0] if best else None}


def main():
    summary = json.load(open(os.path.join(NS, "summary.json")))
    out = {"note": "m0 = 2p - 1; K10 = xi_ref(tau_10)/xi_ref(tau_form); data from modules/network_sim results/runs "
                   "(read only). Fits are over the tested range |m0| = 0.1-0.6 only."}
    theory = {"model_A_3D_lambda_bar_1.60": 1 / (3 - 1.60), "model_A_3D_lambda_bar_1.673_LiuMazenko": 1 / (3 - 1.673),
              "gaussian_OJK_Hindmarsh_3D_lambda_bar_1.5": 1 / 1.5,
              "model_A_2D_lambda_bar_1.25": 1 / (2 - 1.25), "model_A_2D_LiuMazenko_1.289": 1 / (2 - 1.289),
              "gaussian_2D_lambda_bar_1": 1.0}
    out["alpha_theory = 1/(d - lambda_bar)"] = theory
    for prefix, d in GROUPS.items():
        g = analyse_group(prefix, d, summary)
        pts = fit_points(g)
        res = {"tau_form": g["tau_form"], "xi_ref_form": g["xi_ref_form"], "K_box_end": g["K_box_end"],
               "points": pts,
               "per_run": [{k: v for k, v in r.items() if k != "traj"} for r in g["rows"]]}
        if len(pts) >= 2:
            m0 = np.array([q["m0"] for q in pts])
            K = np.array([q["K10_mean"] for q in pts])
            s = np.array([q["sig_lnK"] for q in pts])
            res["fit_all"] = wls_powerlaw(m0, K, s)
            # first-output K10 (network_sim's definition) for comparison
            Kf = np.array([q["K10_first_output_mean"] for q in pts])
            res["fit_all_first_output_def"] = wls_powerlaw(m0, Kf, s)
            # jackknife over m0 points
            jk = []
            for i in range(len(pts)):
                sel = np.arange(len(pts)) != i
                if sel.sum() >= 2:
                    f = wls_powerlaw(m0[sel], K[sel], s[sel])
                    jk.append({"dropped_m0": float(m0[i]), "alpha": f["alpha"], "c": f["c"]})
            res["jackknife"] = jk
            # smallest-bias pair: local slope (least affected by the K10 >= 1 floor)
            res["local_slope_smallest_two"] = float(math.log(K[0] / K[1]) / math.log(m0[1] / m0[0]))
            res["sig_local_slope_smallest_two"] = float(math.hypot(s[0], s[1]) / math.log(m0[1] / m0[0]))
            # alpha fixed at theory values: fit c only, using the smallest-bias point(s) (m0 <= 0.2)
            fixed = {}
            for lab, al in theory.items():
                if (d == 3) != ("3D" in lab):
                    continue
                sel = m0 <= 0.2 + 1e-9
                lnc = np.sum((np.log(K[sel]) + al * np.log(m0[sel])) / s[sel] ** 2) / np.sum(1 / s[sel] ** 2)
                fixed[lab] = {"alpha": al, "c": float(math.exp(lnc)),
                              "max_abs_lnK_residual_all_points": float(np.max(np.abs(np.log(K) - (lnc - al * np.log(m0)))))}
            res["fit_c_at_theory_alpha"] = fixed
        res["collapse"] = collapse_gamma(g, np.round(np.arange(0.5, 2.01, 0.05), 3))
        res["formation_fit"] = formation_fit(g)
        out[prefix] = res
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "fit_K10.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    o = main()
    for prefix in GROUPS:
        r = o[prefix]
        print(prefix, "tau_form", r["tau_form"], "K_box_end %.2f" % r["K_box_end"])
        for q in r["points"]:
            print("   m0 %.2f K10 %.3f (first-output %.3f) sig_lnK %.3f n %d" % (q["m0"], q["K10_mean"],
                  q["K10_first_output_mean"], q["sig_lnK"], q["n"]))
        if "fit_all" in r:
            f = r["fit_all"]
            print("   fit: c = %.3f (sig_lnc %.3f), alpha = %.3f +- %.3f, chi2/dof %.2f/%d" % (
                f["c"], f["sig_lnc"], f["alpha"], f["sig_alpha"], f["chi2"], f["dof"]))
            f = r["fit_all_first_output_def"]
            print("   first-output def: c = %.3f alpha = %.3f" % (f["c"], f["alpha"]))
            print("   jackknife:", [(j["dropped_m0"], round(j["alpha"], 3), round(j["c"], 3)) for j in r["jackknife"]])
            print("   local slope smallest two: %.3f +- %.3f" % (r["local_slope_smallest_two"], r["sig_local_slope_smallest_two"]))
            print("   c at theory alpha:", {k: (round(v["alpha"], 3), round(v["c"], 3), round(v["max_abs_lnK_residual_all_points"], 3)) for k, v in r["fit_c_at_theory_alpha"].items()})
        print("   collapse best gamma:", r["collapse"]["best_gamma"])
        ff = r["formation_fit"]
        if ff:
            print("   formation fit: gamma = %.3f +- %.3f, x10 = %.3f (sig_ln %.3f), chi2/dof %.2f/%d, n %d" % (
                ff["gamma"], ff["sig_gamma"], ff["x10"], ff["sig_lnx10"], ff["chi2"], ff["dof"], ff["n_runs"]))
            print("   y_f, K10:", [(round(q["m0"], 2), round(q["y_f"], 3), round(q["K10"], 3)) for q in ff["points"]])
            print("   x10 at fixed gamma:", ff["x10_at_fixed_gamma"])
