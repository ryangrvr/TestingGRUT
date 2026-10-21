"""analyze_lattice.py -- growth exponent gamma = d - lambda_bar of a small systematic bias, from lattice.py runs.

For each (group, seed): reference run m0 = 0 gives xi_ref(t) = V/A_c and the formation time; the paired runs +-m0
give the response.  Two response estimators:
  R_m(t) = (m_+ - m_-)/(2 m0)                         (raw; saturates once |m| is not small)
  R_y(t) = (erfinv m_+ - erfinv m_-)/(2 erfinv m0)    (erf-linearized; Bray §5.2 eq. (175) form <s> = erf(y))
and the local autocorrelation C(t) = <sign phi(x,t) phi0(x)>/std(phi0) of the unbiased run, C ~ K^-lambda_bar.
Fits: ln R = ln B + gamma ln K and ln C = c - lambda_bar ln K over K >= K_min after formation, in a box-safe window
(t <= L/4 for PRS; xi_ref <= N/6 for model A) and an extended window (t <= L/2; xi_ref <= N/3).
Writes results/lattice_response.json.  [C] computed.
"""
from __future__ import annotations

import glob
import json
import math
import os
import re

import numpy as np
from scipy.special import erfinv

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "results", "lattice_runs")
OUT = os.path.join(HERE, "..", "results")


def load():
    runs = {}
    for fn in glob.glob(os.path.join(RUNS, "*.json")):
        r = json.load(open(fn))
        c = r["config"]
        m = re.match(r"(.+)_s(\d+)_(m0|mp[0-9.e-]+|mm[0-9.e-]+)$", c["name"])
        runs[(c["group"], c["seed"], c["m0"])] = r
    return runs


def series(r, key):
    return np.array([e["t"] for e in r["series"]]), np.array([np.nan if e[key] is None else e[key] for e in r["series"]])


def formation_index(r):
    for i, e in enumerate(r["series"]):
        if e["bulk_core_frac"] < 0.05 and e["phi_min"] < -0.9 and e["phi_max"] > 0.9:
            return i
    return None


def linfit(x, y):
    if len(x) < 3:
        return None
    A = np.vstack([np.ones_like(x), x]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    s2 = float(np.sum((y - yhat) ** 2) / max(len(x) - 2, 1))
    cov = s2 * np.linalg.inv(A.T @ A)
    return {"intercept": float(coef[0]), "slope": float(coef[1]), "sig_slope_formal": float(math.sqrt(cov[1, 1])),
            "n": int(len(x)), "x_range": [float(np.exp(x.min())), float(np.exp(x.max()))]}


def analyse():
    runs = load()
    groups = sorted(set(k[0] for k in runs))
    out = {}
    for g in groups:
        seeds = sorted(set(k[1] for k in runs if k[0] == g))
        res_seeds = []
        for sd in seeds:
            ref = runs.get((g, sd, 0.0))
            if ref is None:
                continue
            c = ref["config"]
            d, N, dyn = c["d"], c["N"], c["dyn"]
            t, xi = series(ref, "xi")
            _, C = series(ref, "C_s_phi0")
            _, Css = series(ref, "C_s_s0")
            _, mref = series(ref, "m")
            iform = formation_index(ref)
            t_f, xi_f = t[iform], xi[iform]
            K = xi / xi_f
            if dyn == "prs":
                safe = (t >= t_f) & (t <= N / 4 + 1e-9)
                ext = (t >= t_f) & (t <= N / 2 + 1e-9)
            else:
                safe = (t >= t_f) & (xi <= N / 6)
                ext = (t >= t_f) & (xi <= N / 3)
            row = {"seed": sd, "d": d, "N": N, "dyn": dyn, "t_form": float(t_f), "xi_form": float(xi_f),
                   "K_end_safe": float(K[safe].max()) if safe.any() else None,
                   "K_end_ext": float(K[ext].max()) if ext.any() else None,
                   "m_ref_end": float(mref[-1]), "pairs": []}
            # autocorrelation exponent
            for lab, win in (("safe", safe), ("ext", ext)):
                sel = win & (K >= 1.5) & (C > 0)
                f = linfit(np.log(K[sel]), np.log(C[sel]))
                row[f"lambda_bar_{lab}"] = (-f["slope"] if f else None)
                row[f"gamma_from_lambda_{lab}"] = (d + f["slope"] if f else None)
                sel = win & (K >= 1.5) & (Css > 0)
                f2 = linfit(np.log(K[sel]), np.log(Css[sel]))
                row[f"lambda_bar_ss_{lab}"] = (-f2["slope"] if f2 else None)
            m0s = sorted(set(k[2] for k in runs if k[0] == g and k[1] == sd and k[2] > 0))
            for m0 in m0s:
                rp, rm = runs.get((g, sd, m0)), runs.get((g, sd, -m0))
                if rp is None or rm is None:
                    continue
                tp, mp = series(rp, "m")
                tm, mm = series(rm, "m")
                n = min(len(tp), len(tm), len(t))
                Rm = (mp[:n] - mm[:n]) / (2 * m0)
                yp = erfinv(np.clip(mp[:n], -0.999999, 0.999999))
                ym = erfinv(np.clip(mm[:n], -0.999999, 0.999999))
                Ry = (yp - ym) / (2 * erfinv(m0))
                xi_i = float(xi[0])
                pr = {"m0": m0, "R_y_at_formation": float(Ry[iform]), "R_m_at_formation": float(Rm[iform]),
                      "xi_initial": xi_i, "xi_form_over_xi_initial": float(xi_f / xi_i),
                      "gamma_pre_formation = ln R_y(form)/ln(xi_form/xi_i)": float(math.log(Ry[iform]) / math.log(xi_f / xi_i)),
                      "max_abs_m": float(max(np.abs(mp[:n]).max(), np.abs(mm[:n]).max())),
                      "t": [float(x) for x in t[:n]], "K": [float(x) for x in K[:n]],
                      "R_y": [float(x) for x in Ry], "R_m": [float(x) for x in Rm]}
                for lab, win in (("safe", safe[:n]), ("ext", ext[:n])):
                    sel = win & (K[:n] >= 1.5) & (Ry > 0)
                    f = linfit(np.log(K[:n][sel]), np.log(Ry[sel] / Ry[iform]))
                    pr[f"gamma_Ry_{lab}"] = f["slope"] if f else None
                    pr[f"fit_Ry_{lab}"] = f
                    sel = win & (K[:n] >= 1.5) & (Rm > 0)
                    f = linfit(np.log(K[:n][sel]), np.log(Rm[sel] / Rm[iform]))
                    pr[f"gamma_Rm_{lab}"] = f["slope"] if f else None
                row["pairs"].append(pr)
            res_seeds.append(row)
        summ = {"seeds": res_seeds}
        for key in ("lambda_bar_safe", "lambda_bar_ext", "gamma_from_lambda_safe", "gamma_from_lambda_ext",
                    "lambda_bar_ss_safe", "lambda_bar_ss_ext"):
            vals = [r[key] for r in res_seeds if r.get(key) is not None]
            if vals:
                summ[key] = {"mean": float(np.mean(vals)), "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else None,
                             "values": vals}
        m0all = sorted(set(p["m0"] for r in res_seeds for p in r["pairs"]))
        for m0 in m0all:
            for key in ("gamma_Ry_safe", "gamma_Ry_ext", "gamma_Rm_safe", "gamma_Rm_ext", "R_y_at_formation",
                        "gamma_pre_formation = ln R_y(form)/ln(xi_form/xi_i)", "xi_form_over_xi_initial"):
                vals = [p[key] for r in res_seeds for p in r["pairs"] if p["m0"] == m0 and p.get(key) is not None]
                if vals:
                    summ[f"{key}_m0={m0:g}"] = {"mean": float(np.mean(vals)),
                                                 "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else None,
                                                 "values": vals}
        out[g] = summ
    return out


def main():
    out = analyse()
    out["_note"] = ("gamma = d - lambda_bar is the growth exponent of a small systematic magnetization, "
                    "m(K) = m_f K^gamma (K = xi_ref/xi_ref(formation)); R_y is erf-linearized.  Literature: "
                    "model A 2D lambda_bar = 1.25 (gamma 0.75), Liu-Mazenko 3D 1.673 (gamma 1.327), Gaussian/OJK "
                    "lambda_bar = d/2 (gamma = d/2, the maximum allowed by the Fisher-Huse bound).")
    with open(os.path.join(OUT, "lattice_response.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    o = main()
    for g, s in o.items():
        if g.startswith("_"):
            continue
        print(g)
        for r in s["seeds"]:
            print("  seed %d t_form %.1f xi_form %.2f K_end safe/ext %s/%s m_ref_end %.3f lam_bar safe/ext %s/%s (ss %s/%s)" % (
                r["seed"], r["t_form"], r["xi_form"], r["K_end_safe"] and round(r["K_end_safe"], 2), r["K_end_ext"] and round(r["K_end_ext"], 2), r["m_ref_end"],
                r["lambda_bar_safe"] and round(r["lambda_bar_safe"], 3), r["lambda_bar_ext"] and round(r["lambda_bar_ext"], 3),
                r["lambda_bar_ss_safe"] and round(r["lambda_bar_ss_safe"], 3), r["lambda_bar_ss_ext"] and round(r["lambda_bar_ss_ext"], 3)))
            for p in r["pairs"]:
                print("     m0 %g: R_y(form) %.2f gamma_Ry safe/ext %s/%s  gamma_Rm safe/ext %s/%s  max|m| %.2f" % (
                    p["m0"], p["R_y_at_formation"],
                    p["gamma_Ry_safe"] and round(p["gamma_Ry_safe"], 3), p["gamma_Ry_ext"] and round(p["gamma_Ry_ext"], 3),
                    p["gamma_Rm_safe"] and round(p["gamma_Rm_safe"], 3), p["gamma_Rm_ext"] and round(p["gamma_Rm_ext"], 3),
                    p["max_abs_m"]))
        for k, v in s.items():
            if k != "seeds":
                print("   ", k, round(v["mean"], 3), v["std"] and round(v["std"], 3))
