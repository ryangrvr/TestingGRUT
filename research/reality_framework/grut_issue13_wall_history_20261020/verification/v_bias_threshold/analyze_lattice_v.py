"""analyze_lattice_v.py -- gamma (small-bias response) and lambda_bar (unbiased autocorrelation) from lattice_v runs.

Per seed: runs mu = 0, +mu, -mu from one noise realization.
  K(tau)   = nl(tau_form)/nl(tau) of the unbiased run (wall area ~ number of sign-changing links; xi = V/area ~ 1/nl)
  R(tau)   = (erfinv m_+ - erfinv m_-)/2               (paired, erf-linearized response)
  C(tau)   = <sign phi(tau) phi0>/std(phi0)            (unbiased run)
  gamma    = d ln R / d ln K,  lambda_bar = -d ln C / d ln K   (least squares over K >= K_min within a window)
  formation: first tau where the unbiased run's fraction of cells with |phi| < 1/2 falls below f_form
             (interpolated in log tau); f_form = 0.20 central, 0.30 and 0.12 as variants.
Windows: box-safe tau <= N/4; extended tau <= N/2.
Writes results/lattice_v_analysis.json.
"""
import glob
import json
import math
import os
import re

import numpy as np
from scipy.special import erfinv

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "results", "lattice_v")


def load():
    runs = {}
    for fn in glob.glob(os.path.join(D, "lattice_v_N*_s*_mu*.json")):
        mm = re.search(r"_N(\d+)_s(\d+)_mu([+-][0-9.]+)\.json$", fn)
        r = json.load(open(fn))
        runs[(int(mm.group(1)), int(mm.group(2)), float(mm.group(3)))] = r
    return runs


def arr(r, k):
    return np.array([e[k] for e in r["series"]], dtype=float)


def slope(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    c, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    n = len(x)
    if n > 2:
        s2 = float(np.sum((y - A @ c) ** 2) / (n - 2))
        cov = s2 * np.linalg.inv(A.T @ A)
        return float(c[1]), float(math.sqrt(cov[1, 1]))
    return float(c[1]), None


def analyse_seed(r0, rp, rm, f_form, K_min=1.5):
    N = r0["N"]
    tau = arr(r0, "tau")
    nl = arr(r0, "nl")
    frac = arr(r0, "frac_lt_half")
    C = arr(r0, "C")
    mp, mm_ = arr(rp, "m"), arr(rm, "m")
    R = (erfinv(np.clip(mp, -0.999999, 0.999999)) - erfinv(np.clip(mm_, -0.999999, 0.999999))) / 2
    m0 = math.erf(rp["mu"] / (math.sqrt(2) * 0.1))
    y0 = float(erfinv(m0))
    i = int(np.nonzero(frac < f_form)[0][0])
    lt = math.log(tau[i - 1]) + (math.log(f_form) - math.log(frac[i - 1])) / (math.log(frac[i]) - math.log(frac[i - 1])) * (
        math.log(tau[i]) - math.log(tau[i - 1]))
    tf = math.exp(lt)
    nlf = math.exp(np.interp(lt, np.log(tau), np.log(nl)))
    Rf = math.exp(np.interp(lt, np.log(tau), np.log(np.abs(R))))
    K = nlf / nl
    out = {"N": N, "seed": r0["seed"], "tau_form": tf, "f_form": f_form, "quench_amplification_R_form_over_y0": Rf / y0,
           "m_plus_end": float(mp[-1]), "m_unbiased_end": float(arr(r0, "m")[-1])}
    for name, tmax in (("safe", N / 4), ("ext", N / 2)):
        sel = (tau > tf) & (tau <= tmax + 1e-9) & (K >= K_min)
        out[f"K_max_{name}"] = float(K[(tau <= tmax + 1e-9)].max())
        if sel.sum() >= 3:
            g, sg = slope(np.log(K[sel]), np.log(np.abs(R[sel])))
            lb, slb = slope(np.log(K[sel]), np.log(np.abs(C[sel])))
            out[f"gamma_{name}"] = g
            out[f"gamma_{name}_fit_err"] = sg
            out[f"lambda_bar_{name}"] = -lb
            out[f"lambda_bar_{name}_fit_err"] = slb
            out[f"n_pts_{name}"] = int(sel.sum())
            # same quantities against conformal time instead of K (length-measure systematic)
            g_t, _ = slope(np.log(tau[sel] / tf), np.log(np.abs(R[sel])))
            l_t, _ = slope(np.log(tau[sel] / tf), np.log(np.abs(C[sel])))
            out[f"gamma_{name}_vs_tau"] = g_t
            out[f"lambda_bar_{name}_vs_tau"] = -l_t
    # pre-formation (linear stage): R vs K over tau < tf (K from sign-change count of the Gaussian field)
    pre = (tau < tf) & (tau > tau[0])
    if pre.sum() >= 3:
        out["gamma_pre_formation"] = slope(np.log(K[pre]), np.log(np.abs(R[pre])))[0]
    # local exponents between successive outputs after formation
    loc = []
    for a in range(len(tau) - 1):
        if tau[a] >= tf and tau[a + 1] <= N / 2 + 1e-9 and K[a + 1] > K[a] * 1.03:
            loc.append({"K_mid": float(math.sqrt(K[a] * K[a + 1])), "tau_mid": float(math.sqrt(tau[a] * tau[a + 1])),
                        "gamma_loc": float(math.log(abs(R[a + 1] / R[a])) / math.log(K[a + 1] / K[a])),
                        "lambda_bar_loc": float(-math.log(abs(C[a + 1] / C[a])) / math.log(K[a + 1] / K[a]))})
    out["local"] = loc
    out["series"] = {"tau": tau.tolist(), "K": K.tolist(), "R": R.tolist(), "C": C.tolist(), "frac": frac.tolist()}
    return out


def main():
    runs = load()
    keys = sorted({(N, s) for (N, s, mu) in runs})
    res = {"per_seed": {}, "aggregate": {}}
    for f_form in (0.20, 0.30, 0.12):
        per = []
        for (N, s) in keys:
            mus = sorted(mu for (n2, s2, mu) in runs if n2 == N and s2 == s)
            if len(mus) < 3:
                continue
            r0, rm, rp = runs[(N, s, 0.0)], runs[(N, s, mus[0])], runs[(N, s, mus[-1])]
            a = analyse_seed(r0, rp, rm, f_form)
            per.append(a)
            res["per_seed"][f"N{N}_s{s}_f{f_form}"] = {k: v for k, v in a.items() if k != "series"}
            if f_form == 0.20:
                res["per_seed"][f"N{N}_s{s}_f{f_form}"]["series"] = a["series"]
        agg = {}
        for k in ("gamma_safe", "gamma_ext", "lambda_bar_safe", "lambda_bar_ext", "gamma_safe_vs_tau",
                  "lambda_bar_safe_vs_tau", "gamma_ext_vs_tau", "lambda_bar_ext_vs_tau", "gamma_pre_formation",
                  "quench_amplification_R_form_over_y0", "K_max_safe", "K_max_ext"):
            v = [p[k] for p in per if k in p]
            if v:
                agg[k] = {"mean": float(np.mean(v)), "sem": float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else None,
                          "values": [round(x, 4) for x in v]}
        res["aggregate"][f"f_form={f_form}"] = agg
    json.dump(res, open(os.path.join(HERE, "results", "lattice_v_analysis.json"), "w"), indent=1)
    for k, v in res["aggregate"].items():
        print("==", k)
        for kk, vv in v.items():
            print("   %-38s mean %.3f sem %s  %s" % (kk, vv["mean"], None if vv["sem"] is None else "%.3f" % vv["sem"], vv["values"]))
    for k, v in res["per_seed"].items():
        if k.endswith("f0.2"):
            print(k, "tau_form %.2f" % v["tau_form"], "m_unb_end %.3f" % v["m_unbiased_end"])
            print("   local:", [(round(q["K_mid"], 2), round(q["gamma_loc"], 2), round(q["lambda_bar_loc"], 2)) for q in v["local"]])


if __name__ == "__main__":
    main()
