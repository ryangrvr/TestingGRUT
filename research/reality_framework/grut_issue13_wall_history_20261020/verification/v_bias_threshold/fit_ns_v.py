"""fit_ns_v.py -- independent re-analysis (verifier) of modules/network_sim biased-decay runs (read only).

Written by the verifier without importing or copying modules/bias_threshold code.
Differences from the module on purpose:
  * the unbiased reference is PAIRED per seed (network_sim uses the same white noise for every p at a given seed,
    shifted by a quantile), not the seed mean;
  * tau_10 is located by linear interpolation of log(A/A_ref) in log(tau) on the paired reference;
  * gamma is estimated three ways:
      G1  ln K10 = (ln x10 - ln y_ref)/gamma  (regression of ln K10 on ln y_ref, y = erfinv|m| at the reference time),
          and the inverse regression (errors-in-variables bracket);
      G2  trajectories: per run, slope of ln erfinv|m(K)| vs ln K between the reference time and |m| = 0.97;
          local slopes between successive outputs, binned in K;
      G3  sensitivity of G1 to the choice of reference ("formation") output, and to the length measure
          (xi_ref = tau/(2 A_ref)  vs  conformal time tau).
  * K10 = c m0^-alpha refit (m0 = initial bias), unweighted OLS on seed means, plus seed bootstrap.
Length measures: "xi" (tau/2A, PRS area), "xi_iso" (tau/2A_iso), "nl" (1/raw link count), "tau".
Outputs results/fit_ns_v.json.  numpy/scipy only.
"""
import glob
import json
import math
import os
import re

import numpy as np
from scipy.special import erfinv

HERE = os.path.dirname(os.path.abspath(__file__))
NS = os.path.join(HERE, "..", "..", "modules", "network_sim", "results")
OUT = os.path.join(HERE, "results")

GROUPS = ["prs3d_bias_N128_L128", "prs2d_bias_a2_N1024_L1024", "prs2d_biasdamp_a2_N1024_L1024",
          "prs2d_biaslate_a2_N1024_L1024"]


def load(prefix):
    runs = {}
    for fn in sorted(glob.glob(os.path.join(NS, "runs", prefix + "_p*_s*.json"))):
        mm = re.search(r"_p([0-9.]+)_s(\d+)\.json$", fn)
        r = json.load(open(fn))
        S = r["series"]
        runs[(float(mm.group(1)), int(mm.group(2)))] = {
            "tau": np.array([e["tau"] for e in S]), "A": np.array([e["A"] for e in S]),
            "m": np.array([2 * e["f_plus"] - 1 for e in S]),
            "bulk_lt_half": np.array([e["frac_absphi_lt_half"] for e in S]),
            "nl": np.array([e["nlinks"] for e in S], dtype=float),
            "A_iso": np.array([e["A_iso"] for e in S])}
    return runs


def interp_log(x, xs, ys):
    return float(np.interp(math.log(x), np.log(xs), ys))


def analyse(prefix, tau_form_ns):
    runs = load(prefix)
    seeds = sorted({s for (p, s) in runs})
    ps = sorted({p for (p, s) in runs if p < 0.5 - 1e-9}, reverse=True)
    refs = {s: runs[(0.5, s)] for s in seeds}
    taus = refs[seeds[0]]["tau"]
    # candidate reference ("formation") outputs: network_sim's tau_form and its neighbours
    i0 = int(np.argmin(np.abs(taus - tau_form_ns)))
    ref_idx = [i for i in (i0 - 2, i0 - 1, i0, i0 + 1, i0 + 2) if 0 <= i < len(taus)]
    res = {"tau_form_networksim": tau_form_ns, "reference_outputs": [float(taus[i]) for i in ref_idx], "by_ref": {}}

    def xi(s, t, which="A"):   # paired reference length V/Area (comoving) = tau / (2 A_ref)
        R = refs[s]
        if which == "nl":          # raw sign-changing link count: xi ~ 1/nl
            return 1.0 / math.exp(interp_log(t, R["tau"], np.log(R["nl"])))
        return t / (2 * interp_log(t, R["tau"], R[which]))

    for i in ref_idx:
        t_ref = float(taus[i])
        for measure in ("xi", "tau", "nl", "xi_iso"):
            rows = []
            for p in ps:
                for s in seeds:
                    r = runs.get((p, s))
                    if r is None:
                        continue
                    R = refs[s]
                    n = min(len(r["tau"]), len(R["tau"]))
                    t = r["tau"][:n]
                    ratio = r["A"][:n] / R["A"][:n]
                    mabs = np.abs(r["m"][:n])
                    m_ref = abs(interp_log(t_ref, r["tau"], r["m"]))
                    sel = np.nonzero((t > t_ref + 1e-9) & (ratio < 0.1))[0]
                    if m_ref > 0.995 or interp_log(t_ref, t, ratio) < 0.1 or sel.size == 0:
                        continue
                    j = int(sel[0])
                    lr0, lr1 = math.log(max(ratio[j - 1], 1e-12)), math.log(max(ratio[j], 1e-12))
                    x = (math.log(0.1) - lr0) / (lr1 - lr0)
                    t10 = math.exp(math.log(t[j - 1]) + x * (math.log(t[j]) - math.log(t[j - 1])))
                    if measure in ("xi", "nl", "xi_iso"):
                        wh = {"xi": "A", "nl": "nl", "xi_iso": "A_iso"}[measure]
                        K10 = xi(s, t10, wh) / xi(s, t_ref, wh)
                        Kj = [xi(s, tt, wh) / xi(s, t_ref, wh) for tt in t]
                    else:
                        K10 = t10 / t_ref
                        Kj = list(t / t_ref)
                    m10 = abs(interp_log(t10, r["tau"], r["m"]))
                    # trajectory: points from t_ref (inclusive) until |m| > 0.97
                    traj = [(Kj[k], float(erfinv(mabs[k]))) for k in range(n)
                            if t[k] >= t_ref - 1e-9 and mabs[k] < 0.97 and t[k] <= R["tau"][-1]]
                    rows.append({"p": p, "seed": s, "m0": 2 * p - 1, "m_ref": m_ref, "y_ref": float(erfinv(m_ref)),
                                 "K10": K10, "m_at_tau10": m10, "traj": traj})
            if len(rows) < 3 or len({r["p"] for r in rows}) < 2:
                res["by_ref"][f"{t_ref:.3f}_{measure}"] = None
                continue
            ly = np.log([r["y_ref"] for r in rows])
            lK = np.log([r["K10"] for r in rows])
            ok = lK > 1e-3
            # G1: regress ln K10 on ln y   (slope = -1/gamma)
            b1, b0 = np.polyfit(ly[ok], lK[ok], 1)
            g_direct = -1.0 / b1
            # inverse regression: ln y on ln K10 (slope = -gamma)
            c1, c0 = np.polyfit(lK[ok], ly[ok], 1)
            g_inverse = -c1
            # seed-bootstrap of the direct estimate (resample seeds within each p)
            rng = np.random.default_rng(12345)
            bylab = {}
            for k, r in enumerate(rows):
                if ok[k]:
                    bylab.setdefault(r["p"], []).append(k)
            gb = []
            for _ in range(2000):
                idx = [rng.choice(v) for v in bylab.values() for _ in v]
                if len({rows[k]["p"] for k in idx}) < 2:
                    continue
                bb1, _ = np.polyfit(ly[idx], lK[idx], 1)
                gb.append(-1.0 / bb1)
            # G2: per-run trajectory slopes
            tr_slopes, local = [], []
            for r in rows:
                tr = np.array(r["traj"])
                if len(tr) >= 3 and tr[-1, 0] > 1.2:
                    sl = np.polyfit(np.log(tr[:, 0]), np.log(tr[:, 1]), 1)[0]
                    tr_slopes.append(float(sl))
                for a, b in zip(tr[:-1], tr[1:]):
                    if b[0] > a[0] * 1.02 and a[1] > 0:
                        local.append((math.sqrt(a[0] * b[0]), math.log(b[1] / a[1]) / math.log(b[0] / a[0])))
            local = np.array(local) if local else np.zeros((0, 2))
            bins = [1.0, 1.3, 1.7, 2.3, 3.2, 4.5, 6.5, 10, 15]
            lb = []
            for lo, hi in zip(bins[:-1], bins[1:]):
                q = local[(local[:, 0] >= lo) & (local[:, 0] < hi), 1] if len(local) else []
                if len(q):
                    lb.append({"K_bin": [lo, hi], "n": int(len(q)), "mean_local_gamma": float(np.mean(q)),
                               "sem": float(np.std(q, ddof=1) / math.sqrt(len(q))) if len(q) > 1 else None})
            x10 = math.exp(b0 * g_direct)   # ln K10 = b0 + b1 ln y => ln x10 = -b0/b1 = b0*gamma
            res["by_ref"][f"{t_ref:.3f}_{measure}"] = {
                "t_ref": t_ref, "measure": measure, "n_runs": int(ok.sum()), "n_bias_values": len(bylab),
                "gamma_direct": float(g_direct), "gamma_inverse_regression": float(g_inverse),
                "gamma_seed_bootstrap_16_84": [float(np.percentile(gb, 16)), float(np.percentile(gb, 84))] if gb else None,
                "x10": float(x10),
                "trajectory_slopes_mean": float(np.mean(tr_slopes)) if tr_slopes else None,
                "trajectory_slopes_range": [float(min(tr_slopes)), float(max(tr_slopes))] if tr_slopes else None,
                "local_gamma_vs_K": lb,
                "rows": [{k: v for k, v in r.items() if k != "traj"} for r in rows]}
    # K10 = c m0^-alpha (initial bias), at network_sim's formation output, xi measure
    key = f"{float(taus[i0]):.3f}_xi"
    rr = res["by_ref"][key]["rows"]
    byp = {}
    for r in rr:
        byp.setdefault(abs(r["m0"]), []).append(math.log(r["K10"]))
    m0s = np.array(sorted(byp))
    lKm = np.array([np.mean(byp[m]) for m in m0s])
    a1, a0 = np.polyfit(np.log(m0s), lKm, 1)
    res["K10_vs_initial_m0"] = {"m0": m0s.tolist(), "K10_seedmean_geo": np.exp(lKm).tolist(),
                                "alpha": float(-a1), "c": float(math.exp(a0))}
    if len(m0s) >= 3:
        res["K10_vs_initial_m0"]["alpha_smallest_two"] = float(-(lKm[1] - lKm[0]) / math.log(m0s[1] / m0s[0]))
    return res


def main():
    summ = json.load(open(os.path.join(NS, "summary.json")))
    out = {}
    for g in GROUPS:
        tf = summ["biased"][g]["tau_formation_ref"]
        out[g] = analyse(g, tf)
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, "fit_ns_v.json"), "w"), indent=1)
    for g, r in out.items():
        print("==", g, "tau_form(ns)", r["tau_form_networksim"])
        for k, v in r["by_ref"].items():
            if v is None:
                print("  ", k, "insufficient")
                continue
            print("   ref %-12s n=%2d nb=%d gamma_dir %.3f inv %.3f boot %s x10 %.3f traj %.3f %s" % (
                k, v["n_runs"], v["n_bias_values"], v["gamma_direct"], v["gamma_inverse_regression"],
                None if v["gamma_seed_bootstrap_16_84"] is None else "[%.3f,%.3f]" % tuple(v["gamma_seed_bootstrap_16_84"]),
                v["x10"], v["trajectory_slopes_mean"] or float("nan"),
                v["trajectory_slopes_range"] and "[%.2f,%.2f]" % tuple(v["trajectory_slopes_range"])))
        k0 = f"{r['tau_form_networksim']:.3f}_xi"
        k0 = [k for k in r["by_ref"] if k.endswith("_xi") and abs(float(k.split('_')[0]) - r['tau_form_networksim']) < 0.6]
        for k in k0:
            v = r["by_ref"][k]
            print("   local gamma vs K (ref %s):" % k, [(b["K_bin"], b["n"], round(b["mean_local_gamma"], 3)) for b in v["local_gamma_vs_K"]])
            print("   m at tau10:", sorted({round(x["m_at_tau10"], 3) for x in v["rows"]}))
        print("   K10 vs m0:", r["K10_vs_initial_m0"])


if __name__ == "__main__":
    main()
