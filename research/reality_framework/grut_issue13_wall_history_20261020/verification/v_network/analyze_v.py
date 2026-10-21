#!/usr/bin/env python3
"""Aggregate the verifier's own runs into results/summary_v.json (no module code used)."""
import json, glob, os
import numpy as np

R = "results/runs"


def load(name):
    p = os.path.join(R, name + ".json")
    return json.load(open(p)) if os.path.exists(p) else None


def window_mean(series, key, t1, t2):
    v = [r[key] for r in series if t1 - 1e-9 <= r["tau"] <= t2 + 1e-9]
    return (float(np.mean(v)) if v else None), len(v)


def slope(series, key, t1, t2):
    t = np.array([r["tau"] for r in series if t1 <= r["tau"] <= t2])
    s = np.array([r[key] for r in series if t1 <= r["tau"] <= t2])
    if len(t) < 3 or np.any(s <= 0):
        return None
    return float(np.polyfit(np.log(t), np.log(s), 1)[0])


out = {"A_definition": "A = rho_w t / sigma = (S_c/V_c) tau / 2 (a ~ tau, t = a tau/2); S_c comoving length (2D) or area (3D)"}

# (a) 2D PRS alpha=2 2048^2
g = {}
for s in (1, 2, 11, 12):
    r = load(f"prs2d_a2_N2048_s{s}")
    if r is None:
        continue
    L = r["L"]
    ser = r["series"]
    d = {"L": L, "tau_end": r["tau_end_reached"], "wall_seconds": r["wall_seconds"]}
    for key in ("A_grad", "A_ms", "A_iso", "A_raw"):
        d[key + "_win_h1/8-1/4"], d["n_win"] = window_mean(ser, key, L / 8, L / 4)
        d[key + "_box_h1/4-1/2"], _ = window_mean(ser, key, L / 4, L / 2)
    d["mu_S_ms_L/32-L/4"] = slope(ser, "S_ms", L / 32, L / 4)
    d["fplus_end"] = ser[-1]["fplus"]
    d["v2_win"], _ = window_mean(ser, "v2", L / 8, L / 4)
    d["A_ms_series"] = [(round(x["tau"], 2), round(x["A_ms"], 4), round(x["fplus"], 3)) for x in ser]
    g[f"seed{s}"] = d
if g:
    for key in ("A_grad_win_h1/8-1/4", "A_ms_win_h1/8-1/4", "A_iso_win_h1/8-1/4", "A_ms_box_h1/4-1/2", "mu_S_ms_L/32-L/4"):
        vals = [g[k][key] for k in g if k.startswith("seed") and g[k].get(key) is not None]
        g["mean_" + key] = float(np.mean(vals))
        g["std_" + key] = float(np.std(vals, ddof=1)) if len(vals) > 1 else None
if g:
    ind = [g[k]["A_ms_win_h1/8-1/4"] for k in ("seed11", "seed12") if k in g]
    same = [g[k]["A_ms_win_h1/8-1/4"] for k in ("seed1", "seed2") if k in g]
    g["note"] = "seeds 1,2 reproduce the module's seed-1,2 realizations (same RNG call); seeds 11,12 are independent"
    if ind:
        g["mean_A_ms_win_independent_seeds"] = float(np.mean(ind))
    allv = ind + same
    g["mean_A_ms_win_all4"] = float(np.mean(allv))
    g["std_A_ms_win_all4"] = float(np.std(allv, ddof=1)) if len(allv) > 1 else None
    # per-output scatter inside the window (all seeds)
    pts = []
    for k in g:
        if k.startswith("seed"):
            L = g[k]["L"]
            pts += [a for (t, a, f) in g[k]["A_ms_series"] if L / 8 <= t <= L / 4]
    g["window_per_output_min_max"] = [float(min(pts)), float(max(pts))]
out["a_prs2d_a2_N2048_L1024w"] = g

# (b) 3D PRS alpha=3 160^3
g = {}
for s in (1, 2):
    r = load(f"prs3d_a3_N160_s{s}")
    if r is None:
        continue
    L = r["L"]
    ser = r["series"]
    d = {"L": L, "tau_end": r["tau_end_reached"], "wall_seconds": r["wall_seconds"]}
    for key in ("A_grad", "A_iso", "A_raw"):
        d[key + "_win_h1/8-1/4"], d["n_win"] = window_mean(ser, key, L / 8, L / 4)
        d[key + "_box_h1/4-1/2"], _ = window_mean(ser, key, L / 4, L / 2)
        d[key + "_tau16-20"], _ = window_mean(ser, key, 16, 20)
    d["mu_S_grad_L/16-L/4"] = slope(ser, "S_grad", L / 16, L / 4)
    d["v2_win"], _ = window_mean(ser, "v2", L / 8, L / 4)
    d["A_grad_series"] = [(round(x["tau"], 2), round(x["A_grad"], 4), round(x["A_iso"], 4)) for x in ser]
    g[f"seed{s}"] = d
if g:
    for key in ("A_grad_win_h1/8-1/4", "A_iso_win_h1/8-1/4", "A_grad_tau16-20", "A_grad_box_h1/4-1/2"):
        vals = [g[k][key] for k in g if k.startswith("seed") and g[k].get(key) is not None]
        g["mean_" + key] = float(np.mean(vals))
        g["std_" + key] = float(np.std(vals, ddof=1)) if len(vals) > 1 else None
out["b_prs3d_a3_N160_L80w"] = g

# (c) biased runs 1024^2
g = {}
ref = load("prs2d_a2_N1024_p050_s1")
for name in ("prs2d_a2_N1024_p035_s1", "prs2d_a2_N1024_p035_s2", "prs2d_a2_N1024_p035_s11", "prs2d_a2_N1024_p045_s1", "prs2d_a2_N1024_p050_s1"):
    r = load(name)
    if r is None:
        continue
    ser = r["series"]
    d = {"series_tau_Ams_fplus": [(round(x["tau"], 2), round(x["A_ms"], 5), round(x["fplus"], 4),
                                   round(x["S_ms"], 2)) for x in ser]}
    # first output with no sign change at all (S_raw == 0)
    t0 = next((x["tau"] for x in ser if x["S_raw"] == 0), None)
    d["tau0_no_sign_change"] = t0
    if ref is not None and name != "prs2d_a2_N1024_p050_s1":
        refA = {round(x["tau"], 6): x["A_ms"] for x in ref["series"]}
        t10 = None
        for x in ser:
            a_ref = refA.get(round(x["tau"], 6))
            if a_ref and x["A_ms"] < 0.1 * a_ref:
                t10 = x["tau"]
                break
        d["tau10_A_below_0.1_ref"] = t10
    d["fplus_end"] = ser[-1]["fplus"]
    d["S_end"] = ser[-1]["S_ms"]
    g[name] = d
out["c_biased"] = g

if os.path.exists("results/bubble_v.json"):
    b = json.load(open("results/bubble_v.json"))
    out["d_bubbles"] = [{k: v for k, v in x.items() if k != "R_of_tau_sample"} for x in b]
if os.path.exists("results/ng_v.json"):
    out["ng"] = json.load(open("results/ng_v.json"))["c_inf"]
json.dump(out, open("results/summary_v.json", "w"), indent=1)
for k, v in out.items():
    if isinstance(v, dict):
        print(k, {kk: vv for kk, vv in v.items() if not isinstance(vv, (dict, list))})
