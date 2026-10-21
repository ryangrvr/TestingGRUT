"""Aggregate results/runs/*.json (+ validation.json, bubbles.json) into results/summary.json.
Every headline number in NOTES.md is taken from summary.json."""
import glob
import json
import math
import os
import re
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
RUNS = os.path.join(RES, "runs")


def load_runs():
    runs = {}
    for fn in sorted(glob.glob(os.path.join(RUNS, "*.json"))):
        with open(fn) as f:
            r = json.load(f)
        runs[r["config"]["name"]] = r
    return runs


def series(r, key):
    t = np.array([m["tau"] for m in r["series"]])
    y = np.array([np.nan if m.get(key) is None else m[key] for m in r["series"]], dtype=float)
    return t, y


def formation_time(r, thresh=0.05):
    """First output time at which the volume fraction with |phi| < 1/2 outside walls is below thresh:
    fh_bulk = fh - 1.0986 * delta_com * (A_c/V_c), where 1.0986 delta = width of the |tanh| < 1/2 band of a planar
    wall and delta_com = 2 / a^(beta/2) (a = tau / tau_a)."""
    beta = r["config"]["beta"]
    tau_a = r["config"].get("tau_a", 1.0)
    for m in r["series"]:
        dcom = 2.0 / (m["tau"] / tau_a) ** (beta / 2) if beta else 2.0
        bulk = m["frac_absphi_lt_half"] - 2 * math.atanh(0.5) * dcom * m["area_density"]
        if bulk < thresh and m["phi_min"] < -0.9 and m["phi_max"] > 0.9:
            return float(m["tau"])
    return None


def window_stats(r, win=(1 / 8, 1 / 4), fitwin=(1 / 32, 1 / 4), boxwin=(1 / 4, 1 / 2)):
    """Statistics of one run in windows defined by the comoving horizon over the box size, h = tau/L.
    Main window h in [1/8, 1/4] (box >= 4 comoving horizons across); slope fit h in [1/32, 1/4];
    box-affected window h in [1/4, 1/2] reported separately."""
    L = r["config"]["L"]
    t, A = series(r, "A")
    _, dens = series(r, "area_density")
    h = t / L
    out = {"tau_last": float(t[-1]), "h_last": float(h[-1]), "A_end": float(A[-1])}

    def wmean(key, lo, hi):
        _, y = series(r, key)
        m = (h >= lo * (1 - 1e-9)) & (h <= hi * (1 + 1e-9)) & np.isfinite(y)
        return (float(np.mean(y[m])) if m.any() else None), int(m.sum())

    out["A_win"], out["n_win"] = wmean("A", *win)
    tf0 = formation_time(r)
    if tf0 is not None:
        hpf = 1.25 * tf0 / L
        out["A_postform"], out["n_postform"] = wmean("A", hpf, h[-1])
        out["A_postform_window_tau"] = [1.25 * tf0, float(t[-1])]
    out["A_box"], _ = wmean("A", *boxwin)
    tf = formation_time(r)
    out["tau_formation"] = tf
    lo = max(fitwin[0] * L, 2 * tf if tf else 0.0)
    wf = (t >= lo * (1 - 1e-9)) & (h <= fitwin[1] * (1 + 1e-9)) & (dens > 0)
    if wf.sum() >= 3:
        mu, c = np.polyfit(np.log(t[wf]), np.log(dens[wf]), 1)
        out["mu_area_density_slope"] = float(mu)   # scaling: mu = -1
        out["fit_window_tau"] = [float(t[wf][0]), float(t[wf][-1])]
    for key in ("spanning_cluster_area_frac", "largest_cluster_area_frac", "n_domains", "f_plus",
                "n_wall_clusters", "xi_over_tau", "A_iso"):
        out[key + "_win"], _ = wmean(key, *win)
        out[key + "_box"], _ = wmean(key, *boxwin)
        _, y = series(r, key)
        out[key + "_end"] = None if not np.isfinite(y[-1]) else float(y[-1])
    v2, _ = wmean("v2_wall", *win)
    out["v_rms_win"] = float(np.sqrt(v2)) if v2 is not None else None
    _, wc = series(r, "wall_cells_2delta")
    out["wall_cells_2delta_min"] = float(np.nanmin(wc))
    out["stopped_early"] = r["stopped_early"]
    # persistence: wall links present at every output, and box-spanning wall cluster present at h <= 1/4
    _, nl = series(r, "nlinks")
    _, sp = series(r, "n_spanning_clusters")
    out["walls_present_at_all_outputs"] = bool(np.all(nl > 0))
    m4 = h <= 0.25 * (1 + 1e-9)
    out["spanning_wall_present_all_h_le_quarter"] = bool(np.all(sp[m4] > 0))
    out["A_series"] = [[float(x), float(y)] for x, y in zip(t, A)]
    return out


def agg(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if v.size == 0:
        return None
    return {"mean": float(v.mean()), "std": float(v.std(ddof=1)) if v.size > 1 else None,
            "sem": float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else None,
            "n": int(v.size), "values": [float(x) for x in v]}


def compare_resolution(r_lo, r_hi, tau_min, tau_max=None):
    """Relative difference of A(tau) between two runs with the same initial continuum field."""
    t1, A1 = series(r_lo, "A")
    t2, A2 = series(r_hi, "A")
    tmax = min(t1[-1], t2[-1]) if tau_max is None else tau_max
    grid = np.geomspace(tau_min, tmax, 24)
    a1 = np.interp(np.log(grid), np.log(t1), A1)
    a2 = np.interp(np.log(grid), np.log(t2), A2)
    rel = (a1 - a2) / a2
    late = grid >= 0.5 * tmax
    return {"tau_range": [float(tau_min), float(tmax)], "max_abs_rel_diff": float(np.max(np.abs(rel))),
            "late_mean_rel_diff": float(np.mean(rel[late])),
            "A_lo_late": float(np.mean(a1[late])), "A_hi_late": float(np.mean(a2[late])),
            "rel_diff_series": [[float(g), float(x)] for g, x in zip(grid, rel)]}


def biased(runs, prefix, ref_prefix=None):
    """Decay statistics for biased runs named <prefix>_p{p}_s{seed}."""
    groups = defaultdict(dict)
    for name, r in runs.items():
        m = re.match(re.escape(prefix) + r"_p([0-9.]+)_s(\d+)$", name)
        if m:
            groups[float(m.group(1))][int(m.group(2))] = r
    if not groups:
        return None
    # reference unbiased A(tau): mean over p = 0.5 seeds
    ref = groups.get(0.5, {})
    tref = None
    if ref:
        ts = [series(r, "A") for r in ref.values()]
        tmin = max(t[0] for t, _ in ts)
        tmax = min(t[-1] for t, _ in ts)
        tref = np.geomspace(tmin, tmax, 64)
        Aref = np.mean([np.interp(np.log(tref), np.log(t), A) for t, A in ts], axis=0)
    tfs = [formation_time(r) for r in ref.values()]
    tfs = [x for x in tfs if x is not None]
    tform_ref = float(np.median(tfs)) if tfs else None
    out = {"tau_formation_ref": tform_ref}
    for p in sorted(groups, reverse=True):
        rows = []
        for s, r in sorted(groups[p].items()):
            t, A = series(r, "A")
            _, nl = series(r, "nlinks")
            _, fp = series(r, "f_plus")
            row = {"seed": s, "tau_last": float(t[-1]), "stopped_early": r["stopped_early"],
                   "f_plus_initial": float(fp[0]), "f_plus_end": float(fp[-1])}
            # tau_0: first output time with no wall links (network completely gone)
            z = np.nonzero(nl == 0)[0]
            row["tau_zero"] = float(t[z[0]]) if z.size else None
            # formation time: median formation time of the unbiased (p = 0.5) runs of the same group
            tau_form = tform_ref if tform_ref is not None else float(t[0])
            k = int(np.argmin(np.abs(t - tau_form)))
            row["network_gone_before_formation"] = bool(row["tau_zero"] is not None and row["tau_zero"] <= tau_form)
            row["tau_formation"] = tau_form
            # tau_10: first time A < 0.1 * A_ref(tau) after formation
            if tref is not None:
                ok = (t >= tau_form) & (t <= tref[-1])
                ar = np.interp(np.log(np.clip(t, tref[0], tref[-1])), np.log(tref), Aref)
                hit = np.nonzero(ok & (A < 0.1 * ar))[0]
                row["tau_10"] = float(t[hit[0]]) if hit.size else None
                hit5 = np.nonzero(ok & (A < 0.5 * ar))[0]
                row["tau_50"] = float(t[hit5[0]]) if hit5.size else None
                # coarsening factor of the unbiased reference network between formation and tau_10:
                # K = xi_ref(tau_10) / xi_ref(tau_form), xi = V / area = tau / (2 A)
                if row["tau_10"] is not None and tau_form >= tref[0]:
                    def xi(tt):
                        return tt / (2 * np.interp(np.log(tt), np.log(tref), Aref))
                    row["coarsening_factor_ref_at_tau_10"] = float(xi(min(row["tau_10"], tref[-1])) / xi(tau_form))
            # minority-phase percolation right after formation
            mser = r["series"][k]
            minority = "plus" if p < 0.5 else "minus"
            row["minority_spans_at_formation"] = mser.get(f"largest_{minority}_spans")
            row["tau_formation_sample"] = float(t[k])
            row["minority_spans_initial"] = r["series"][0].get(f"largest_{minority}_spans")
            rows.append(row)
        def coll(key):
            vals = [x[key] for x in rows if x.get(key) is not None]
            return agg(vals) if vals else None
        out[f"{p:.2f}"] = {"runs": rows,
                           "n_seeds_no_network_formed": int(sum(1 for x in rows if x["network_gone_before_formation"])), "tau_zero": coll("tau_zero"), "tau_10": coll("tau_10"),
                           "tau_50": coll("tau_50"), "tau_formation": coll("tau_formation"),
                           "coarsening_factor_ref_at_tau_10": coll("coarsening_factor_ref_at_tau_10"),
                           "minority_spans_after_formation": [x["minority_spans_at_formation"] for x in rows],
                           "n_seeds_gone_by_end": int(sum(1 for x in rows if x["tau_zero"] is not None)),
                           "n_seeds": len(rows)}
    # power-law fit tau_10 ~ (0.5 - p)^(-q) over p with all seeds decayed
    def usable(p, v, key):
        """all seeds decayed, no seed lost its walls before formation, decay later than formation"""
        return (p[0].isdigit() and float(p) < 0.5 and v[key] is not None and v[key]["n"] == v["n_seeds"]
                and v["n_seeds_no_network_formed"] == 0
                and (tform_ref is None or min(v[key]["values"]) > tform_ref * (1 + 1e-9)))
    pts = [(0.5 - float(p), v["tau_10"]["mean"]) for p, v in out.items() if usable(p, v, "tau_10")]
    if len(pts) >= 3:
        x = np.log([q for q, _ in pts])
        y = np.log([tv for _, tv in pts])
        q, c = np.polyfit(x, y, 1)
        out["fit_tau10_vs_bias"] = {"law": "tau_10 = C (0.5 - p)^q", "q": float(q), "C": float(math.exp(c)),
                                    "points": pts}
    ptsK = [(0.5 - float(p), v["coarsening_factor_ref_at_tau_10"]["mean"]) for p, v in out.items()
            if usable(p, v, "tau_10") and v["coarsening_factor_ref_at_tau_10"] is not None
            and v["coarsening_factor_ref_at_tau_10"]["n"] == v["n_seeds"]]
    if len(ptsK) >= 3:
        x = np.log([q for q, _ in ptsK])
        y = np.log([kv for _, kv in ptsK])
        q, c = np.polyfit(x, y, 1)
        out["fit_K10_vs_bias"] = {"law": "K10 = xi_ref(tau_10)/xi_ref(tau_form) = C (0.5 - p)^q", "q": float(q),
                                  "C": float(math.exp(c)), "points": ptsK,
                                  "note": "fit over tested p only; extrapolation outside 0.05 <= 0.5 - p <= 0.3 not done"}
    pts0 = [(0.5 - float(p), v["tau_zero"]["mean"]) for p, v in out.items() if usable(p, v, "tau_zero")]
    if len(pts0) >= 3:
        x = np.log([q for q, _ in pts0])
        y = np.log([tv for _, tv in pts0])
        q, c = np.polyfit(x, y, 1)
        out["fit_tau0_vs_bias"] = {"law": "tau_zero = C (0.5 - p)^q", "q": float(q), "C": float(math.exp(c)),
                                   "points": pts0}
    return out


def write_headline(S):
    """Compact machine-readable headline (all values copied from summary.json)."""
    G = S["groups"]

    def ms(g, k):
        d = G[g][k]
        return None if d is None else {"mean": d["mean"], "seed_std": d["std"], "n_seeds": d["n"]}
    H = {
        "module": "network_sim",
        "scope": "classical lattice simulation of the X = D = 0 sector of the unchanged Candidate 1 potential, "
                 "radiation era a ~ tau, no backreaction",
        "A_definition": "A = rho_w t / sigma = (A_c/V_c) tau / 2 (Hiramatsu-Kawasaki-Saikawa)",
        "A_3D_PRS_256_box_safe": ms("prs3d_N256_L256", "A_win"),
        "A_3D_PRS_384_box_safe": ms("prs3d_N384_L384", "A_win"),
        "A_3D_PRS_256_box_affected": ms("prs3d_N256_L256", "A_box"),
        "A_3D_phys_256_postformation_box_affected": ms("phys3d_N256_L24", "A_postform"),
        "A_2D_phys_2048_box_safe": ms("phys2d_N2048_L64", "A_win"),
        "A_2D_PRS_alpha2_2048_box_safe": ms("prs2d_a2_N2048_L2048", "A_win"),
        "A_2D_PRS_alpha2_4096_L4096_box_safe": ms("prs2d_a2_N4096_L4096", "A_win"),
        "A_2D_PRS_alpha3_2048_box_safe": ms("prs2d_a3_N2048_L2048", "A_win"),
        "mu_3D_PRS_256": ms("prs3d_N256_L256", "mu_area_density_slope"),
        "mu_3D_PRS_384": ms("prs3d_N384_L384", "mu_area_density_slope"),
        "mu_2D_PRS_alpha2": ms("prs2d_a2_N2048_L2048", "mu_area_density_slope"),
        "mu_2D_PRS_alpha3": ms("prs2d_a3_N2048_L2048", "mu_area_density_slope"),
        "mu_2D_phys": ms("phys2d_N2048_L64", "mu_area_density_slope"),
        "spanning_area_frac_3D_PRS": ms("prs3d_N256_L256", "spanning_cluster_area_frac_win"),
        "spanning_area_frac_2D_phys": ms("phys2d_N2048_L64", "spanning_cluster_area_frac_win"),
        "spanning_area_frac_2D_PRS_alpha2": ms("prs2d_a2_N2048_L2048", "spanning_cluster_area_frac_win"),
        "v_rms_3D_PRS": ms("prs3d_N256_L256", "v_rms_win"),
        "xi_over_tau_3D_PRS": ms("prs3d_N256_L256", "xi_over_tau_win"),
        "resolution_checks": {k: {"late_mean_rel_diff": v["late_mean_rel_diff"], "max_abs_rel_diff": v["max_abs_rel_diff"]}
                              for k, v in S["resolution_checks"].items()},
        "ic_sensitivity_late_start_rel_diff": S.get("ic_sensitivity_late_start", {}).get("late_mean_rel_diff"),
        "biased_fits": {pref: {k: {"q": v["q"], "C": v["C"], "law": v["law"]} for k, v in b.items() if k.startswith("fit")}
                        for pref, b in S["biased"].items()},
        "biased_K10": {pref: {p: (v["coarsening_factor_ref_at_tau_10"]["mean"]
                                  if v.get("coarsening_factor_ref_at_tau_10") else None)
                              for p, v in b.items() if p[0].isdigit()} for pref, b in S["biased"].items()},
        "biased_tau10": {pref: {p: (v["tau_10"]["mean"] if v.get("tau_10") else None)
                                for p, v in b.items() if p[0].isdigit()} for pref, b in S["biased"].items()},
        "biased_failed_setup": "prs2d_biaslate_a2_N1024_L1024 (undamped post-quench oscillations; not wall decay)",
    }
    if "bubbles" in S:
        B = S["bubbles"]
        H["bubble_NG_flat_collapse_over_R0"] = B["ng_flat_collapse_over_R0"]
        H["bubble_NG_superhorizon"] = B["ng_superhorizon"]
        H["bubble_field_vs_NG_tau10_max_abs_rel_dev"] = {
            k: S["bubble_summary"][k]["max_abs_rel_dev_tau10_field_vs_NG"] for k in ("radial_physical", "radial_prs")}
    sc = os.path.join(RES, "scales.json")
    if os.path.exists(sc):
        with open(sc) as f:
            H["scales"] = json.load(f)
    with open(os.path.join(RES, "headline.json"), "w") as f:
        json.dump(H, f, indent=1)


def main():
    runs = load_runs()
    S = {"n_runs": len(runs), "groups": {}}
    groups = defaultdict(list)
    for name in runs:
        if "_bias" in name:
            continue
        g = re.sub(r"_s\d+$", "", name)
        groups[g].append(name)
    for g, names in sorted(groups.items()):
        st = {n: window_stats(runs[n]) for n in sorted(names)}
        keys = ["A_win", "A_box", "A_postform", "A_end", "mu_area_density_slope", "tau_formation", "spanning_cluster_area_frac_win",
                "largest_cluster_area_frac_win", "spanning_cluster_area_frac_box", "v_rms_win", "n_domains_win",
                "n_domains_end", "xi_over_tau_win", "A_iso_win", "f_plus_end"]
        S["groups"][g] = {"per_seed": st, **{k: agg([s.get(k) for s in st.values()]) for k in keys},
                          "config": {k: runs[names[0]]["config"][k] for k in
                                     ("d", "N", "L", "alpha", "beta", "tau_i", "tau_end", "s_c", "amp", "courant")},
                          "dx": runs[names[0]]["dx"]}
    # resolution checks (same continuum initial field)
    rc = {}
    pairs = [("phys2d_N1024_L64_s1", "phys2d_N2048_L64_s1", 4.0),
             ("phys2d_N2048_L64_s1", "phys2d_N4096_L64_s1", 4.0),
             ("prs2d_a2_N2048_L2048_s1", "prs2d_a2_N4096_L2048_s1", 8.0),
             ("prs3d_bias_N128_L128_p0.50_s1", "prs3d_N256_L128_s1", 8.0),
             ("phys3d_N128_L24_s1", "phys3d_N256_L24_s1", 4.0)]
    for lo, hi, tmin in pairs:
        if lo in runs and runs.get(hi):
            rc[f"{lo}__vs__{hi}"] = compare_resolution(runs[lo], runs[hi], tmin)
    S["resolution_checks"] = rc
    # IC sensitivity: late start
    if "prs2d_a2_N2048_L2048_late_s1" in runs and "prs2d_a2_N2048_L2048_s1" in runs:
        S["ic_sensitivity_late_start"] = compare_resolution(runs["prs2d_a2_N2048_L2048_late_s1"],
                                                            runs["prs2d_a2_N2048_L2048_s1"], 40.0, 512.0)
        S["ic_sensitivity_late_start"]["note"] = "tau_i = 20 vs tau_i = 2, same seed and grid, h = tau/L <= 1/4"
    # biased
    S["biased"] = {}
    for pref in ("phys2d_bias_N1024_L32", "prs2d_bias_a2_N1024_L1024", "prs2d_biaslate_a2_N1024_L1024",
                 "prs2d_biasdamp_a2_N1024_L1024", "prs3d_bias_N128_L128"):
        b = biased(runs, pref)
        if b:
            S["biased"][pref] = b
    for extra in ("validation.json", "bubbles.json"):
        fn = os.path.join(RES, extra)
        if os.path.exists(fn):
            with open(fn) as f:
                S[extra.replace(".json", "")] = json.load(f)
    if "bubbles" in S:
        B = S["bubbles"]
        bs = {}
        for key in ("radial_physical", "radial_prs"):
            rows = [r for r in B[key] if not r.get("resolution_check")]
            rel = [r["tau_10"] / r["tau_10_NG"] - 1 for r in rows if r["tau_10"] is not None]
            bs[key] = {
                "n_runs": len(rows),
                "max_abs_rel_dev_tau10_field_vs_NG": float(max(abs(x) for x in rel)),
                "mean_rel_dev_tau10_field_vs_NG": float(np.mean(rel)),
                "all_collapsed": bool(all(r["tau_10"] is not None for r in rows)),
                "numerical_failures_after_collapse": [
                    {"kappa": r["kappa"], "alpha": r["alpha"], "x": r["x"], "tau_i": r["tau_i"],
                     "tau_fail": r["numerical_failure_tau"], "tau_c_NG": r["tau_c_NG"]}
                    for r in rows if r["numerical_failure_tau"] is not None],
                "debris_valid_runs": [
                    {"kappa": r["kappa"], "alpha": r["alpha"], "x": r["x"], "tau_i": r["tau_i"],
                     "debris_max_R_over_R0": r["debris_max_R_over_R0"],
                     "debris_lifetime_over_R0": r.get("debris_lifetime_over_R0"),
                     "negative_cells_at_end": r["negative_cells_at_end"]}
                    for r in rows if r["numerical_failure_tau"] is None],
                "frozen_R_over_R0_at_tau_eq_R0": [
                    {"kappa": r["kappa"], "alpha": r["alpha"], "x": r["x"], "tau_i": r["tau_i"],
                     "value": r.get("R_over_R0_at_tau_R0")} for r in rows if r.get("R_over_R0_at_tau_R0") is not None]}
            # all debris-valid runs end without negative cells?
            bs[key]["debris_valid_all_vanished_by_end"] = bool(all(not d["negative_cells_at_end"]
                                                                   for d in bs[key]["debris_valid_runs"]))
        res_rows = [r for r in B["radial_physical"] if r.get("resolution_check")]
        bs["radial_resolution_check"] = [{"kappa": r["kappa"], "cells_per_delta_end": r["wall_cells_per_delta_end"],
                                          "tau_10": r["tau_10"], "tau_10_NG": r["tau_10_NG"]} for r in res_rows]
        bs["cartesian_vs_radial"] = [{"d": c["d"], "R0": c["R0"], "tau_10_cart": c["tau_10"],
                                      "tau_10_radial": c["tau_10_radial"], "tau_10_NG": c["tau_10_NG"],
                                      "rel_cart_vs_radial": c["tau_10"] / c["tau_10_radial"] - 1,
                                      "cells_per_delta_end": c["wall_cells_per_delta_end"],
                                      "tau_last_negative": c["tau_last_negative"]} for c in B["cartesian_checks"]]
        S["bubble_summary"] = bs
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(S, f, indent=1)
    write_headline(S)
    # console table
    for g, v in S["groups"].items():
        A = v["A_win"]
        Ab = v["A_box"]
        mu = v["mu_area_density_slope"]
        sp = v["spanning_cluster_area_frac_win"]
        f = lambda d, k="mean": (d[k] if d and d.get(k) is not None else float("nan"))
        print(f"{g:32s} A_win={f(A):.3f}+-{f(A,'std'):.3f} (n={f(A,'n')}) A_box={f(Ab):.3f} A_pf={f(v['A_postform']):.3f} "
              f"mu={f(mu):.3f}+-{f(mu,'std'):.3f} span={f(sp):.3f} v_rms={f(v['v_rms_win']):.3f} "
              f"nd={f(v['n_domains_win']):.1f}")
    for k, v in rc.items():
        print("RES", k, f"max|rel|={v['max_abs_rel_diff']:.3f} late_rel={v['late_mean_rel_diff']:.3f}")
    for pref, b in S["biased"].items():
        for p, v in b.items():
            if p.startswith("fit"):
                print("BIAS", pref, p, v["q"], v["C"])
                continue
            if not p[0].isdigit():
                print("BIAS", pref, p, v)
                continue
            print("BIAS", pref, p, "tau10", v["tau_10"] and round(v["tau_10"]["mean"], 2),
                  "tau0", v["tau_zero"] and round(v["tau_zero"]["mean"], 2), "gone", v["n_seeds_gone_by_end"],
                  "/", v["n_seeds"], "noform", v["n_seeds_no_network_formed"],
                  "K10", v["coarsening_factor_ref_at_tau_10"] and round(v["coarsening_factor_ref_at_tau_10"]["mean"], 2))


if __name__ == "__main__":
    main()
