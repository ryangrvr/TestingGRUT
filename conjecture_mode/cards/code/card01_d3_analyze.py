#!/usr/bin/env python3
"""CARD-01 v1 mechanical analysis: applies the owner-locked S1-S5 logic (card01_run_config.decide) to the
D3 minimization outputs.  No threshold, dataset or logic is chosen here.
usage: card01_d3_analyze.py OUTDIR SUMMARY_JSON
"""
import glob, json, os, sys
import numpy as np
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_run_config as RC  # noqa: E402

OUT, SUMMARY = sys.argv[1], sys.argv[2]


def load(combo):
    rows = [json.load(open(f)) for f in glob.glob(os.path.join(OUT, f"{combo}__*.json"))]
    grid, free, cpl = {}, [], []
    for r in rows:
        if r["model"] == "card":
            grid.setdefault(r["eps_input"], []).append(r)
        elif r["model"] == "free":
            free.append(r)
        elif r["model"] == "cpl":
            cpl.append(r)
    return grid, free, cpl


def profile_set(eps, chi, cmin, thr):
    """eps sorted ascending; return [lo, hi] of {chi - cmin <= thr} by linear interpolation."""
    ok = chi - cmin <= thr
    idx = np.where(ok)[0]
    lo_i, hi_i = idx.min(), idx.max()

    def cross(i, j):
        d1, d2 = chi[i] - cmin - thr, chi[j] - cmin - thr
        return eps[i] + (eps[j] - eps[i]) * (0 - d1) / (d2 - d1)
    lo = eps[0] if lo_i == 0 else cross(lo_i - 1, lo_i)
    hi = eps[-1] if hi_i == len(eps) - 1 else cross(hi_i, hi_i + 1)
    return float(lo), float(hi), bool(lo_i == 0), bool(hi_i == len(eps) - 1)


def analyze(combo, branch):
    grid, free, cpl = load(combo)
    sign = -1 if branch == "primary" else +1
    pts = sorted(e for e in grid if (e <= 0 if branch == "primary" else e >= 0))
    table, flags = [], []
    for e in pts:
        c = [r["chi2_eff"] for r in grid[e]]
        best = min(grid[e], key=lambda r: r["chi2_eff"])
        spread = max(c) - min(c)
        if spread > RC.CHI2_START_DISAGREEMENT_FLAG:
            flags.append(e)
        table.append({"eps": e, "chi2": min(c), "start_chi2": sorted(c), "n_starts": len(c), "start_spread": spread,
                      "flag_spread_gt_0p2": spread > RC.CHI2_START_DISAGREEMENT_FLAG,
                      "H0": best["point"].get("H0"), "omch2": best["point"].get("omch2")})
    if 0.0 not in grid:
        return {"combo": combo, "status": "INCOMPLETE (no eps=0)"}
    chi_L = min(r["chi2_eff"] for r in grid[0.0])
    eps = np.array([t["eps"] for t in table]); chi = np.array([t["chi2"] for t in table])
    i = int(np.argmin(chi)); eps_hat, chi_hat = float(eps[i]), float(chi[i])
    free_best = None
    if branch == "primary" and free:
        fb = min(free, key=lambda r: r["chi2_eff"])
        free_best = {"eps": fb["point"]["card01_eps"], "chi2": fb["chi2_eff"],
                     "start_chi2": sorted(r["chi2_eff"] for r in free),
                     "start_eps": [r["point"]["card01_eps"] for r in sorted(free, key=lambda r: r["chi2_eff"])],
                     "spread": max(r["chi2_eff"] for r in free) - min(r["chi2_eff"] for r in free), "n": len(free)}
        if fb["chi2_eff"] < chi_hat:
            eps_hat, chi_hat = float(fb["point"]["card01_eps"]), float(fb["chi2_eff"])
    order = np.argsort(eps)
    lo, hi, lo_edge, hi_edge = profile_set(eps[order], chi[order], chi_hat, RC.PROFILE_SET_DCHI2)
    res = {"combo": combo, "branch": branch, "profile": table, "start_disagreement_flags": flags,
           "free_eps_minimization": free_best, "chi2_Lambda": chi_L, "eps_hat": eps_hat, "chi2_card_best": chi_hat,
           "I_card": chi_L - chi_hat, "profile_set_95": [lo, hi], "profile_set_hits_grid_edge": [lo_edge, hi_edge]}
    q = max(chi_L - chi_hat, 0.0)
    res["boundary_null_chernoff_p"] = float(0.5 * chi2dist.sf(q, 1)) if q > 0 else 0.5
    if cpl:
        cb = min(cpl, key=lambda r: r["chi2_eff"])
        res["chi2_w0wa"] = cb["chi2_eff"]
        res["w0wa_best"] = {"w0": cb["point"].get("w"), "wa": cb["point"].get("wa"),
                            "start_chi2": sorted(r["chi2_eff"] for r in cpl),
                            "spread": max(r["chi2_eff"] for r in cpl) - min(r["chi2_eff"] for r in cpl), "n": len(cpl)}
        res["I_CPL"] = chi_L - cb["chi2_eff"]
        res["F"] = res["I_card"] / res["I_CPL"] if res["I_CPL"] > 0 else None
        res["tension_branch_vs_CPL_dchi2"] = chi_hat - cb["chi2_eff"]
    if branch == "primary" and "I_CPL" in res:
        at_bound = eps_hat <= RC.EPS_RANGE_BOUNDARY + 1e-3
        mech = RC.decide(res["I_CPL"], res["I_card"], eps_hat, lo, hi, at_bound)
        # Owner convergence ruling: load-bearing minimizations = LambdaCDM null, CPL comparator, branch best fit
        # (best grid point and free-eps fit), and the grid points bracketing each 95% profile-set crossing.
        F = RC.CHI2_START_DISAGREEMENT_FLAG
        sp = {t["eps"]: t for t in table}
        srt = sorted(sp)
        lb = {"LambdaCDM_null(eps=0)": sp[0.0]["start_spread"],
              "CPL_comparator": res["w0wa_best"]["spread"],
              f"branch_best_grid(eps={float(eps[i])})": sp[float(eps[i])]["start_spread"]}
        if free_best:
            lb["free_eps_fit"] = free_best["spread"]
        for edge, val, hit in (("lo", lo, lo_edge), ("hi", hi, hi_edge)):
            if hit:
                continue
            below = max([e for e in srt if e <= val], default=None); above = min([e for e in srt if e >= val], default=None)
            for e in (below, above):
                if e is not None:
                    lb[f"profile_set_{edge}_bracket(eps={e})"] = sp[e]["start_spread"]
        incomplete = [k for k, t in sp.items() if t["n_starts"] < RC.N_STARTS]
        if free_best and free_best["n"] < RC.N_STARTS:
            incomplete.append("free")
        if res["w0wa_best"]["n"] < RC.N_STARTS:
            incomplete.append("cpl")
        flagged = {k: v for k, v in lb.items() if v > F}
        res["load_bearing_start_spreads"] = lb
        res["load_bearing_flagged_gt_0p2"] = flagged
        res["incomplete_points"] = incomplete
        res["mechanical_S5_output"] = mech
        if incomplete:
            res["result_state"] = "INCOMPLETE - NO VERDICT"
        elif flagged:
            res["result_state"] = "CARD-01-v1 - NUMERICALLY UNRESOLVED / NO ACCEPTED SCIENTIFIC VERDICT"
            res["mechanical_S5_output_label"] = ("MECHANICAL OUTPUT - NOT SCIENTIFICALLY ACCEPTED DUE TO FROZEN "
                                                 "CONVERGENCE FAILURE")
        else:
            res["result_state"] = mech
    elif branch == "control":
        res["result_state"] = "CONTROL - NOT A GRUT CLAIM (no result state)"
    return res


def main():
    out = {"primary": analyze(RC.PRIMARY, "primary")}
    for c in ["EXT_PANTHEONPLUS", "EXT_UNION3", "EXT_DESY5"]:
        if glob.glob(os.path.join(OUT, f"{c}__*.json")):
            out[c] = analyze(c, "primary")
    if any(k for k in glob.glob(os.path.join(OUT, f"{RC.PRIMARY}__card__+0.[0-9]*.json")) if "+0.000" not in k):
        out["control"] = analyze(RC.PRIMARY, "control")
    json.dump(out, open(SUMMARY, "w"), indent=1, default=float)
    for k, v in out.items():
        print(k, {kk: v.get(kk) for kk in ("result_state", "mechanical_S5_output", "load_bearing_flagged_gt_0p2", "eps_hat", "I_card", "I_CPL", "F", "profile_set_95",
                                            "boundary_null_chernoff_p", "start_disagreement_flags")})


if __name__ == "__main__":
    main()
