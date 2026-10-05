#!/usr/bin/env python3
"""CARD #1 v1R analysis (POST-DATA VERIFICATION-REPAIR).  Applies the UNCHANGED frozen S5 decide() only if every
load-bearing target passes the v1R convergence criteria.  usage: card01_v1r_analyze.py V1R_OUT SUMMARY_JSON
Also importable: load_bearing(), summarize().
"""
import glob, json, os, sys
import numpy as np
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_run_config as RC   # noqa: E402
import card01_v1r_config as C    # noqa: E402
from card01_d3_analyze import profile_set  # noqa: E402  (registered descriptive-crossing interpolation, unchanged)


def load(outdir):
    T = {}
    for f in glob.glob(os.path.join(outdir, "V1R__*.json")):
        r = json.load(open(f))
        key = "cpl" if r["model"] == "cpl" else ("free" if r["model"] == "free" else f"card:{r['eps_target']}")
        T.setdefault(key, {})[r["arm"]] = r
    return T


def best(t):
    return min((r for r in t.values()), key=lambda r: r["chi2_eff"])


def profile(T):
    pts = sorted(float(k.split(":")[1]) for k in T if k.startswith("card:"))
    return np.array(pts), np.array([best(T[f"card:{e}"])["chi2_eff"] for e in pts])


def load_bearing(T):
    eps, chi = profile(T)
    i = int(np.argmin(chi)); cmin = chi[i]; ehat = eps[i]
    if "free" in T and best(T["free"])["chi2_eff"] < cmin:
        cmin = best(T["free"])["chi2_eff"]; ehat = best(T["free"])["point"]["card01_eps"]
    lb = {"card:0.0", "cpl", "free", f"card:{eps[i]}"}
    lo, hi, lo_edge, hi_edge = profile_set(eps, chi, cmin, RC.PROFILE_SET_DCHI2)
    for val, hit in ((lo, lo_edge), (hi, hi_edge)):
        if hit:
            continue
        below = max([e for e in eps if e <= val], default=None); above = min([e for e in eps if e >= val], default=None)
        lb |= {f"card:{e}" for e in (below, above) if e is not None}
    return lb, (eps, chi, cmin, ehat, lo, hi, lo_edge, hi_edge)


def check(t):
    a, b = t.get("A"), t.get("B")
    if not a or not b:
        return False, "missing arm " + ("A" if not a else "B")
    if not (a["converged_cycles"] and b["converged_cycles"]):
        return False, "refinement cycles not converged"
    d = abs(a["chi2_eff"] - b["chi2_eff"])
    return d <= C.ACCEPT_AB, f"|A-B| = {d:.3f}"


def summarize(outdir):
    T = load(outdir)
    lb, (eps, chi, cmin, ehat, lo, hi, lo_edge, hi_edge) = load_bearing(T)
    table = []
    for k in sorted(T, key=lambda k: (k[:4] != "card", float(k.split(":")[1]) if k.startswith("card:") else 0)):
        t = T[k]; ok, why = check(t)
        table.append({"target": k, "load_bearing": k in lb,
                      "A_cycles": t.get("A", {}).get("cycle_chi2"), "A_converged": t.get("A", {}).get("converged_cycles"),
                      "B_cycles": t.get("B", {}).get("cycle_chi2"), "B_converged": t.get("B", {}).get("converged_cycles"),
                      "best_chi2": best(t)["chi2_eff"], "check_pass": ok, "check": why,
                      "best_point": {p: best(t)["point"].get(p) for p in ("H0", "omch2", "card01_eps", "w", "wa") if p in best(t)["point"]}})
    lb_fail = {r["target"]: r["check"] for r in table if r["load_bearing"] and not r["check_pass"]}
    chi_L = best(T["card:0.0"])["chi2_eff"]; chi_cpl = best(T["cpl"])["chi2_eff"]
    I_card = chi_L - cmin; I_cpl = chi_L - chi_cpl
    res = {"label": "POST-DATA VERIFICATION-REPAIR (v1R) - NOT THE PREREGISTERED v1 RESULT",
           "table": table, "load_bearing": sorted(lb), "load_bearing_failures": lb_fail,
           "chi2_Lambda": chi_L, "chi2_card_best": cmin, "eps_hat": float(ehat), "chi2_w0wa": chi_cpl,
           "w0wa_best": best(T["cpl"])["point"], "I_card": I_card, "I_CPL": I_cpl,
           "F": I_card / I_cpl if I_cpl > 0 else None, "profile_set_95": [lo, hi],
           "profile_set_hits_grid_edge": [lo_edge, hi_edge],
           "boundary_null_chernoff_p": float(0.5 * chi2dist.sf(max(I_card, 0), 1)) if I_card > 0 else 0.5}
    if lb_fail:
        res["numerical_status"] = "CARD-01-v1R - STILL NUMERICALLY UNRESOLVED"
        res["mechanical_S5_state"] = None
    else:
        res["numerical_status"] = "CARD-01-v1R - NUMERICALLY RESOLVED"
        st = RC.decide(I_cpl, I_card, ehat, lo, hi, ehat <= RC.EPS_RANGE_BOUNDARY + 1e-3)
        res["mechanical_S5_state"] = st + " [POST-DATA VERIFICATION-REPAIR RESULT - NOT THE PREREGISTERED v1 RESULT]"
    return res


if __name__ == "__main__":
    r = summarize(sys.argv[1])
    json.dump(r, open(sys.argv[2], "w"), indent=1, default=float)
    print({k: r[k] for k in ("numerical_status", "mechanical_S5_state", "load_bearing_failures", "chi2_Lambda",
                             "chi2_card_best", "eps_hat", "chi2_w0wa", "I_card", "I_CPL", "F", "profile_set_95")})
