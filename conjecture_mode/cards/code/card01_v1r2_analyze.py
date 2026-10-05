#!/usr/bin/env python3
"""CARD #1 v1R2 analysis (POST-DATA NUMERICAL VERIFICATION).  Lists every executed start (v1R + v1R2) per target,
applies the lower-minimum reproducibility criterion, and applies the UNCHANGED frozen decide() only if every
load-bearing target is reproduced.  usage: card01_v1r2_analyze.py V1R_OUT V1R2_OUT DRIVER_LOG SUMMARY_JSON
"""
import json, sys, os
import numpy as np
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import card01_run_config as RC                # noqa: E402
import card01_v1r2_lib as L                   # noqa: E402
from card01_d3_analyze import profile_set     # noqa: E402

LABEL = "POST-DATA NUMERICAL VERIFICATION RESULT - NOT PREREGISTERED v1"


def main():
    v1r, out, drv, summ = sys.argv[1:5]
    P = L.pool(v1r, out)
    log = open(drv).read() if os.path.exists(drv) else ""
    stageA_fail = "Stage A failed" in log
    core = ["card:0.0", "cpl", "card:-0.08", "free"]
    table = {}
    for t, rows in sorted(P.items()):
        ok, lo, why = L.reproduced(rows)
        table[t] = {"reproduced": ok, "lowest": lo, "detail": why,
                    "starts": sorted([{"file": r["_file"], "chi2": r["chi2_eff"], "converged": r.get("converged_cycles"),
                                       "eps": r["point"].get("card01_eps"), "w": r["point"].get("w"), "wa": r["point"].get("wa"),
                                       "H0": r["point"].get("H0")} for r in rows], key=lambda x: x["chi2"])}
    res = {"label": LABEL, "targets": table}
    if stageA_fail:
        res["numerical_status"] = "CARD-01-v1R2 - COMPUTATIONALLY UNRESOLVED AT AVAILABLE MINIMIZATION METHOD"
        res["core_status"] = {t: table[t]["detail"] + (" REPRODUCED" if table[t]["reproduced"] else " NOT REPRODUCED") for t in core}
        res["mechanical_S5_state"] = None
    else:
        grid = sorted(RC.EPS_GRID_PRIMARY)
        prof = {e: table[f"card:{round(e, 3)}"]["lowest"] for e in grid}
        eps = np.array(grid); chi = np.array([prof[e] for e in grid])
        i = int(np.argmin(chi)); fb = table["free"]["lowest"]
        free_pt = L.best(P["free"])["point"]["card01_eps"]
        cmin, ehat = (chi[i], float(eps[i])) if chi[i] <= fb else (fb, free_pt)
        lo, hi, lo_e, hi_e = profile_set(eps, chi, cmin, RC.PROFILE_SET_DCHI2)
        lb = {f"card:{round(float(eps[i]), 3)}", "card:0.0", "cpl", "free"}
        for val, hit in ((lo, lo_e), (hi, hi_e)):
            if not hit:
                lb |= {f"card:{round(float(e), 3)}" for e in (max([x for x in grid if x <= val], default=None),
                                                              min([x for x in grid if x >= val], default=None)) if e is not None}
        fails = {t: table[t]["detail"] for t in sorted(lb) if not table[t]["reproduced"]}
        if chi[i] < fb - L.TOL:
            fails["free"] = f"fixed-eps grid point {eps[i]} is {fb - chi[i]:.3f} below the free fit"
        chiL, chiC = table["card:0.0"]["lowest"], table["cpl"]["lowest"]
        I_card, I_cpl = chiL - cmin, chiL - chiC
        res.update({"profile": [{"eps": e, "chi2": prof[e], "reproduced": table[f"card:{round(e, 3)}"]["reproduced"]} for e in grid],
                    "load_bearing": sorted(lb), "load_bearing_failures": fails,
                    "chi2_Lambda": chiL, "chi2_card_best": cmin, "eps_hat": ehat, "free_eps": free_pt,
                    "chi2_w0wa": chiC, "w0wa_best": {k: L.best(P["cpl"])["point"][k] for k in ("w", "wa")},
                    "I_card": I_card, "I_CPL": I_cpl, "F": I_card / I_cpl if I_cpl > 0 else None,
                    "profile_set_95": [lo, hi], "profile_set_hits_grid_edge": [lo_e, hi_e],
                    "boundary_null_chernoff_p": float(0.5 * chi2dist.sf(max(I_card, 0), 1)) if I_card > 0 else 0.5})
        if fails:
            res["numerical_status"] = "CARD-01-v1R2 - PRIMARY NOT NUMERICALLY RESOLVED (load-bearing failure in Stage B)"
            res["mechanical_S5_state"] = None
        else:
            res["numerical_status"] = "CARD-01-v1R2 - PRIMARY NUMERICALLY RESOLVED"
            st = RC.decide(I_cpl, I_card, ehat, lo, hi, ehat <= RC.EPS_RANGE_BOUNDARY + 1e-3)
            res["mechanical_S5_state"] = f"{st} [{LABEL}]"
    json.dump(res, open(summ, "w"), indent=1, default=float)
    print({k: res.get(k) for k in ("numerical_status", "mechanical_S5_state", "load_bearing_failures", "core_status",
                                   "chi2_Lambda", "chi2_card_best", "eps_hat", "chi2_w0wa", "I_card", "I_CPL", "F",
                                   "profile_set_95")})


if __name__ == "__main__":
    main()
