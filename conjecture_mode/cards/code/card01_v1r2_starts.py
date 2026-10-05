#!/usr/bin/env python3
"""CARD #1 v1R2 Stage-A start points (frozen before execution; no likelihood evaluation).

Routes per target (route 1 = retained lowest v1R solution, not re-run; routes 2-4 = new runs, max 4 starts total):
  card:0.0 (LCDM): 1 v1R best | 2 highest-likelihood official DESI LCDM sample (overall best) |
                   3 best sample of a distinct official LCDM chain file | 4 best sample of a third chain file
  cpl            : 1 v1R best (~(-0.425,-1.759)) | 2,3 best samples of two OTHER official w0wa chain files
                   (excluding the file that seeded v1R arm B) | 4 best sample of the remaining file
  card:-0.08     : 1 v1R best | 2 v1R free-eps low solution with eps fixed to -0.08 |
                   3 v1R best at eps=-0.06 continued to -0.08 | 4 v1R best at eps=-0.10 continued to -0.08
  free           : 1 v1R free best | 2 v1R eps=-0.08 best with eps released (start -0.08) |
                   3 v1R eps=-0.10 best released (start -0.10) | 4 v1R eps=-0.12 best released (start -0.12)
Official-chain samples are STARTING LOCATIONS ONLY (owner ruling v1R2); never results.
usage: card01_v1r2_starts.py V1R_OUT CHAINS_ROOT OUT_JSON
"""
import glob, json, os, sys
import numpy as np

PRIMARY = "desi-bao-all_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing"
NUIS = ["A_planck", "amp_143", "amp_217", "amp_143x217", "n_143", "n_217", "n_143x217", "calTE", "calEE"]
COSMO = ["logA", "ns", "H0", "ombh2", "omch2", "tau"]


def v1r_best(v1r, key):
    rows = [json.load(open(f)) for f in glob.glob(os.path.join(v1r, f"V1R__{key}__*.json")) if not f.endswith("tmp.json")]
    b = min(rows, key=lambda r: r["chi2_eff"])
    return b, {k: v for k, v in b["point"].items()}


def chain_best(root, model, k, params):
    fn = os.path.join(root, "cobaya", model, PRIMARY, f"chain.{k}.txt")
    h = open(fn).readline().lstrip("#").split(); A = np.loadtxt(fn)
    i = int(np.argmin(A[:, h.index("minuslogpost")]))
    return {"source": f"official DESI {model}/{PRIMARY}/chain.{k}.txt row {i}",
            "point": {p: float(A[i, h.index(p)]) for p in params}}


def fixed(pt, eps):
    return {k: v for k, v in pt.items() if k != "card01_eps"}


def released(pt, eps):
    d = {k: v for k, v in pt.items() if k != "card01_eps"}; d["card01_eps"] = eps
    return d


def main():
    v1r, root, out = sys.argv[1], sys.argv[2], sys.argv[3]
    S = {}
    bL, pL = v1r_best(v1r, "card__+0.000")
    S["card:0.0"] = [{"route": 1, "source": f"retained v1R {bL['tag']}", "retained_chi2": bL["chi2_eff"], "point": pL},
                     {"route": 2, **chain_best(root, "base", 1, COSMO + NUIS)},
                     {"route": 3, **chain_best(root, "base", 2, COSMO + NUIS)},
                     {"route": 4, **chain_best(root, "base", 3, COSMO + NUIS)}]
    bC, pC = v1r_best(v1r, "cpl__+0.000")
    S["cpl"] = [{"route": 1, "source": f"retained v1R {bC['tag']}", "retained_chi2": bC["chi2_eff"], "point": pC},
                {"route": 2, **chain_best(root, "base_w_wa", 4, COSMO + NUIS + ["w", "wa"])},
                {"route": 3, **chain_best(root, "base_w_wa", 1, COSMO + NUIS + ["w", "wa"])},
                {"route": 4, **chain_best(root, "base_w_wa", 3, COSMO + NUIS + ["w", "wa"])}]
    b8, p8 = v1r_best(v1r, "card__-0.080"); bF, pF = v1r_best(v1r, "free__+0.000")
    _, p6 = v1r_best(v1r, "card__-0.060"); _, p10 = v1r_best(v1r, "card__-0.100"); _, p12 = v1r_best(v1r, "card__-0.120")
    S["card:-0.08"] = [{"route": 1, "source": f"retained v1R {b8['tag']}", "retained_chi2": b8["chi2_eff"], "point": fixed(p8, -0.08)},
                       {"route": 2, "source": f"v1R free-eps low solution {bF['tag']} projected to eps=-0.08", "point": fixed(pF, -0.08)},
                       {"route": 3, "source": "v1R best eps=-0.06 solution continued to eps=-0.08", "point": fixed(p6, -0.08)},
                       {"route": 4, "source": "v1R best eps=-0.10 solution continued to eps=-0.08", "point": fixed(p10, -0.08)}]
    S["free"] = [{"route": 1, "source": f"retained v1R {bF['tag']}", "retained_chi2": bF["chi2_eff"], "point": pF},
                 {"route": 2, "source": "v1R best eps=-0.08 solution, eps released (start -0.08)", "point": released(p8, -0.08)},
                 {"route": 3, "source": "v1R best eps=-0.10 solution, eps released (start -0.10)", "point": released(p10, -0.10)},
                 {"route": 4, "source": "v1R best eps=-0.12 solution, eps released (start -0.12)", "point": released(p12, -0.12)}]
    json.dump(S, open(out, "w"), indent=1)
    for k, v in S.items():
        print(k, [(r["route"], r["source"][:70]) for r in v])


if __name__ == "__main__":
    main()
