"""Frozen recipe (chosen on m=0 data only, before the m/g=0.25 holdout finished):
 1. For each x, fit omega1(N)^2 = M^2 + b/N^2 + c/N^4 on all N >= 14.   (primary)
    Variants for the error band: (b/N^2 only), (b/N^2 + c/N^3); N>=14 and N>=16.
 2. Fit M(x) - 2m = B + alpha/sqrt(x) + beta/x over x in {2,4,8,16,32}.  (primary)
    Variants: drop x=32; drop x=2; linear in 1/sqrt(x) on x in {8,16,32}.
 Report primary value and the min-max spread over all variant combinations.
"""
import json, sys, itertools
import numpy as np
def load(path, m):
    r = json.load(open(path)); out = {}
    for v in r.values():
        if abs(v["m"] - m) < 1e-12: out.setdefault(v["x"], []).append((v["N"], v["omega1"]))
    return {x: sorted(l) for x, l in out.items()}
FS = {"N2N4": [2, 4], "N2": [2], "N2N3": [2, 3]}
def finite_size(data, form, nmin):
    N = np.array([a for a, b in data], float); w = np.array([b for a, b in data]); s = N >= nmin
    A = np.c_[np.ones(s.sum()), *[1 / N[s] ** p for p in FS[form]]]
    c = np.linalg.lstsq(A, w[s] ** 2, rcond=None)[0]; return np.sqrt(c[0])
XF = {"quad_all": ([2, 4, 8, 16, 32], 2), "quad_no32": ([2, 4, 8, 16], 2),
      "quad_no2": ([4, 8, 16, 32], 2), "lin_fine": ([8, 16, 32], 1)}
def cont(Mx, m, xf):
    xs, deg = XF[xf]; xs = [x for x in xs if x in Mx]
    u = 1 / np.sqrt(np.array(xs)); y = np.array([Mx[x] for x in xs]) - 2 * m
    A = np.c_[[u ** p for p in range(deg + 1)]].T
    return np.linalg.lstsq(A, y, rcond=None)[0][0]
def analyse(path, m):
    d = load(path, m); vals = {}
    for form, nmin, xf in itertools.product(FS, [14, 16], XF):
        Mx = {x: finite_size(d[x], form, nmin) for x in d if x in [2, 4, 8, 16, 32] and len([1 for N, _ in d[x] if N >= nmin]) >= len(FS[form]) + 2}
        vals[(form, nmin, xf)] = cont(Mx, m, xf)
    prim = vals[("N2N4", 14, "quad_all")]
    Mx = {x: float(finite_size(d[x], "N2N4", 14)) for x in sorted(d) }
    return {"bulk_omega_minus_2m_by_x": {x: round(v - 2*m, 5) for x, v in Mx.items()},
            "primary": float(prim), "spread": [float(min(vals.values())), float(max(vals.values()))]}
if __name__ == "__main__":
    res = {"m0": analyse("spectrum_cache.json", 0.0)}
    try: res["m0.25"] = analyse("spectrum_cache_m025.json", 0.25)
    except Exception as e: res["m0.25"] = f"not available: {e}"
    targets = {"m0": 0.5641896, "m0.25": 0.51918}
    for k, t in targets.items():
        if isinstance(res[k], dict):
            res[k]["target"] = t; res[k]["rel_err_primary"] = (res[k]["primary"] - t) / t
            res[k]["pass_3pct"] = abs(res[k]["rel_err_primary"]) <= 0.03
    json.dump(res, open("extrapolation_results.json", "w"), indent=1); print(json.dumps(res, indent=1))
