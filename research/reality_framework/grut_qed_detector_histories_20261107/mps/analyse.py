"""Implements PREREGISTRATION_AMENDMENT_01.json recipe verbatim. Written before production results were read."""
import json, itertools
import numpy as np

FS = {"L2L3_all(primary)": ([2, 3], 12), "L2_all": ([2], 12), "L2L4_all": ([2, 4], 12), "L2_Lge16": ([2], 16)}
CX = {"quad_all(primary)": ([9, 16, 36, 64], 2), "lin_16_36_64": ([16, 36, 64], 1), "quad_16_36_64": ([16, 36, 64], 2)}


def fs_fit(rows, powers, Lmin):
    r = [v for v in rows if v["L"] >= Lmin]
    L = np.array([v["L"] for v in r], float)
    w2 = np.array([v["omega1"] for v in r]) ** 2
    A = np.c_[np.ones(len(L)), *[L ** -p for p in powers]]
    return float(np.sqrt(np.linalg.lstsq(A, w2, rcond=None)[0][0]))


def cx_fit(Mx, m, xs, deg):
    u = 1 / np.sqrt(np.array(xs, float))
    y = np.array([Mx[x] for x in xs]) - 2 * m
    return float(np.polyfit(u, y, deg)[-1])


def analyse(m, target):
    d = json.load(open(f"production_m{m}.json"))
    by_x = {}
    for k, v in d.items():
        if k.startswith("x"):
            by_x.setdefault(v["x"], []).append(v)
    vals = {}
    for (fn, (pw, Lmin)), (cn, (xs, deg)) in itertools.product(FS.items(), CX.items()):
        Mx = {x: fs_fit(by_x[x], pw, Lmin) for x in by_x}
        vals[f"{fn} | {cn}"] = cx_fit(Mx, m, xs, deg)
    prim = vals["L2L3_all(primary) | quad_all(primary)"]
    Mx = {x: fs_fit(by_x[x], [2, 3], 12) - 2 * m for x in sorted(by_x)}
    lo, hi = min(vals.values()), max(vals.values())
    return {"m_over_g": m, "target": target, "primary": prim, "rel_err": (prim - target) / target,
            "verdict": "PASS" if abs(prim - target) / target <= 0.03 else "MISS",
            "spread": [lo, hi], "target_in_spread": lo <= target <= hi,
            "bulk_binding_by_x_primary_fs": Mx, "all_variants": vals,
            "chi160_check": d.get("chi160_check_x64_L24")}


if __name__ == "__main__":
    out = {"P-A": analyse(0.0, 0.5641896), "P-B": analyse(0.25, 0.51918)}
    json.dump(out, open("DMRG_RESULT.json", "w"), indent=1)
    for k, v in out.items():
        print(k, v["verdict"], round(v["primary"], 5), f"{100 * v['rel_err']:+.2f}%",
              [round(s, 4) for s in v["spread"]],
              "in spread" if v["target_in_spread"] else "outside spread",
              {x: round(b, 4) for x, b in v["bulk_binding_by_x_primary_fs"].items()})
