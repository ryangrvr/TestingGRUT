"""Where inside the rigorous interval [(1-s)^-1/2, (1-s)^-1] does R = E/(sigma A_area) fall for smooth, unimodal
speed distributions with energy-weighted s_E = <v^2>_E = 0.456^2?  Families (area measure, normal speed v):
  half-normal in gamma*v, uniform in v on [0, vmax], beta(a, b) in v, and a two-component 'rest + fast' mixture
  (the bound-saturating shape).  For each family the scale is solved so that s_E = target.  [E]-level illustration:
  the true 3D PRS distribution was not measured.

Usage: python3 -I realistic_R.py <out_json>
"""
import json
import math
import sys

import numpy as np
from scipy import optimize, stats

TARGET = 0.456 ** 2
V = np.linspace(0.0, 0.999999, 200001)


def moments(wA):
    g = 1 / np.sqrt(1 - V ** 2)
    wA = wA / np.trapezoid(wA, V)
    E = np.trapezoid(wA * g, V)
    sE = np.trapezoid(wA * g * V ** 2, V) / E
    sA = np.trapezoid(wA * V ** 2, V)
    return sE, sA, E      # E = <gamma>_A = R


def solve(family, lo, hi):
    f = lambda p: moments(family(p))[0] - TARGET
    p = optimize.brentq(f, lo, hi)
    sE, sA, R = moments(family(p))
    return {"param": p, "s_E": sE, "s_A": sA, "R": R, "R_minus_lower_over_width":
            (R - (1 - TARGET) ** -0.5) / ((1 - TARGET) ** -1 - (1 - TARGET) ** -0.5)}


def main(out):
    res = {"lower": (1 - TARGET) ** -0.5, "upper": 1 / (1 - TARGET)}
    u = V / np.sqrt(1 - V ** 2)                                   # gamma v
    jac = (1 - V ** 2) ** -1.5                                    # d(gamma v)/dv
    res["halfnormal_in_gamma_v"] = solve(lambda s: np.exp(-0.5 * (u / s) ** 2) * jac, 0.05, 3.0)
    res["maxwell_like_3D_in_gamma_v"] = solve(lambda s: u ** 2 * np.exp(-0.5 * (u / s) ** 2) * jac, 0.05, 2.0)
    res["uniform_in_v"] = solve(lambda vm: (V <= vm).astype(float) + 1e-300, 0.3, 0.999)
    for a, b in ((2, 2), (2, 4), (1, 3), (0.5, 0.5)):
        def fam(c, a=a, b=b):
            x = np.clip(V / c, 0, 1)
            return np.where(V < c, stats.beta.pdf(np.clip(x, 1e-9, 1 - 1e-9), a, b), 0.0) + 1e-300
        try:
            res[f"beta_{a}_{b}_scaled"] = solve(fam, 0.3, 0.999)
        except ValueError as e:
            res[f"beta_{a}_{b}_scaled"] = str(e)
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
