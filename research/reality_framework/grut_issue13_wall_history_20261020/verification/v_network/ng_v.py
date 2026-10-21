#!/usr/bin/env python3
"""Verifier's own Nambu-Goto collapse ODE (thin wall, radiation era a ~ tau, from rest at tau_i).
Action S = -sigma int dtau a^Dh R^kappa sqrt(1 - R'^2)  ->  R'' = -(1 - R'^2)[kappa/R + Dh R'/tau].
kappa = 1 circle/cylinder, 2 sphere.  Derived independently (E-L equation re-derived by hand).
Returns collapse time tau_c (R = 0) and tau_10 (R = R0/10)."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma


def collapse(kappa, Dh, x, R0=1.0):
    ti = R0 / x

    def rhs(t, y):
        R, Rp = y
        return [Rp, -(1 - Rp * Rp) * (kappa / R + Dh * Rp / t)]

    def ev0(t, y):
        return y[0] - 1e-6 * R0
    ev0.terminal = True
    ev0.direction = -1

    def ev10(t, y):
        return y[0] - 0.1 * R0
    ev10.direction = -1
    sol = solve_ivp(rhs, [ti, ti + 50 * R0 + 50 * ti], [R0, 0.0], events=[ev0, ev10],
                    rtol=1e-11, atol=1e-13, method="DOP853")
    tc = sol.t_events[0][0]
    t10 = sol.t_events[1][0]
    return tc, t10, ti


if __name__ == "__main__":
    out = {"flat_exact": {"circle_pi_over_2": np.pi / 2, "sphere_int": gamma(0.25) ** 2 / (4 * np.sqrt(2 * np.pi)),
                          "circle_tau10_arccos": float(np.arccos(0.1))}}
    rows = []
    for name, kappa, Dh in [("circle_Dh2", 1, 2), ("cylinder_Dh3", 1, 3), ("sphere_Dh3", 2, 3)]:
        for x in (1e-3, 1e-2, 0.178, 1.0, 5.6, 31.6, 1000.0):
            tc, t10, ti = collapse(kappa, Dh, x)
            rows.append({"shape": name, "x": x, "tau_c_over_tau_i": tc / ti, "dtau_c_over_R0": tc - ti,
                         "tau_c_over_R0": tc, "dtau_10_over_R0": t10 - ti, "tau_10_over_R0": t10})
    out["ng"] = rows
    # c_inf at x = 1000 and cosmic-time ratio c_inf^2
    out["c_inf"] = {r["shape"]: r["tau_c_over_R0"] for r in rows if r["x"] == 1000.0}
    out["t_c_over_t_R_sphere"] = out["c_inf"]["sphere_Dh3"] ** 2
    json.dump(out, open("results/ng_v.json", "w"), indent=1)
    for r in rows:
        print(r)
    print(out["c_inf"], out["t_c_over_t_R_sphere"], out["flat_exact"])
