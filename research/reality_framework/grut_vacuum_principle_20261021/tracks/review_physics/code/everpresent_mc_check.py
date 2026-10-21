"""Independent spot check of the everpresent-Lambda Model 1 Monte Carlo (unimodular track sec. 4.4).

Own implementation: x = ln a grid (n steps per e-fold), RK4 for (tau, J0..J3) with S held fixed
inside a step, V = 8 pi J3, S += alpha xi sqrt(dV) after each step, stationary start
S_i ~ N(0, alpha^2 V_i) at a_start, termination when H^2 <= 0 at a grid point or RK stage.
Units: rho / (M^2 H_m^2), V H_m^4, H_m = Lambda-free H today (omega_m = 0.1430, omega_r = 4.18e-5).
"""
import json
import math
import sys

import numpy as np

OMH2, ORH2 = 0.1430, 4.18e-5
HM2 = OMH2 + ORH2
Om, Or = OMH2 / HM2, ORH2 / HM2


def run(alpha, nreal=8000, n_per_efold=80, a_start=1e-10, seed=1):
    rng = np.random.default_rng(seed)
    x0, x1 = math.log(a_start), 0.0
    nsteps = int(round((x1 - x0) * n_per_efold))
    dx = (x1 - x0) / nsteps
    A = math.sqrt(Or)
    a = a_start
    tau = a / A
    J = np.array([A**4 * tau**5 / 5, A**4 * tau**6 / 30, A**4 * tau**7 / 210, A**4 * tau**8 / 1680])
    Jall = np.tile(J[:, None], (1, nreal))
    tauall = np.full(nreal, tau)
    V = 8 * math.pi * Jall[3]
    S = rng.normal(0, alpha * np.sqrt(V))
    alive = np.ones(nreal, bool)

    def deriv(x, tau_, J_, S_):
        a_ = math.exp(x)
        rL = S_ / (8 * math.pi * J_[3])
        H2 = (3 * Om / a_**3 + 3 * Or / a_**4 + rL) / 3
        ok = H2 > 0
        H = np.sqrt(np.where(ok, H2, 1.0))
        dtau = 1 / (a_ * H)
        return dtau, np.vstack([a_**4 * dtau, J_[0] * dtau, J_[1] * dtau, J_[2] * dtau]), ok

    x = x0
    rec = {}
    for i in range(nsteps):
        k1t, k1, o1 = deriv(x, tauall, Jall, S)
        k2t, k2, o2 = deriv(x + dx / 2, tauall + dx / 2 * k1t, Jall + dx / 2 * k1, S)
        k3t, k3, o3 = deriv(x + dx / 2, tauall + dx / 2 * k2t, Jall + dx / 2 * k2, S)
        k4t, k4, o4 = deriv(x + dx, tauall + dx * k3t, Jall + dx * k3, S)
        alive &= o1 & o2 & o3 & o4
        Vold = 8 * math.pi * Jall[3]
        tauall = tauall + dx / 6 * (k1t + 2 * k2t + 2 * k3t + k4t)
        Jall = Jall + dx / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        x += dx
        Vnew = 8 * math.pi * Jall[3]
        S = S + alpha * rng.standard_normal(nreal) * np.sqrt(np.maximum(Vnew - Vold, 0))
        a_ = math.exp(x)
        rL = S / Vnew
        H2 = (3 * Om / a_**3 + 3 * Or / a_**4 + rL) / 3
        alive &= H2 > 0
    a_ = 1.0
    rL = S / (8 * math.pi * Jall[3])
    H2 = (3 * Om + 3 * Or + rL) / 3
    OmL = rL / (3 * np.where(H2 > 0, H2, np.nan))
    surv = alive
    out = {"alpha_red": alpha, "nreal": nreal, "terminated_frac": float(1 - surv.mean())}
    if surv.sum() > 0:
        o = OmL[surv]
        out["Omega0_16_50_84"] = [float(v) for v in np.percentile(o, [16, 50, 84])]
        out["P_pos_given_surv"] = float((o > 0).mean())
        out["P_today_window_all"] = float(((o > 0.64) & (o < 0.74)).sum() / nreal)
        out["H0_median_kms"] = float(100 * math.sqrt(HM2) * np.median(np.sqrt(H2[surv])))
    out["V_today_noLambda_check_H^4V"] = None
    return out


if __name__ == "__main__":
    res = {"runs": [run(a) for a in (0.02, 0.10, 0.15, 0.20, 0.30)]}
    res["linear_check_alpha0.02_std_vs_pred"] = None
    path = sys.argv[1] if len(sys.argv) > 1 else "everpresent_mc_check.json"
    with open(path, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
