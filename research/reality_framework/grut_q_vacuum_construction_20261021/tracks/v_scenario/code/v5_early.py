"""v_scenario / v5: independent check of the scenario's Planck-era relaxation numbers (NOTES §3).

Units M = 1 (K0 = 1, chi = 1, q0 = 3), KV quartic (v1 check C): m^2 = 2/3.
My formulation (as in v2): H solved from the Friedmann constraint 3H^2(1+Y) + 3H Ydot = rho_m + rho_V(Y) at every
step; Ydd = [S + 4 rho_V - (1+Y) R_tot]/3 - 3 H Ydot with R_tot = R(Y) + 3 Gamma Ydot/(1+Y)^2 (dissipative
extension; v1 check B shows matter must then gain Q = (3/2) Gamma Ydot^2/(1+Y)^2).
(a) KV 0806.2805 Figs 1-2 data (scenario mapping): Y_i = -0.5, H_i = 0.5/sqrt(3), rho_M = 0.05 (w_M = 1/3), t_i =
    sqrt(3) (KV tau = 1); t<H> and t^2<rho_V> averaged over KV tau = 1800-3000 (t = sqrt(3) tau).
(b) displaced q at rest, no matter, Gamma = 4 m^3/(192 pi) (Gorbunov-Panin (T4+), N_s = 4, read in TeX by this
    verifier): f_inj = int Q a^3 dt / (rho_V,i a_i^3), T_rh from rho_r at t = 1/Gamma (g* = 106.75), E-decay rate.
Output: ../results/v5_early.json
"""
import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp

m2 = 2.0 / 3.0
MRED = 2.435e18


def rhoV(Y):
    return 0.5 * Y * Y * (1 + 4 * Y / 3 + Y * Y / 3)


def Rq(Y):
    return 2 * Y * (1 + 2 * Y + 2 * Y * Y / 3)


def Hc(Y, V, rho):
    A = 3 * (1 + Y)
    return (-3 * V + math.sqrt(9 * V * V + 4 * A * (rho + rhoV(Y)))) / (2 * A)


def rhs_factory(Gam, wM):
    def rhs(t, s):
        Y, V, lna, rM, rr, Qint = s
        rho = rM + rr
        H = Hc(Y, V, rho)
        S = (1 - 3 * wM) * rM
        R = Rq(Y) + 3 * Gam * V / (1 + Y) ** 2
        Q = 1.5 * Gam * V * V / (1 + Y) ** 2
        Vd = (S + 4 * rhoV(Y) - (1 + Y) * R) / 3 - 3 * H * V
        return [V, Vd, H, -3 * (1 + wM) * H * rM, -4 * H * rr + Q, Q * math.exp(3 * lna)]
    return rhs


out = {}
# ---------------------------------------------------------------- (a) KV relaxation
Yi, Hi, rM = -0.5, 0.5 / math.sqrt(3), 0.05
Vi = (rM + rhoV(Yi) - 3 * Hi * Hi * (1 + Yi)) / (3 * Hi)
ti = math.sqrt(3)
sol = solve_ivp(rhs_factory(0.0, 1 / 3), (ti, math.sqrt(3) * 3000 + 10), [Yi, Vi, 0.0, rM, 0.0, 0.0],
                method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True, max_step=0.2)
tt = np.linspace(math.sqrt(3) * 1800, math.sqrt(3) * 3000, 400001)
ss = sol.sol(tt)
Hs = np.array([Hc(a, b, c + d) for a, b, c, d in zip(ss[0], ss[1], ss[3], ss[4])])
out['KV_relaxation'] = dict(H_i=Hi, V_i=Vi, rhoV_i=rhoV(Yi),
                            t_H_mean=float(np.mean(tt * Hs)), t2_rhoV_mean=float(np.mean(tt ** 2 * rhoV(ss[0]))),
                            KV_prediction=2 / 3,
                            rhoV_over_rhoc_ratio_of_means=float(np.mean(tt ** 2 * rhoV(ss[0])) / (3 * np.mean(tt * Hs) ** 2)),
                            scenario=dict(t_H=0.6655, t2_rhoV=0.6620))

# ---------------------------------------------------------------- (b) decay
Gam = 4 * m2 ** 1.5 / (192 * math.pi)
dec = []
for yi in (-0.1, -0.03, 0.03, 0.1, 0.25):
    rV = rhoV(yi)
    Hi = Hc(yi, 0.0, 0.0)
    rhs = rhs_factory(Gam, 0.0)
    t_end = 14 / Gam
    sol = solve_ivp(rhs, (0.0, t_end), [yi, 0.0, 0.0, 0.0, 0.0, 0.0], method='DOP853', rtol=1e-10, atol=1e-16,
                    dense_output=True, max_step=0.3)
    f_inj = sol.y[5, -1] / rV
    s1 = sol.sol(1 / Gam)
    rr1 = s1[4]
    T_rh = (30 * rr1 / (math.pi ** 2 * 106.75)) ** 0.25 * MRED
    # Einstein-frame oscillation energy decay rate from late windows (Y small): E = (3/4)(V^2 + m^2 Y^2) a^3
    tw = np.linspace(4 / Gam, 10 / Gam, 13)
    Ea3 = []
    for tc in tw:
        ts = tc + (2 * math.pi / math.sqrt(m2)) * (np.arange(32) / 32 - 0.5)
        q = sol.sol(ts)
        Ea3.append(np.mean(0.75 * (q[1] ** 2 + m2 * q[0] ** 2) * np.exp(3 * q[2])))
    rate = -np.polyfit(tw, np.log(Ea3), 1)[0] / Gam
    sl = sol.sol(t_end)
    Hl = Hc(sl[0], sl[1], sl[4])
    f_com = sl[4] * math.exp(4 * sl[2]) / rV
    dec.append(dict(y_i=yi, rhoV_i=rV, H_i_over_m=Hi / math.sqrt(m2), f_inj=float(f_inj), f_com=float(f_com),
                    T_rh_GeV=float(T_rh), decay_rate_over_Gamma=float(rate), late_tH=float(t_end * Hl)))
out['decay'] = dict(Gamma_over_M=Gam, Gamma_over_m=Gam / math.sqrt(m2), runs=dec,
                    GP_T5plus_GeV=4.5e-2 * 106.75 ** -0.25 * math.sqrt(4 * (math.sqrt(m2) * MRED) ** 3 / MRED),
                    scenario=dict(f_inj={'-0.10': 1.062, '-0.03': 1.016, '+0.03': 0.986, '+0.10': 0.959, '+0.25': 0.917},
                                  T_rh={'-0.10': 5.05e16, '-0.03': 4.82e16, '+0.03': 4.83e16, '+0.10': 5.06e16,
                                        '+0.25': 5.12e16}))
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, '..', 'results', 'v5_early.json'), 'w') as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
