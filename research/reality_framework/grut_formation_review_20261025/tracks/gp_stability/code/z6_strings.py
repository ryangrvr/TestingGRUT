#!/usr/bin/env python3
"""Tension of the winding-1 Sigma (U(1) -> Z6) string at the stated point (track gp_stability, own code).

L = -F^2/4 + |D Sigma|^2 - lam (|Sigma|^2 - eta^2)^2,  D = d - i g A,  g = 6e (Sigma has charge 6),
eta^2 = h^2/2 (|Sigma| = h/sqrt2).  m_v = g h = 6 e h, m_s = sqrt(2 lam) h, beta = m_s^2/m_v^2.
Nielsen-Olesen ansatz Sigma = eta f(x) e^{i theta}, A_theta = a(x)/(g r), x = m_v r:
    mu = 2 pi eta^2 eps(beta) = pi h^2 eps(beta),
    eps = int x dx [ f'^2 + (1-a)^2 f^2/x^2 + a'^2/x^2 + (beta/4)(f^2-1)^2 ]   (eps(1) = 1, BPS).
Solved in s = ln x so that beta up to 1e13 is resolved:
    f_ss = (1-a)^2 f + (beta/2) e^{2s} (f^2-1) f,   a_ss = 2 a_s - e^{2s} (1-a) f^2.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

E_G, H, LAM_S = 1e-7, 1e11, 1.0
M_PL = 1.220890e19          # GeV (G = 1/M_PL^2)


def solve_string(beta, xmin=None, xmax=40.0, n=4000, init=None, tols=(1e-4, 1e-6, 1e-8)):
    if xmin is None:
        xmin = 1e-5 / np.sqrt(max(beta, 1.0))
    s0, s1 = np.log(xmin), np.log(xmax)
    s = np.linspace(s0, s1, n)
    x = np.exp(s)
    k = np.sqrt(beta / 2)
    f = np.tanh(k * x * 0.8)
    a = x * x / (1 + x * x)
    y0 = np.vstack([f, np.gradient(f, s), a, np.gradient(a, s)])
    if init is not None:          # continuation: previous solution, core shifted by ln of the scale ratio
        y0 = init(s)

    def rhs(s, y):
        f, fs, a, as_ = y
        e2 = np.exp(2 * s)
        return np.vstack([fs, (1 - a) ** 2 * f + 0.5 * beta * e2 * (f * f - 1) * f,
                          as_, 2 * as_ - e2 * (1 - a) * f * f])

    def bc(ya, yb):
        return np.array([ya[1] - ya[0], ya[3] - 2 * ya[2], yb[0] - 1, yb[2] - 1])

    sol = None
    for tol in tols:
        new = solve_bvp(rhs, bc, s if sol is None else sol.x, y0 if sol is None else sol.y,
                        tol=tol, max_nodes=500000)
        if not new.success:
            assert sol is not None, new.message
            break
        sol = new
    # tension integral on a fine log grid (own composite Simpson in s)
    ss = np.linspace(s0, s1, 400001)
    f, fs, a, as_ = sol.sol(ss)
    e2 = np.exp(2 * ss)
    dens = fs ** 2 + (1 - a) ** 2 * f ** 2 + as_ ** 2 / e2 + 0.25 * beta * e2 * (f * f - 1) ** 2
    from scipy.integrate import simpson
    eps = simpson(dens, x=ss)
    # small-x remainder below xmin (f ~ c x, a ~ d x^2) is O(xmin^2 beta) and neglected
    parts = dict(grad=float(simpson(fs ** 2, x=ss)), cov=float(simpson((1 - a) ** 2 * f ** 2, x=ss)),
                 mag=float(simpson(as_ ** 2 / e2, x=ss)), pot=float(simpson(0.25 * beta * e2 * (f * f - 1) ** 2, x=ss)))
    return float(eps), parts, int(sol.x.size), sol


def shifted(sol, beta_old, beta_new):
    ds = 0.5 * np.log(beta_new / beta_old)
    def init(s):
        f = sol.sol(np.clip(s + ds, sol.x[0], sol.x[-1]))[0:2]
        a = sol.sol(np.clip(s, sol.x[0], sol.x[-1]))[2:4]
        return np.vstack([f, a])
    return init


def main():
    rows = []
    last = None
    for beta in [0.5, 1.0, 2.0, 10.0, 100.0, 1e4, 1e6, 1e8, 1e10]:
        eps, parts, nodes, last = solve_string(beta)
        rows.append(dict(beta=beta, eps=eps, eps_minus_half_ln_beta=eps - 0.5 * np.log(beta) if beta > 1 else None,
                         HHT_fit_1p04_beta0p195=1.04 * beta ** 0.195 if 1e-2 <= beta <= 1e2 else None,
                         parts=parts, nodes=nodes))
        print(rows[-1])
    beta_pt = 2 * LAM_S / 36 / E_G ** 2          # m_s^2/m_v^2 = 2 lam h^2 /(36 e^2 h^2), independent of h
    # Large-beta extrapolation: eps - ln(beta)/2 converges (rows for 1e8 and 1e10 differ by ~1e-6).
    c_hi = rows[-1]['eps'] - 0.5 * np.log(rows[-1]['beta'])
    c_lo = rows[-2]['eps'] - 0.5 * np.log(rows[-2]['beta'])
    eps_asym = 0.5 * np.log(beta_pt) + c_hi
    eps_unc = abs(c_hi - c_lo)
    # Direct attempt at the point (collocation becomes very stiff; recorded only if it converges)
    try:
        eps_direct, parts_pt, nodes_pt, _ = solve_string(beta_pt, init=shifted(last, rows[-1]['beta'], beta_pt),
                                                         tols=(1e-4, 1e-6))
    except AssertionError as exc:
        eps_direct, parts_pt = None, dict(direct_solve='did not converge: ' + str(exc))
    eps_pt = eps_asym
    eps_ctrl = eps_direct
    mu = np.pi * H * H * eps_pt
    Gmu = mu / M_PL ** 2
    Gmu_bps = np.pi * H * H / M_PL ** 2
    limits = dict(
        Planck2013_XXV_NG_95=1.5e-7, Planck2013_XXV_AH_95=3.2e-7,
        LVK_O3_excluded_above=dict(model_A_BOS=[9.6e-9, 1e-6], model_B_LRS=[4.0e-15, 6.3e-15],
                                   model_C1=[2.1e-15, 4.5e-15], model_C2=[4.2e-15, 7.0e-15]),
        note='[Q] Planck 2013 XXV (arXiv:1303.5085); LVK O3 stochastic search (arXiv:2101.12248), ranges over kink '
             'number with one cusp. PTA bounds near 1e-10 [R]. LISA forecast reach about 1e-17 [R]. Retrospective only.')
    out = dict(point=dict(e=E_G, h_GeV=H, lambda_Sigma=LAM_S, m_v_GeV=6 * E_G * H, m_s_GeV=np.sqrt(2 * LAM_S) * H,
                          beta=beta_pt, eps=eps_pt, eps_extrapolation_uncertainty=eps_unc, eps_direct_solve_if_converged=eps_ctrl, asymptotic_constant=c_hi, parts=parts_pt,
                          mu_GeV2=mu, Gmu=Gmu, Gmu_if_BPS_pi_h2=Gmu_bps,
                          flux_tube_radius_GeVinv=1 / (6 * E_G * H), scalar_core_GeVinv=1 / (np.sqrt(2 * LAM_S) * H),
                          Gmu_over_LVK_O3_modelB_low=Gmu / 4.0e-15, Gmu_over_LVK_O3_modelC1_low=Gmu / 2.1e-15, Gmu_over_Planck_NG=Gmu / 1.5e-7),
               scan=rows, limits=limits,
               scope='Classical Abelian-Higgs winding-1 string tension. Whether a string network forms (post-inflation '
                     'U(1) restoration and re-breaking at T ~ h) is not computed; GW bounds assume Nambu-Goto-like loops.')
    path = Path(__file__).resolve().parents[1] / 'results' / 'z6_strings.json'
    path.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out['point'], indent=1))


if __name__ == '__main__':
    main()
