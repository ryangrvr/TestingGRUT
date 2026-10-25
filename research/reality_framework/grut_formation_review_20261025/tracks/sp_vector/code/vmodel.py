#!/usr/bin/env python3
"""Track sp_vector (Claude, own code): inflationary background and vector bookkeeping at the 373d37e leading point.

Reuses our earlier gp_matching model (tracks/gp_matching/code/model.py, canonical real components, H_i = 1 units)
for the inflationary minimum and the 6x6 Hessian. From these:
  - Higgs norm F_g^2 = 2 sum q_i^2 |z_i|^2 (q = 1, 3, 6) and m_A = e F_g during inflation and in vacuum (6 e h);
  - the phase-sector mass matrix, its Goldstone (zero) mode and the two physical phase modes;
  - the BD gauge angle alpha = H/(2 pi F_g) and the light physical phase spread;
  - the GMR scales a_* (H = m_vac), k_* = a_* m_vac and the e-folds between k_* exit and the end of inflation.
Usage: python3 vmodel.py ../results/background.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'gp_matching' / 'code'))
import model as gpm  # our earlier track's model (reuse permitted by the brief)

HI = 1.1211167304108865e13          # GeV, Codex leading point
CH = 1.0
E_GAUGE = 1e-7
QCH = np.array([1.0, 3.0, 6.0])
MRED = 2.435e18
YQ = 4.8294009076988175e-9
GSTAR = 106.75
PAIR_ALLOW = 0.0049296413205678555   # Codex / gp_leftovers 1% allowance per net charge
NX_CODEX = -238.64868                 # Codex comoving N_X (H_i^3 units), confirmed by sp_yledger
YX_OVER_YQ = -65.98227759048056       # Codex raw yield


def minimum(m):
    best = None
    rng = np.random.default_rng(20261025)
    for trial in range(48):
        x0 = rng.normal(size=6) * np.array([25, 25, 2, 2, 0.5, 0.5])
        r = minimize(m['V'], x0, jac=m['grad'], method='BFGS', options=dict(gtol=1e-13, maxiter=20000))
        if best is None or r.fun < best.fun - 1e-12:
            best = r
    z = best.x
    for _ in range(30):  # Newton polish on the gradient
        g = m['grad'](z); Hs = m['hess'](z)
        step = np.linalg.lstsq(Hs, g, rcond=None)[0]
        z = z - step
    return z, m['V'](z), np.linalg.norm(m['grad'](z))


def main(out):
    m = gpm.build(c_H=CH, Hinf=HI)
    z, V0, gnorm = minimum(m)
    X = (z[0] + 1j * z[1]) / np.sqrt(2); Y = (z[2] + 1j * z[3]) / np.sqrt(2); S = (z[4] + 1j * z[5]) / np.sqrt(2)
    mods = np.abs([X, Y, S])
    F2 = 2 * np.sum(QCH**2 * mods**2)
    Fg = np.sqrt(F2)
    m_inf = E_GAUGE * Fg                       # H_i units
    h = m['h']
    m_vac = E_GAUGE * 6 * h                     # = 6 e h, H_i units
    # Phase sector: rotate each canonical pair so the field is real positive, then take the angular block.
    Hs = m['hess'](z)
    phases = np.angle([X, Y, S])
    R = np.zeros((6, 6))
    for i, ph in enumerate(phases):
        c, s = np.cos(ph), np.sin(ph)
        R[2 * i:2 * i + 2, 2 * i:2 * i + 2] = [[c, -s], [s, c]]   # (radial, angular) -> (re, im)
    Hrot = R.T @ Hs @ R
    ang = [1, 3, 5]; rad = [0, 2, 4]
    Haa = Hrot[np.ix_(ang, ang)]
    Hrr = Hrot[np.ix_(rad, rad)]
    wa, va = np.linalg.eigh(Haa)
    wall, vall = np.linalg.eigh(Hs)
    # Goldstone direction in canonical angular coordinates: sqrt2|z_i| q_i
    gdir = np.sqrt(2) * mods * QCH; gdir /= np.linalg.norm(gdir)
    overlaps = [float(abs(va[:, j] @ gdir)) for j in range(3)]
    # Physical phase modes: canonical angular amplitude u_j (unit vector) -> phase theta_i = u_ji/(sqrt2|z_i|)
    modes = []
    for j in range(3):
        u = va[:, j]
        dtheta = u / (np.sqrt(2) * mods)       # rad per unit canonical amplitude (H_i units)
        psi2 = 2 * dtheta[1] - dtheta[2]        # source phase psi2 = 2 theta_Y - theta_Sigma
        psi1 = dtheta[1] - 3 * dtheta[0]        # kappa phase theta_Y - 3 theta_X
        modes.append(dict(m2_over_H2=float(wa[j]), goldstone_overlap=overlaps[j],
                          dtheta_per_canonical=[float(t) for t in dtheta],
                          dpsi2_per_canonical=float(psi2), dpsi1_per_canonical=float(psi1)))
    light = min((md for md in modes if md['goldstone_overlap'] < 0.5), key=lambda d: d['m2_over_H2'])
    nu = np.sqrt(9 / 4 - light['m2_over_H2'])
    # GMR scales (a = 1 at the end of inflation = release; RD with H = H_i/a^2)
    a_star = 1 / np.sqrt(m_vac)
    k_star = a_star * m_vac
    N_star = np.log(1 / k_star)
    alpha_BD = 1 / (2 * np.pi * Fg)
    dpsi_exit = abs(light['dpsi2_per_canonical']) / (2 * np.pi)      # per ln k at horizon exit
    dpsi_kstar_end = dpsi_exit * np.exp(-(1.5 - nu) * N_star)
    stoch_rms = abs(light['dpsi2_per_canonical']) * np.sqrt(3 / (8 * np.pi**2 * light['m2_over_H2']))
    # small-scale variance of psi2 from modes between k_* and the end of inflation (per-ln-k power summed)
    var_small = dpsi_exit**2 * (1 - np.exp(-2 * (1.5 - nu) * N_star)) / (2 * (1.5 - nu))
    res = dict(
        H_i_GeV=HI, c_H=CH, e=E_GAUGE,
        minimum_abs_over_H=[float(v) for v in mods], minimum_V=V0, grad_norm=gnorm,
        F_g_over_H=float(Fg), F_g_GeV=float(Fg * HI), F_vac_GeV=6e11, F_inf_over_F_vac=float(Fg / (6 * h)),
        m_A_inf_over_H=float(m_inf), m_A_inf_GeV=float(m_inf * HI), m_A_vac_GeV=float(m_vac * HI),
        m_A_vac_over_H=float(m_vac), mass_ratio_sq_vac_over_inf=float((m_vac / m_inf)**2),
        e_heavy_during_inflation=float(1 / Fg),
        goldstone_weights_wprime=[float(v) for v in QCH * 2 * mods**2 / F2],
        angular_eigen_m2_over_H2=[float(v) for v in wa], full_hessian_eigen=[float(v) for v in wall],
        radial_eigen_m2_over_H2=[float(v) for v in np.linalg.eigvalsh(Hrr)],
        phase_modes=modes, light_mode=light, light_nu=float(nu),
        alpha_BD_rad=float(alpha_BD),
        a_star=float(a_star), k_star_over_Hi=float(k_star), efolds_kstar_to_end=float(N_star),
        dpsi2_per_lnk_at_exit=float(dpsi_exit), dpsi2_per_lnk_kstar_end=float(dpsi_kstar_end),
        psi2_stochastic_rms=float(stoch_rms), psi2_smallscale_rms_above_kstar=float(np.sqrt(var_small)),
        GMR_equivalent_goldstone_angle=float(1 / (2 * np.pi * 6 * h)),
    )
    print(json.dumps(res, indent=1))
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else str(HERE.parent / 'results' / 'background.json'))
