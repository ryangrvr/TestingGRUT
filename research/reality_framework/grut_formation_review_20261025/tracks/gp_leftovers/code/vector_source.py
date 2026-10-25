#!/usr/bin/env python3
"""Flag from tracks/gp_matching: longitudinal A' produced during inflation (vector light, m_A/H = 3.5e-6 at the
inflationary point) and its decay A' -> X Xbar (open: m_A' = 6e4 GeV >> 2 m_X). [E] (own code)

GMR scaling (Graham, Mardon & Rajendran 2016) [R]: Omega_A/Omega_DM = sqrt(m_A/6e-6 eV) (H_I/1e14 GeV)^2,
assumed to hold with the vacuum mass and instant reheating (not checked for this Higgsing history).
Width Gamma(A'->X Xbar) = e^2 m_A beta^3/(48 pi). The decay plasma is then fed to the maximal-equilibrium
model of equilibrium.py (case H) with the net charge n_Q present, to obtain the frozen Xbar per net charge.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from common import *
import equilibrium as eq

OUT = Path(__file__).resolve().parents[1] / 'results'
RHO_DM_S = 4.37e-10  # GeV per entropy for Omega_c h^2 = 0.12 [R]


def plasma_relic(Tv0, rho_h):
    nQ0 = YQ * entropy(Tv0)
    T0, mu0, nc0, th0 = eq.solve_state(nQ0, MX, energy=rho_h + MX * nQ0, eps_c=MX)
    sigma = th0['s'] / nQ0
    Ym = th0['nmi'] / nQ0
    A = eq.A
    lna_prev = None; frozen = None
    for Tv in np.geomspace(Tv0, 1e-3, 260):
        nQ = YQ * entropy(Tv)
        lna = -np.log(gstar_s(Tv) ** (1 / 3) * Tv)
        T, mu, nc, th = eq.solve_state(nQ, MX, sigma=sigma)
        p_eq = th['nmi'] / nQ
        nmi = Ym * nQ; npl = nmi + nQ
        G = A['a'] * npl**2 * nmi / 4 + A['b'] * npl**3 / 6 + A['b'] * npl * nmi**2 / 6
        H = hubble(Tv)
        if lna_prev is not None:
            Ym = p_eq + (Ym - p_eq) * np.exp(-min(G / H * (lna - lna_prev), 700))
        lna_prev = lna
        if frozen is None and G / H < 1:
            frozen = (Tv, T / MX)
        if frozen and Tv < frozen[0] / 30:
            break
    return dict(T_h0_over_mX=T0 / MX, sigma=sigma, p_initial=th0['nmi'] / nQ0, freeze=frozen, p_final=Ym,
                f_extra=2 * Ym * MX * YQ / BDM)


def main():
    rows = []
    for HI in (4.7e13, 1e13, 1e12, 1e11):
        omega_ratio = np.sqrt(M_AP / 6e-15) * (HI / 1e14) ** 2
        rho_s = omega_ratio * RHO_DM_S
        width = E_GAUGE**2 * M_AP / (48 * np.pi)
        Tdec = brentq(lambda T: hubble(T) - width, 1e-3, 1e8)
        rho_h = rho_s * entropy(Tdec)
        pairs = rho_s / M_AP / YQ
        rel = plasma_relic(Tdec, rho_h)
        rows.append(dict(H_inf=HI, Omega_A_over_DM=omega_ratio, rho_A_over_s_GeV=rho_s, width_GeV=width,
                         T_vis_decay_GeV=Tdec, pairs_per_net_charge_at_decay=pairs,
                         rho_h_over_rho_rad=rho_h / (np.pi**2 / 30 * gstar(Tdec) * Tdec**4), **rel))
        print(rows[-1], flush=True)
    (OUT / 'vector_source.json').write_text(json.dumps(rows, indent=1, default=float))


if __name__ == '__main__':
    main()
