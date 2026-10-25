#!/usr/bin/env python3
"""Maximal-equilibrium benchmark for the symmetric remainder (own code, track gp_leftovers).

Assumption (optimistic for erasure): right after the phi-transition onset the hidden X/Xbar population
thermalises completely (elastic + number changing) at fixed energy and charge, with a condensed phase at
chemical potential mu_c that absorbs charge beyond the thermal capacity:
  case H : homogeneous Bose condensate, mu_c = eps_c = m_X; Xbar destroyed only by 4->2 (amp42 coefficients,
           condensate quanta counted as reactants, Bose enhancement ignored);
  case L : lumps (Q=200000 proxy), mu_c = omega_Q = 0.8612 m_X, eps_c = E/Q = 0.9233 m_X; Xbar destroyed by
           4->2 among free quanta plus the bb211d1 conversion on lumps (n_lump = n_c/Q, quoted cold sigma v).
Then isentropic evolution (comoving S and N_Q conserved) with Bose-Einstein thermodynamics, and
Xbar relaxation dY/dlna = -(Gamma/H)(Y - Y_eq) integrated with exact exponential steps.
Initial energy per net charge: sudden onset (rho_c + V0)/n_Q = 0.8672 m_X/eps (criterion A of budget.py)
or adiabatic count m_X/eps_NR (budget.py map).  Output: frozen Xbar per net charge p_f(eps), f_extra.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from common import *

OUT = Path(__file__).resolve().parents[1] / 'results'
A = json.loads((OUT / 'amplitudes.json').read_text())['threshold_4to2']['A_GeV_minus8']
BUD = json.loads((OUT / 'budget.json').read_text())
SV_CONV = 1.13039e-16 / GEV2_CM3S
m = MX


def species(T, mu):
    """n, rho, P for one complex-scalar species with chemical potential mu (|mu| <= m), BE statistics."""
    def integ(k, which):
        def f(y):  # y = p/T ; E - mu computed without cancellation
            p = y * T; E = np.sqrt(p * p + m * m)
            x = (p * p / (E + m) + (m - mu)) / T
            w = {0: 1.0, 1: E, 2: p * p / (3 * E)}[which]
            if x >= 700:
                return 0.0
            if x < 1e-10:
                return (y * y / x if x > 0 else (E + m) / T) * w
            return y * y / np.expm1(x) * w
        # integrate over y with a split near the occupancy scale
        val = 0.0
        ys = np.sqrt(1 + 2 * m / T)
        for a_, b_ in ((0, 0.1 * ys), (0.1 * ys, ys), (ys, 5 * ys), (5 * ys, 20 * ys), (20 * ys, 80 * ys)):
            val += quad(f, a_, b_, limit=200, epsabs=0, epsrel=1e-9)[0]
        return val * T**3 / (2 * np.pi**2)
    return integ(0, 0), integ(0, 1), integ(0, 2)


def thermo(T, mu):
    npl, rpl, ppl = species(T, mu)
    nmi, rmi, pmi = species(T, -mu)
    rho, P = rpl + rmi, ppl + pmi
    nq = npl - nmi
    s = (rho + P - mu * nq) / T
    return dict(npl=npl, nmi=nmi, nq=nq, rho=rho, P=P, s=s)


def solve_state(nQ, mu_c, sigma=None, energy=None, eps_c=None):
    """Return (T, mu, n_cond, thermo) for given charge and either entropy or energy density.
    T_b: temperature at which the thermal gas at mu = mu_c holds exactly nQ (onset of the condensed phase)."""
    gb = lambda lT: thermo(np.exp(lT), mu_c)['nq'] - nQ
    T_b = 1e-14 * m if gb(np.log(1e-14 * m)) >= 0 else np.exp(brentq(gb, np.log(1e-14 * m), np.log(1e4 * m), xtol=1e-12))
    def target2(T):
        th = thermo(T, mu_c)
        if sigma is not None:
            return th['s'] - sigma * nQ, th
        return eps_c * (nQ - th['nq']) + th['rho'] - energy, th
    tb, thb = target2(T_b)
    if tb >= 0:
        lo = np.log(1e-15 * m)
        if target2(np.exp(lo))[0] >= 0:   # extremely cold, fully condensed: below the resolved range
            T = np.exp(lo); th = thermo(T, mu_c)
            return T, mu_c, max(nQ - th['nq'], 0.0), th
        lT = brentq(lambda lT: target2(np.exp(lT))[0], lo, np.log(T_b), xtol=1e-12)
        T = np.exp(lT); th = thermo(T, mu_c)
        return T, mu_c, max(nQ - th['nq'], 0.0), th
    def mu_of_T(T):
        g = lambda mu: thermo(T, mu)['nq'] - nQ
        if g(mu_c) <= 0:
            return mu_c
        return brentq(g, 0.0, mu_c, xtol=1e-14 * m)
    def f1(lT):
        T = np.exp(lT); th = thermo(T, mu_of_T(T))
        return th['s'] - sigma * nQ if sigma is not None else th['rho'] - energy
    lT = brentq(f1, np.log(T_b), np.log(1e4 * m), xtol=1e-11)
    T = np.exp(lT); mu = mu_of_T(T)
    return T, mu, 0.0, thermo(T, mu)


def run(eps, case, energy_rule):
    T_on = BUD['transition']['onset']
    row = next(r for r in T_on if abs(r['eps'] - eps) < 1e-9) if any(abs(r['eps'] - eps) < 1e-9 for r in T_on) else None
    rho_c, V0, n_c = BUD['transition']['rho_c'], BUD['transition']['V0_GeV4'], BUD['transition']['n_c']
    nQ0 = eps * n_c
    Tv0 = brentq(lambda T: YQ * entropy(T) - nQ0, 1, 1e3)
    if energy_rule == 'sudden':
        E0 = rho_c + V0
    else:
        from budget import quartic_map
        eN = quartic_map(eps)[1] if eps < 1 else 1.0
        E0 = m * nQ0 / eN
    mu_c, eps_c = (m, m) if case == 'H' else (OMEGA_REF, E_PER_FREE_MASS * m)
    if E0 < eps_c * nQ0:
        return dict(eps=eps, case=case, energy_rule=energy_rule, status='energy per charge below eps_c: no free-gas state',
                    E_per_charge_over_mX=E0 / nQ0 / m)
    T0, mu0, nc0, th0 = solve_state(nQ0, mu_c, energy=E0, eps_c=eps_c)
    sigma = th0['s'] / nQ0
    Ym = th0['nmi'] / (YQ * entropy(Tv0)) * YQ / nQ0 * entropy(Tv0) * 0 + th0['nmi'] / nQ0  # Xbar per net charge
    hist = []
    Tvs = np.geomspace(Tv0, 1e-6, 240)
    lna_prev = None
    frozen_at = None
    for Tv in Tvs:
        nQ = YQ * entropy(Tv)
        lna = -np.log(gstar_s(Tv) ** (1 / 3) * Tv)
        T, mu, nc, th = solve_state(nQ, mu_c, sigma=sigma)
        p_eq = th['nmi'] / nQ
        nmi = Ym * nQ
        npl = nmi + nQ - nc if case == 'L' else nmi + nQ
        # relaxation rate of n- (per Xbar): linearised (n+ n- - n_eq^2)[A_a n+ n-/4 + A_b (n+^2 + n-^2)/6]
        G = A['a'] * npl**2 * nmi / 4 + A['b'] * npl**3 / 6 + A['b'] * npl * nmi**2 / 6
        if case == 'L':
            G += nc / Q_REF * SV_CONV
        H = hubble(Tv)
        if lna_prev is not None:
            x = G / H * (lna - lna_prev)
            Ym = p_eq + (Ym - p_eq) * np.exp(-min(x, 700))
        lna_prev = lna
        if G / H < 1 and (not hist or hist[-1]['G_over_H'] >= 1):
            frozen_at = Tv          # last crossing below 1 (lumps may re-form and restart conversion)
        hist.append(dict(T_vis=Tv, T_h_over_m=T / m, mu_over_m=mu / m, cond_frac=nc / nQ, p_eq=p_eq, p=Ym, G_over_H=G / H))
        if case == 'H' and frozen_at is not None and Tv < frozen_at / 30:
            break
    c = MX * YQ / BDM
    return dict(eps=eps, case=case, energy_rule=energy_rule, T_vis_onset=Tv0, E_per_charge_over_mX=E0 / nQ0 / m,
                T_h0_over_m=T0 / m, mu0_over_m=mu0 / m, cond_frac0=nc0 / nQ0, sigma=sigma,
                p_initial_equilibrium=th0['nmi'] / nQ0, T_vis_freeze=frozen_at, p_final=Ym,
                f_extra_pairs=2 * Ym * c,
                free_charge_fraction_at_freeze=(1 - next(h['cond_frac'] for h in hist if h['T_vis'] <= frozen_at)) if frozen_at else None,
                T_h_over_m_at_freeze=next(h['T_h_over_m'] for h in hist if h['T_vis'] <= frozen_at) if frozen_at else None,
                min_cond_frac=min(h['cond_frac'] for h in hist), final_cond_frac=hist[-1]['cond_frac'],
                history=hist[::8] + [hist[-1]])


def _job(args):
    import warnings
    warnings.filterwarnings('ignore')
    try:
        return run(*args)
    except Exception as exc:  # report, do not hide
        return dict(eps=args[0], case=args[1], energy_rule=args[2], status='solver failure: %r' % exc)


def main():
    from multiprocessing import Pool
    import sys
    grid = (0.95, 0.9, 0.7, 0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0.03, 0.01, 0.003, 0.001)
    rules = ('sudden', 'adiabatic') if 'full' in sys.argv else ('sudden',)
    jobs = [(eps, case, rule) for rule in rules for case in ('H', 'L') for eps in grid]
    with Pool(4) as pool:
        out = pool.map(_job, jobs)
    for r in out:
        print({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items() if k != 'history'}, flush=True)
    (OUT / 'equilibrium.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
