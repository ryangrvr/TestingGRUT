#!/usr/bin/env python3
"""Leftover budget, remainder laws and the transition-energy ledger (own code, track gp_leftovers).

1. Reproduce Codex's cold-pair budget (b5174e5 POINT_RESULTS symmetric_remainder).
2. Remainder laws f_extra(eps):
   (a) quoted review law 0.9(eps^-4/3 - 1);
   (b) NR particle count, eps = (N+ - N-)/(N+ + N-): f = (m_X Y_Q/B)(1/eps - 1);
   (c) adiabatic map from the quartic-era AD fraction eps_q = n/n_max(rho) to the NR count, using the
       adiabatic invariants a^3 n and a^3 J_r of the planar orbit (exact for a slowly varying central potential);
   (d) energy cap at the phi-transition onset: late cold matter per net charge <= (rho_X + V0)/n_Q.
3. Ledger content of the calibrated 4.4e-10 GeV (objects + binding debris), Codex's quench benchmark,
   its self-consistent version, the double-counting test, the homogeneous first-order test and the
   zero-pressure coexistence point where no vacuum energy is released.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from numpy.polynomial.legendre import leggauss
from common import *

OUT = Path(__file__).resolve().parents[1] / 'results'
XG, WG = leggauss(400)
TT = 0.5 * np.pi * (XG + 1)
WT = 0.5 * np.pi * WG


def orbit(Vf, rho, n):
    """Radial turning points, J_r = (1/2pi) oint p_R dR with p_R = 2 Rdot, and time-averaged R^2.
    Lagrangian density Rdot^2 + R^2 thdot^2 - V(R); n = 2 R^2 thdot."""
    f = lambda R: rho - n * n / (4 * R * R) - Vf(R)
    Rs = np.geomspace(1e-8, 10, 20001)
    fs = np.array([f(R) for R in Rs])
    i = np.argmax(fs)
    R2 = brentq(f, Rs[i], 10.0) if fs[-1] < 0 else None
    R1 = brentq(f, 1e-12, Rs[i]) if n > 0 else 0.0
    R = 0.5 * (R1 + R2) - 0.5 * (R2 - R1) * np.cos(TT)
    dR = 0.5 * (R2 - R1) * np.sin(TT)
    fv = np.maximum([f(r) for r in R], 0)
    Jr = (2 / np.pi) * np.sum(WT * np.sqrt(fv) * dR)
    with np.errstate(divide='ignore', invalid='ignore'):
        dt = np.where(fv > 0, dR / np.sqrt(fv), 0.0)
    R2avg = np.sum(WT * R * R * dt) / np.sum(WT * dt)
    return R1, R2, Jr, R2avg


def quartic_map(eps):
    """lambda=1, rho=1 (scale free). Returns j = J_r/n, eps_NR = 1/(1+2j), <R^2>/rho^(1/2)."""
    nmax = 2 * np.sqrt(2) * (1 / 3) ** 0.75
    n = eps * nmax
    _, _, Jr, R2 = orbit(lambda R: R**4, 1.0, n)
    j = Jr / n
    return j, 1 / (1 + 2 * j), R2


def harmonic_control():
    """V = R^2 (m = 1): exact rho = m (2 J_r + n)."""
    out = []
    for n in (0.2, 0.5, 0.9):
        _, _, Jr, _ = orbit(lambda R: R * R, 1.0, n)
        out.append(dict(n=n, Jr=Jr, exact=(1 - n) / 2, rel_err=abs(Jr - (1 - n) / 2) / ((1 - n) / 2)))
    return out


def main():
    res = {}
    # 1. pair budget
    rows = []
    for d in (0.01, 0.1, 1.0):
        yp = d * BDM / (2 * MX)
        rows.append(dict(allowed=d, Y_pair_max=yp, pairs_per_charge=yp / YQ,
                         eps_quoted_law=(1 + d / 0.9) ** (-0.75)))
    res['pair_budget'] = rows
    # ledger decomposition (portal-debris witness, Q=200000, lambda_phi=0.03) [Q] -> recomputed
    E_obj = E_PER_FREE_MASS * Q_REF * MX
    release = Q_REF * MX - E_obj
    gx = G**4 * V * V / (4 * np.pi * MS) * np.sqrt(1 - 2 * G**2 / LP)
    gd = HD**4 * V * V / (8 * np.pi * MS) * np.sqrt(1 - 2 * (G**2 + HD**2) / LP)
    ns = release / MS
    late = E_obj + 2 * ns * (gx * MX + gd * MD) / (gx + gd)
    res['ledger'] = dict(object_energy_GeV=E_obj, binding_release_GeV=release, scalars_per_object=ns,
                         ledger_pairs_per_net_charge=ns * gx / (gx + gd) / Q_REF,
                         late_mass_per_object_GeV=late, calibrated_per_charge_over_mX=B_PER_CHARGE / MX,
                         mX_YQ_over_B=MX * YQ / BDM, recomputed_B=YQ / Q_REF * late,
                         note='No phi-transition (V0) energy is in this ledger; objects assumed assembled from free cold X in the +v vacuum.')
    # 2. remainder laws
    c = MX * YQ / BDM
    ctrl = harmonic_control()
    epsg = [0.999, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01]
    lam = []
    for e in epsg:
        j, eN, R2 = quartic_map(e)
        lam.append(dict(eps=e, quoted=0.9 * (e ** (-4 / 3) - 1), NR_count=c * (1 / e - 1),
                        adiabatic_j=j, adiabatic_epsNR=eN, adiabatic=c * (1 / eN - 1), R2avg_over_sqrt_rho=R2))
    # slope near 1
    j1, _, _ = quartic_map(1 - 1e-4)
    slope_j = j1 / 1e-4
    # thresholds
    thr = []
    for d in (0.01, 0.1, 1.0):
        eNR = 1 / (1 + d / c)
        ea = brentq(lambda e: c * (1 / quartic_map(e)[1] - 1) - d, 0.05, 0.99999)
        thr.append(dict(allowed=d, eps_quoted=(1 + d / 0.9) ** (-0.75), eps_NR_count=eNR, eps_adiabatic_quartic=ea))
    res['remainder_laws'] = dict(harmonic_control=ctrl, table=lam, j_slope_near_eps1=slope_j,
                                 analytic_slope=(4 / 3) / (3 * np.sqrt(6) / (4 / 3) * (4 / 3)) if False else None,
                                 thresholds=thr, mX_YQ_over_B=c)
    # 3. transition energy
    V0 = LP * V**4 / 4
    Xc = np.sqrt(LP / 2) * V / G
    rho_c = 3 * LX * Xc**4
    n_c = 2 * np.sqrt(2 * LX) * Xc**3
    nmax = lambda rho: 2 * np.sqrt(2 * LX) * (rho / (3 * LX)) ** 0.75
    onset = []
    for e in (1.0, 0.9, 0.8794, 0.7, 0.5, 0.3, 0.1, 0.03):
        _, _, R2 = quartic_map(e) if e < 1 else (0, 1, 1 / np.sqrt(3))
        # criterion A (tracks/affleck_dine): rho_X = rho_c ; criterion B: <R^2>_t = Xc^2
        rB = (Xc**2 / (R2 / np.sqrt(LX))) ** 2 if e < 1 else rho_c  # R^2 = R2 * (rho/lambda)^(1/2)
        row = dict(eps=e)
        for tag, rho in (('A', rho_c), ('B', rB)):
            nq = e * nmax(rho)
            T = brentq(lambda T: YQ * entropy(T) - nq, 1, 1e4)
            row[tag] = dict(rho_X=rho, n_Q=nq, T_onset_GeV=T,
                            E_per_charge_over_mX=(rho + V0) / nq / MX,
                            cap_over_calibrated=(rho + V0) / nq / B_PER_CHARGE,
                            quench_over_calibrated=(2 * MX / MS) * V0 / nq / B_PER_CHARGE,
                            ledger_plus_quench=1 + (2 * MX / MS) * V0 / nq / B_PER_CHARGE)
        onset.append(row)
    codex = []
    for T in (40., 69., 86.):
        for gs in (86.25, 106.75):
            s = 2 * np.pi**2 * gs * T**3 / 45
            codex.append(dict(T=T, gs=gs, ratio=(2 * MX / MS) * V0 / s / BDM))
    eps_dc = brentq(lambda e: 1 + (2 * MX / MS) * V0 / (e * n_c) / B_PER_CHARGE - (rho_c + V0) / (e * n_c) / B_PER_CHARGE, 0.05, 1)
    # homogeneous first-order test at fixed charge density (circular orbit, eps=1)
    def Ehom(u, n):
        fX = lambda lx: n * n / (4 * np.exp(lx)) + G * G * u * np.exp(lx) + LX * np.exp(2 * lx)
        r = minimize_scalar(fX, bounds=(np.log(1e-8), np.log(1e3)), method='bounded', options=dict(xatol=1e-12))
        return LP / 4 * (u - V * V) ** 2 + r.fun, np.exp(r.x)
    us = np.linspace(0, V * V, 4001)
    def jump(n):
        E = np.array([Ehom(u, n)[0] for u in us])
        k = int(np.argmin(E))
        return E[0], E[k], us[k]
    E0, Emin, umin = jump(n_c)
    curvature = LP / 2 - G**4 / (6 * LX)
    # density where homogeneous global minimum leaves phi=0 (first-order point)
    def gap(lnn):
        E0_, Em_, _ = jump(np.exp(lnn))
        return Em_ - E0_ + 1e-30
    ns_ = np.geomspace(n_c, 3 * n_c, 25)
    gaps = [jump(n)[1] - jump(n)[0] for n in ns_]
    n_first = None
    for a, b, ga, gb in zip(ns_[:-1], ns_[1:], gaps[:-1], gaps[1:]):
        if ga < -1e-15 and gb >= -1e-15:
            n_first = float(brentq(lambda x: (lambda t: t[1] - t[0] + 1e-15)(jump(x)), a, b))
    R_coex = (V0 / LX) ** 0.25
    omega_coex = np.sqrt(2 * LX) * R_coex
    res['transition'] = dict(
        V0_GeV4=V0, Xc_GeV=Xc, rho_c=rho_c, n_c=n_c,
        onset=onset, codex_quench_reproduced=codex,
        eps_above_which_ledger_plus_quench_exceeds_onset_energy=eps_dc,
        homogeneous=dict(curvature_d2E_du2_at_spinodal=curvature, E0=E0, Emin=Emin, phi2_at_min=umin,
                         release_per_charge_over_mX=(E0 - Emin) / n_c / MX,
                         n_first_order_over_n_c=None if n_first is None else n_first / n_c,
                         R_first_order_over_Xc=None if n_first is None else (n_first / n_c) ** (1 / 3)),
        coexistence=dict(R_coex_GeV=R_coex, omega_over_mX=omega_coex / MX, E_per_charge_over_mX=omega_coex / MX,
                         a_ratio_coex_to_spinodal=R_coex / Xc,
                         note='P = lambda_X R^4 - V0 = 0: zero-pressure Q-matter coexists with the +v vacuum; quasi-static separation releases no V0.'),
        object_E_per_charge_over_mX=E_PER_FREE_MASS)
    OUT.mkdir(exist_ok=True)
    (OUT / 'budget.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res['pair_budget'], indent=1))
    print(json.dumps(res['ledger'], indent=1))
    print('harmonic control', ctrl)
    print('j slope', slope_j)
    for r in lam:
        print({k: round(v, 5) if isinstance(v, float) else v for k, v in r.items()})
    print(json.dumps(thr, indent=1))
    t = dict(res['transition']); t.pop('onset')
    print(json.dumps(t, indent=1))
    for r in onset:
        print(r['eps'], {k: {kk: round(vv, 4) for kk, vv in r[k].items()} for k in ('A', 'B')})


if __name__ == '__main__':
    main()
