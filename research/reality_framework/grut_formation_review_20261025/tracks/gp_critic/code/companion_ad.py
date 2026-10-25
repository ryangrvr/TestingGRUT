#!/usr/bin/env python3
"""gp_critic: one-field Affleck-Dine release with the compulsory radial companion (own code).

Matched EFT of the 7737487 mediator (valid while kappa|X|^2 < M_Y^2):
  V = -cH H^2 |X|^2 [inflation only] + lam |X|^4 + (alpha6 X^6 + h.c.) + beta6 |X|^6,
  beta6 = -M_Y^2 kappa^2/D < 0,  r = 2|alpha6|/|beta6| = sqrt2 mu h / M_Y^2 < 1 (D > 0).
Units: X = X0 u, X0^2 = cH H^2/(2 lam), conformal time s = sqrt(2 lam) X0 tau, radiation a = h s, h = 1/sqrt(cH),
release at s = 1/h (a = 1). eta = |alpha6| X0^2/lam (tracks/affleck_dine normalisation).
  u'' = -|u|^2 u - 3 eta conj(u)^5 / a^2 + (3 eta / r) |u|^4 u / a^2       (arg alpha6 rotated into the phase)
Initial state: Hubble-held radial minimum at a supplied (DRT-type) phase phi0 = 6 theta0, zero physical velocity.
eps = Im(conj(u) u') / ((2|u'|^2 + |u|^4)/3)^(3/4) at late time; q = 2 Im(conj(u) u') (affleck_dine A_6 = q/(eta sin phi0)).
Usage: python3 companion_ad.py <out.json>
"""
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp


def initial_y(eta, r, phi0):
    k = eta * (2.0 / r - 2.0 * np.cos(phi0))          # b X0^2 / lam, b = |beta6| - 2|alpha6| cos phi0
    if abs(k) < 1e-12:
        return 1.0
    if k <= 0:
        # angular term dominates and is positive along phi0: minimum always exists
        k = abs(k)
        return (np.sqrt(1 + 6 * k) - 1) / (3 * k) if k > 0 else 1.0
    if 1 - 6 * k < 0:
        return None                                    # no Hubble-held minimum
    return (1 - np.sqrt(1 - 6 * k)) / (3 * k) if k > 1e-14 else 1.0


def run(eta, r, phi0, cH=1.0, a_end=1000.0, rtol=1e-10):
    h = 1 / np.sqrt(cH)
    y0 = initial_y(eta, r, phi0)
    if y0 is None:
        return None
    u0 = np.sqrt(y0) * np.exp(1j * phi0 / 6)
    v0 = h * u0
    inv_r = 0.0 if not np.isfinite(r) else 1.0 / r

    def rhs(s, w):
        u = w[0] + 1j * w[1]; v = w[2] + 1j * w[3]
        a = h * s
        acc = -abs(u)**2 * u - 3 * eta * np.conj(u)**5 / a**2 + 3 * eta * inv_r * abs(u)**4 * u / a**2
        return [v.real, v.imag, acc.real, acc.imag]

    def blow(s, w):
        return 50.0 - np.hypot(w[0], w[1])
    blow.terminal = True
    sol = solve_ivp(rhs, (1 / h, a_end / h), [u0.real, u0.imag, v0.real, v0.imag], method='DOP853',
                    rtol=rtol, atol=1e-13, events=blow)
    if sol.status == 1:
        return dict(runaway=True)
    u = sol.y[0, -1] + 1j * sol.y[1, -1]; v = sol.y[2, -1] + 1j * sol.y[3, -1]
    q = 2 * np.imag(np.conj(u) * v)
    E = 2 * abs(v)**2 + abs(u)**4
    return dict(runaway=False, y0=y0, q=q, eps=q / 2 / (E / 3)**0.75)


def main():
    out = {}
    # Control 1: no companion, small eta -> affleck_dine linear law (A_6 = 12.97, eps ~ 4.5 eta) [Q]
    c = run(1e-3, np.inf, np.pi / 2)
    out['control_no_companion_eta1e-3'] = dict(A6=c['q'] / 1e-3, eps_over_eta=c['eps'] / 1e-3)
    c2 = run(1e-3, np.inf, np.pi / 2, rtol=1e-12)
    out['control_tolerance_rel_change'] = abs(c2['q'] - c['q']) / abs(c['q'])
    c3 = run(1e-3, np.inf, np.pi / 2, a_end=1e4)
    out['control_a1e4_rel_change'] = abs(c3['q'] - c['q']) / abs(c['q'])
    m = run(1e-3, np.inf, -np.pi / 2)
    out['control_cp_mirror_sum'] = c['q'] + m['q']
    # Scan at the stated companion ratio
    r = np.sqrt(2) * 5e11 * 1e11 / 1e24
    phis = np.linspace(np.pi / 24, np.pi - np.pi / 24, 47)
    scan = []
    best = dict(eps=0.0)
    for phi0 in phis:
        edge = 1.0 / (12 * (1 / r - np.cos(phi0)))
        for frac in (0.25, 0.5, 0.75, 0.9, 0.99, 0.999):
            res = run(frac * edge, r, phi0)
            if res is None:
                continue
            row = dict(phi0=phi0, eta=frac * edge, frac_of_edge=frac, **res)
            scan.append(row)
            if not res['runaway'] and abs(res['eps']) > abs(best['eps']):
                best = row
        # just beyond the edge: no minimum
        assert run(1.001 * edge, r, phi0) is None
    out['r'] = r
    out['analytic_cap_linear'] = 4.5 * r / (12 * np.sqrt(1 - r**2))
    out['best'] = best
    out['n_runaway'] = int(sum(1 for s in scan if s['runaway']))
    out['n_runs'] = len(scan)
    # same scan with no companion (gp_phase-type family at fixed eta = best eta) for contrast
    out['no_companion_same_eta_phi0'] = run(best['eta'], np.inf, best['phi0'])
    # cH = 10 at the best point
    out['best_cH10'] = run(best['eta'], r, best['phi0'], cH=10.0)
    # r dependence of the cap (numerical) at phi0 = arccos(r), eta = 0.99 edge
    rs = {}
    for rr in (0.0707, 0.2, 0.45, 0.7):
        phi0 = np.arccos(rr)
        edge = 1.0 / (12 * (1 / rr - np.cos(phi0)))
        res = run(0.99 * edge, rr, phi0)
        rs[f'{rr:g}'] = dict(eta=0.99 * edge, **res)
    out['cap_vs_r'] = rs
    # yields at the existence edge, Y = q cH^(3/2) H^(3/2) / (2 lam D M^(3/2)), radiation release, instant reheating [D]
    lam, MRED, gstar = 1e-3, 2.435e18, 106.75
    D = 2 * np.pi**2 / 45 * gstar * (90 / (np.pi**2 * gstar))**0.75
    beta6 = 1e24 * 1e-6 / (1e48 - 2 * (5e11 * 1e11)**2)
    X_edge = np.sqrt(lam / (6 * beta6))
    ylds = {}
    for cH in (1.0, 10.0, 100.0):
        for frac in (0.9, 0.999):
            edge = 1.0 / (12 * (1 / r))
            res = run(frac * edge, r, np.pi / 2, cH=cH)
            H = np.sqrt(frac) * X_edge / np.sqrt(cH / (2 * lam))   # eta ∝ X0^2 ∝ H^2
            Y = res['q'] * cH**1.5 * H**1.5 / (2 * lam * D * MRED**1.5)
            ylds[f'cH={cH:g},frac={frac}'] = dict(H_GeV=H, eps=res['eps'], q=res['q'], Y=Y, Y_over_YQ=Y / 4.8294009e-9)
    out['yields_phi0_pi2'] = ylds
    out['scan'] = scan
    json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != 'scan'}, indent=1, default=float))


if __name__ == '__main__':
    main()
