"""T_c, order and phi_c, independent of the module. Writes results/transition_v.json."""
import json
import os
import numpy as np
from scipy.optimize import brentq
import vt_core as V

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
os.makedirs(OUT, exist_ok=True)
SCHEMES = ['oneloop', 'AE', 'Parwani', 'DR']


def curv0(T, p, mu, scheme):
    # curvature at phi=0; for oneloop the X cubic does not enter the curvature at 0
    return 2 * float(V.dVds(np.array([0.0]), T, p, mu, scheme)[0])


def T0_of(p, mu, scheme):
    return brentq(lambda T: curv0(T, p, mu, scheme), 1.50, 1.60, xtol=1e-14, rtol=1e-15)


def minima_scan(T, p, mu, scheme, phimax=0.5, n=20000):
    """return local minima phi>0 of V along phi (sign changes of dV/ds from - to +)."""
    ph = np.geomspace(1e-7, phimax, n)
    g = V.dVds(ph, T, p, mu, scheme)   # dV/dphi = 2 phi g
    idx = np.where((g[:-1] < 0) & (g[1:] >= 0))[0]
    mins = []
    for i in idx:
        r = brentq(lambda x: float(V.dVds(np.array([x]), T, p, mu, scheme)[0]), ph[i], ph[i + 1], xtol=1e-16)
        mins.append(r)
    # also maxima (+ to -)
    idx2 = np.where((g[:-1] > 0) & (g[1:] <= 0))[0]
    maxs = [brentq(lambda x: float(V.dVds(np.array([x]), T, p, mu, scheme)[0]), ph[i], ph[i + 1], xtol=1e-16) for i in idx2]
    return mins, maxs, g


def analyse(p, mu, scheme):
    T0 = T0_of(p, mu, scheme)
    res = dict(T0=T0)
    # search for coexisting minima above T0
    Ts = T0 + np.concatenate([[1e-9, 1e-8, 1e-7, 1e-6], np.linspace(1e-5, 2e-3, 40)])
    coexist = []
    for T in Ts:
        mins, maxs, g = minima_scan(T, p, mu, scheme)
        if mins:
            coexist.append((T, mins, maxs))
    res['n_T_with_coexisting_min_above_T0'] = len(coexist)
    if not coexist:
        res['order'] = 'second (no coexisting minimum found for T0 < T <= T0+2e-3 on phi in [1e-7,0.5])'
        res['Tc'] = T0
        res['phic_over_Tc'] = 0.0
        # check: just below T0 there is exactly one minimum
        mins, maxs, g = minima_scan(T0 - 1e-6, p, mu, scheme)
        res['minima_at_T0_minus_1e-6'] = mins
        # quartic at phi=0 at T0: d(dV/ds)/ds at s->0, V ~ (1/2) curv phi^2 + (lam_eff/4) phi^4 -> dV/ds = curv/2 + lam_eff s/2
        s = np.array([1e-6, 2e-6])
        gg = V.dVds(np.sqrt(s), T0, p, mu, scheme)
        res['lam_eff_over_lp_at_T0'] = float(2 * (gg[1] - gg[0]) / (s[1] - s[0]) / p['lp'])
        return res
    # first order: T1 = largest T with coexisting minimum; Tc from degenerate minima
    Tlo, Thi = T0, coexist[-1][0]
    # bracket T1 more precisely by bisection
    Tup = Thi + 1e-5
    while minima_scan(Tup, p, mu, scheme)[0]:
        Tup += 1e-5
    a, b = Thi, Tup
    for _ in range(60):
        m = 0.5 * (a + b)
        if minima_scan(m, p, mu, scheme, n=6000)[0]:
            a = m
        else:
            b = m
    T1 = a
    res['T1'] = T1

    def dV(T):
        mins = minima_scan(T, p, mu, scheme, n=6000)[0]
        pb = max(mins)
        return float(V.Veff(np.array([pb]), T, p, mu, scheme)[0] - V.Veff(np.array([0.0]), T, p, mu, scheme)[0]), pb
    Tc = brentq(lambda T: dV(T)[0], T0 + 1e-12, T1, xtol=1e-14)
    pc = dV(Tc)[1]
    res.update(order='first (perturbative)', Tc=Tc, phic_over_Tc=pc / Tc)
    # latent heat L = T dDeltaV/dT at Tc (envelope theorem), DeltaV = V_sym - V_brk
    h = 1e-7
    L = -Tc * (dV(Tc + h)[0] - dV(Tc - h)[0]) / (2 * h)
    rho_h = np.pi**2 / 30 * 4 * Tc**4
    res['latent_heat_over_rho_h'] = L / rho_h
    res['latent_heat_GeV4'] = L * 1e12
    return res


def main():
    out = {}
    for mode in ['fixed', 'running']:
        out[mode] = {}
        for lab, mu in [('mu0/2', 0.05), ('mu0', 0.1), ('2mu0', 0.2), ('v', 1.0)]:
            p = V.BENCH if mode == 'fixed' else V.run_couplings(V.BENCH, mu)
            out[mode][lab] = {sc: analyse(p, mu, sc) for sc in SCHEMES}
            print(mode, lab, {sc: (round(out[mode][lab][sc]['Tc'], 7), out[mode][lab][sc]['phic_over_Tc']) for sc in SCHEMES}, flush=True)
    # leading
    Pi = V.debye(V.BENCH)
    out['leading_Tc_over_v'] = float(np.sqrt(V.BENCH['lp'] / Pi['phi']))
    allres = [out['running'][l][s]['Tc'] for l in out['running'] for s in ['AE', 'Parwani', 'DR']]
    out['range_resummed_running'] = [min(allres), max(allres)]
    # AE phi daisy mass at T0
    p = V.BENCH
    T = out['fixed']['mu0']['AE']['Tc']
    out['AE_mu0_M_over_T'] = float(np.sqrt(p['mp'] + Pi['phi'] * T**2) / T)
    with open(os.path.join(OUT, 'transition_v.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)


if __name__ == '__main__':
    main()
