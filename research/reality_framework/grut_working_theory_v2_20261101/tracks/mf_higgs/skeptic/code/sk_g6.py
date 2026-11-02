#!/usr/bin/env python3
"""mf_higgs SKEPTIC: the hypothetical G6 portal lam_chiH (chi^2 - v^2)|H|^2, lam_chiH < 0 (the benchmark has none).
m_H^2 = Lam H_i^2 (1 - chi^2/v^2), Lam = |lam_chiH| v^2/H_i^2, v = 0.0196 M_P = 3.67e5 H_i.  [D/C, model choices]
  (1) virial check of the period mean <1 - chi^2/v^2> on the exact anharmonic valley oscillation, against the track's
      (A/v)^2 = 2 rho_chi/(m_chi^2 v^2) = f_chi H^2/4;
  (2) Floquet exponents of h'' + [k^2 + Lam (1 - chi(t)^2/v^2)] h = 0 at fixed amplitude (no expansion), and the
      per-mode growth through the first band in an expanding background, mu * (band crossing time ~ q/2H);
  (3) Lam <-> lam_chiH, and the tree-stability and loop sizes.
Output: results/sk_g6.json
"""
import sys
sys.dont_write_bytecode = True
import math

import numpy as np
from scipy.integrate import solve_ivp

import sk_common as C

MCHI = 249.995
V_HI = 0.019595918 * C.MPL / 1.30002426e11


def chi_orbit(Av, n=4000):
    """exact oscillation of chi = v (1 + d) in V = lam/4 (chi^2 - v^2)^2, m_chi^2 = 2 lam v^2 (units m_chi = 1, v = 1):
    d'' = -(1/2)(1 + d)((1 + d)^2 - 1).  Start at the outer turning point d = +Av.  Returns t, d over one period,
    the period, and rho/(m^2 v^2/2) = energy normalised so that the harmonic value is Av^2."""
    f = lambda t, y: [y[1], -0.5 * (1 + y[0]) * ((1 + y[0]) ** 2 - 1)]
    ev = lambda t, y: y[1]
    ev.direction = 1
    sol = solve_ivp(f, (0, 40), [Av, 0.0], rtol=1e-11, atol=1e-13, events=ev, dense_output=True, max_step=0.01)
    Th = sol.t_events[0][0] if sol.t_events[0][0] > 1.0 else sol.t_events[0][1]
    # first upward zero of d' = inner turning point = half period; mirror (time reversal) for the full period
    th = np.linspace(0, Th, n)
    dh = sol.sol(th)[0]
    tt = np.concatenate([th, 2 * Th - th[::-1][1:]])
    d = np.concatenate([dh, dh[::-1][1:]])
    T = 2 * Th
    E = 0.25 * 0.5 * ((1 + Av) ** 2 - 1) ** 2      # V/(m^2 v^2) with lam v^2 = 1/2 -> V = (1/8)(chi^2-1)^2 ... in m^2 v^2 units
    rho_norm = 2 * E                                # rho/(m^2 v^2 / 2)  (harmonic limit: Av^2)
    return tt, d, T, rho_norm


def virial():
    rows = []
    for Av in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
        tt, d, T, rn = chi_orbit(Av)
        g = 1 - (1 + d) ** 2
        mean = float(np.trapezoid(g, tt) / T)
        a = math.sqrt(rn)
        rows.append(dict(A_over_v_turning=Av, rho_over_half_m2v2=rn, mean_1_minus_chi2_over_v2=mean,
                         track_formula_2rho_over_m2v2=rn, ratio=mean / rn, period_over_2pi_m=T / (2 * math.pi),
                         track_rebuilt_min=-2 * a + rn, track_rebuilt_max=2 * a + rn,
                         min_1_minus_chi2=float(g.min()), max_1_minus_chi2=float(g.max())))
    return rows


def floquet(Lam, rn, ks, H=1.0, harmonic=False):
    """Floquet exponents (per H_i time) of h'' + [k^2 + Lam(1 - chi^2/v^2)] h = 0 for the orbit with energy rn
    (= f_chi H^2/4 in the track's normalisation (A/v)^2)."""
    # turning amplitude for that energy: solve ((1+A)^2 - 1)^2/4 = rn
    A = math.sqrt(1 + 2 * math.sqrt(rn)) - 1
    if harmonic:      # the track's rebuilt form: -2 a cos(m t) + a^2, a = sqrt(rn), period 2 pi/m
        T = 2 * math.pi
        tt = np.linspace(0, T, 6000)
        g = -2 * math.sqrt(rn) * np.cos(tt) + rn
    else:
        tt, d, T, _ = chi_orbit(A, n=6000)
        g = 1 - (1 + d) ** 2
    Tphys = T / MCHI
    out = []
    for k in ks:
        f = lambda t, y: [y[1], -(k * k + Lam * np.interp((t % Tphys) * MCHI, tt, g)) * y[0]]
        M = np.zeros((2, 2))
        for j, y0 in enumerate(([1.0, 0.0], [0.0, 1.0])):
            s = solve_ivp(f, (0, Tphys), y0, rtol=1e-10, atol=1e-12, max_step=Tphys / 400)
            M[:, j] = s.y[:, -1]
        tr = np.trace(M) / 2
        mu = math.acosh(abs(tr)) / Tphys if abs(tr) > 1 else 0.0
        out.append(mu)
    return np.array(out), A


def main():
    out = dict(conversion=dict(v_over_Hi=V_HI, lam_chiH_per_Lam=1 / V_HI ** 2), virial=virial())
    print('virial', [(r['A_over_v_turning'], round(r['ratio'], 3)) for r in out['virial']])
    # amplitudes on the two histories: (A/v)^2 = f_chi H^2/4 (track's normalisation)
    cases = {'EXP t=0.25 (f_chi 0.56, H 0.91)': (0.56, 0.91), 'EXP t=1 (f_chi 0.15, H 0.40)': (0.15, 0.40),
             'Pauli t=0.3 (f_chi 0.83, H 0.85)': (0.83, 0.85), 'Pauli t=3 (f_chi 0.86, H 0.19)': (0.86, 0.19)}
    ks = np.concatenate([np.linspace(0, 40, 9), np.linspace(60, 160, 21)])
    fl = {}
    for Lam in (1000.0, 3000.0, 1e4, 3e4):
        for lab, (fc, H) in cases.items():
          rn = fc * H ** 2 / 4
          for harm in (False, True):
            mu, A = floquet(Lam, rn, ks, harmonic=harm)
            j = int(np.argmax(mu))
            q = 2 * Lam * 2 * math.sqrt(rn) / MCHI ** 2     # Mathieu q for the fundamental line (amplitude 2 Lam A/v)
            key = f'Lam={Lam:g} | {lab} | ' + ('harmonic (track form)' if harm else 'exact orbit')
            fl[key] = dict(A_over_v=A, mu_max=float(mu[j]), k_at_max=float(ks[j]), mu_max_over_H=float(mu[j] / H),
                                               q=q, per_mode_band_growth=float(mu[j] * q / (2 * H)) if ks[j] > 0 else None,
                                               tachyonic_k0=float(mu[0]))
            print(key, {k: round(v, 4) if isinstance(v, float) else v for k, v in fl[key].items()}, flush=True)
    out['floquet'] = fl
    out['notes'] = ('mu in H_i units at fixed amplitude, no expansion; per_mode_band_growth = mu_max x (q/2H), the time a mode '
                    'spends in a narrow band of relative width ~q/2 while k/a redshifts; tachyonic_k0 = mu at k = 0')
    C.jdump(out, C.RES / 'sk_g6.json')


if __name__ == '__main__':
    main()
