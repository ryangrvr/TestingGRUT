#!/usr/bin/env python3
"""Two-quantum (anti-X pair) emission estimate, absent from the classical forced response (track xcharge).

Expanding lam X^d around the condensate gives  C(d,2) lam X0^(d-2) dX^2 + h.c., a parametric drive at
lab frequency (d-2) omega that creates anti-X pairs (energy (d-2) omega > 2 m_X for d = 5, 6 at
omega = 0.8612 m_X). It is O(|lam|^2) like the single-quantum rate but carries one power less of the
condensate field per amplitude, and it needs no momentum from the profile edge (pairs can be
back to back), so it is not suppressed by the steep surface form factor.

Free-Born estimate (vacuum dispersion for the emitted pair, finite source, no medium distortion):
  Gamma_events = 2 C^2 |lam|^2 v^(2(d-2)) 2^-(d-2) m_X^-3 J
  J = (1/(16 pi^3)) int dK1 K1^2 K2/eps1 int_{-1}^{1} dc |Gt(|K1+K2|)|^2,  eps1 + eps2 = (d-2) Omega
  Gt(Q) = 4 pi int rho^2 F^(d-2) j0(Q rho) drho
Local-density variants: vacuum threshold, and local effective mass m_eff^2 = S^2 + 2 A F^2.
Core charge loss per pair event: d-2. Total X violation per event: d.
"""
import numpy as np
from math import comb
from scipy.integrate import simpson

from profile import G, VPHI, MX, A as AQ
from forced import HBAR, MRED, T_REF


def Gt_table(bg, n, Qmax, nq=801):
    x = np.linspace(0, min(bg.R, 70.0), 14001)
    F, _ = bg.FS(x)
    Qs = np.linspace(0, Qmax, nq)
    j0 = np.sinc(np.outer(Qs, x) / np.pi)
    return Qs, 4 * np.pi * simpson(x[None, :]**2 * F[None, :]**n * j0, x=x, axis=1)


def pair_rate(bg, d, Q, c=1.0, M=MRED, nk=600, nc=201):
    Om = bg.Om
    Wp = (d - 2) * Om
    if Wp <= 2:
        return dict(open=False)
    Kmax = np.sqrt((Wp - 1)**2 - 1)
    Qs, Gt = Gt_table(bg, d - 2, 2 * Kmax + 1e-9)
    # Gauss-Legendre in K1 and cos
    xk, wk = np.polynomial.legendre.leggauss(nk)
    K1 = 0.5 * Kmax * (xk + 1); wK = 0.5 * Kmax * wk
    xc, wc = np.polynomial.legendre.leggauss(nc)
    e1 = np.sqrt(K1**2 + 1); e2 = Wp - e1
    K2 = np.sqrt(np.maximum(e2**2 - 1, 0))
    Qv = np.sqrt(np.maximum(K1[:, None]**2 + K2[:, None]**2 + 2 * K1[:, None] * K2[:, None] * xc[None, :], 0))
    G2 = np.interp(Qv, Qs, Gt)**2
    inner = (G2 * wc[None, :]).sum(axis=1)
    J = np.sum(wK * K1**2 * K2 / e1 * inner) / (16 * np.pi**3)
    lam = abs(c) / M**(d - 4)
    pref = 2 * comb(d, 2)**2 * lam**2 * VPHI**(2 * (d - 2)) * 2.0**(-(d - 2)) / MX**3  # GeV
    events = pref * J / HBAR
    # local-density variants
    x = np.linspace(0, min(bg.R, 70.0), 14001)
    F, S = bg.FS(x)
    beta_vac = np.sqrt(1 - 4 / Wp**2)
    meff2 = S**2 + 2 * AQ * F**2
    beta_loc = np.sqrt(np.clip(1 - 4 * meff2 / Wp**2, 0, None))
    vol_vac = 4 * np.pi * simpson(x**2 * F**(2 * (d - 2)) * beta_vac, x=x)
    vol_loc = 4 * np.pi * simpson(x**2 * F**(2 * (d - 2)) * beta_loc, x=x)
    # LDA: Gamma = C^2 lam^2 v^(2(d-2)) 2^-(d-2) m_X^-3 * vol * 1/(4 pi)
    lda = lambda vol: comb(d, 2)**2 * lam**2 * VPHI**(2 * (d - 2)) * 2.0**(-(d - 2)) / MX**3 * vol / (4 * np.pi) / HBAR
    out = dict(open=True, pair_energy_over_mX=Wp, K_max=float(Kmax), J=float(J),
               events_per_s_free_born=float(events), events_per_s_lda_vacuum=float(lda(vol_vac)),
               events_per_s_lda_local_mass=float(lda(vol_loc)))
    for k in ['free_born', 'lda_vacuum', 'lda_local_mass']:
        ev = out['events_per_s_' + k]
        out['frac_core_rate_' + k] = float((d - 2) * ev / Q)
    out['frac_times_age_free_born'] = out['frac_core_rate_free_born'] * T_REF
    return out


if __name__ == '__main__':
    import json, sys
    from profile import ground_state
    from forced import Bg
    Om, sol = ground_state(200000.)
    bg = Bg(sol, Om)
    res = {}
    for d in [5, 6]:
        r = pair_rate(bg, d, 200000.)
        res[str(d)] = r
        print(d, json.dumps(r, indent=1))
    json.dump(res, open(sys.argv[1], 'w'), indent=1)
