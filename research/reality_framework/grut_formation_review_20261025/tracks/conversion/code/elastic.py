#!/usr/bin/env python3
"""Post-escape elastic scattering of X on the frozen Q = 200000 object (Claude review, track `conversion`).

After escape (E_X < 1 + 2 Omega) the anti-X channel is closed, so the only open channel is
elastic X scattering.  The v and s channels are closed, and their growing parts are removed
at the matching radius.  Same radial integrator as scatter.py (own code).
  sigma_T = (4 pi / k^2) sum_l (l+1) sin^2(delta_l - delta_{l+1})   (momentum transfer)
  sigma_el = (4 pi / k^2) sum_l (2l+1) sin^2 delta_l
Then Gamma_T/H = n_obj sigma_T v / H along the post-escape redshift of each cohort, the
elastic decoupling scale factor (Gamma_T/H = 1), and the collisionless horizon counted
from decoupling.  This is a sensitivity estimate for the frozen-object model.  It is not a
transport solution with diffusion.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import spherical_jn, spherical_yn, spherical_in, spherical_kn

sys.path.insert(0, str(Path(__file__).resolve().parent))
import scatter as S
import transport_mc as T


def phase_shift(ell, E, rho_m=40.0, rtol=1e-10):
    """single open (u) channel: E = X energy / m_X, nu = E - Omega."""
    OM = S.OM
    nu = E - OM
    K2 = np.array([E ** 2 - 1, (nu - OM) ** 2 - 1, nu ** 2 - 6])
    assert K2[0] > 0 and K2[1] < 0 and K2[2] < 0
    Y = S.propagate(ell, K2, rho_m, rtol)
    U, Up = Y[:3], Y[3:]
    dl = []
    for i in (1, 2):
        kap = np.sqrt(-K2[i]); x = kap * rho_m
        ih = x * spherical_in(ell, x); ihp = kap * (spherical_in(ell, x) + x * spherical_in(ell, x, True))
        kh = x * spherical_kn(ell, x); khp = kap * (spherical_kn(ell, x) + x * spherical_kn(ell, x, True))
        Wc = kh * ihp - khp * ih
        dl.append((kh * Up[i] - khp * U[i]) / Wc)
    _, _, vt = np.linalg.svd(np.array(dl))
    c = vt[-1]                       # null vector (both growing parts removed)
    u, up = U[0] @ c, Up[0] @ c
    k = np.sqrt(K2[0]); x = k * rho_m
    jh = x * spherical_jn(ell, x); jhp = k * (spherical_jn(ell, x) + x * spherical_jn(ell, x, True))
    nh = x * spherical_yn(ell, x); nhp = k * (spherical_yn(ell, x) + x * spherical_yn(ell, x, True))
    W = jh * nhp - jhp * nh
    A = (u * nhp - up * nh) / W; Bc = (jh * up - jhp * u) / W
    # U = A jhat + B nhat ~ A sin - B cos  ->  tan(delta) = -B/A
    return float(np.arctan2(-Bc, A))


def cross_sections(p):
    E = np.sqrt(1 + p * p); k = p
    deltas = []
    ell = 0
    while True:
        deltas.append(phase_shift(ell, E))
        if ell > 6 and abs(np.sin(deltas[-1])) < 1e-6 and abs(np.sin(deltas[-2])) < 1e-6:
            break
        ell += 1
        if ell > 60:
            break
    d = np.array(deltas + [0.0])
    sT = 4 * np.pi / k ** 2 * np.sum((np.arange(len(d) - 1) + 1) * np.sin(d[:-1] - d[1:]) ** 2)
    sel = 4 * np.pi / k ** 2 * np.sum((2 * np.arange(len(d) - 1) + 1) * np.sin(d[:-1]) ** 2)
    return float(sT), float(sel), len(deltas) - 1


def main():
    ps = np.r_[0.01, 0.03, 0.1, 0.2, 0.4, 0.7, 1.0, 1.4, 1.8, 2.2, 2.5]
    rows = []
    for p in ps:
        sT, sel, lmax = cross_sections(p)
        rows.append(dict(p_over_m=float(p), sigma_T_mX2=sT, sigma_el_mX2=sel, ellmax=lmax,
                         sigma_T_cm2=sT * (S.HBARC_GEV_CM / S.M_X) ** 2))
        print(rows[-1], flush=True)
    lp = np.log([r['p_over_m'] for r in rows]); ls = np.log([r['sigma_T_cm2'] for r in rows])
    sig = lambda p: np.exp(np.interp(np.log(np.clip(p, ps[0], ps[-1])), lp, ls))
    out = dict(rows=rows)
    # decoupling for the three worked sources: cohort at the median escape
    res = {}
    for k, a_esc in [('T20MeV', 9.227617887109783e-12), ('T1keV', 2.6014e-07), ('z1e6', 1.1243e-06)]:
        q = T.UESC * a_esc

        def gam(a):
            p = q / a; v = p / np.sqrt(1 + p * p)
            return T.NOBJ0 / a ** 3 * sig(p) * v * 2.99792458e10 / T.hubble(a)
        g0 = gam(a_esc)
        a_dec = brentq(lambda a: gam(a) - 1, a_esc, 1.0) if g0 > 1 else a_esc
        res[k] = dict(a_escape=a_esc, GammaT_over_H_at_escape=g0, a_elastic_decoupling=a_dec,
                      z_elastic_decoupling=1 / a_dec - 1,
                      horizon_from_escape_Mpc=T.horizon(a_esc, q, 1.0),
                      horizon_from_decoupling_Mpc=T.horizon(a_dec, q, 1.0),
                      v_over_c_at_decoupling=float((q / a_dec) / np.sqrt(1 + (q / a_dec) ** 2)))
        print(k, res[k], flush=True)
    out['decoupling'] = res
    (T.RES / 'elastic.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
