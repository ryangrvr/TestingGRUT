"""Skeptic of track higgs_reheat, inflation side.  Own code (Claude, skeptic pass, 2 Nov 2026).

1. xi running, independently: own one-loop SM RGEs (g', g, g3, y_t, lambda; 1-loop only) from the Buttazzo et al.
   1307.3536 matching fits (not the SMDR grid the track uses); factor F = (xi(H_i) - 1/6)/(xi(M_t) - 1/6).
2. The curvature mass the track calls 'one loop' (0.937 H^2): decomposed into 12 xi(mu) e^{2 Gamma} (tree, M_t-normalised
   field) and the genuine one-loop piece, using smpot.V_dS.
3. Super-horizon variance of the canonical field (Bunch-Davies, scipy Hankel, own quadrature) for m^2 = 12 xi
   (canonical) and for the track's m^2 = 0.937 (which mixes normalisations); light-field formula in both units.
4. Hawking-Moss exponent, tree RG-improved, own barrier search on V = xi(mu) R phi^2/2 + lambda(mu) phi^4/4,
   mu^2 = h^2 + R, couplings and Gamma from smpot (SM input); xi_min (tree) by brentq.
Writes ../results/sk_infl.json.
"""
import json
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import hankel1

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'code')))
import smpot as S  # noqa: E402

K16 = 16 * np.pi ** 2


def buttazzo_match(Mt, Mh=125.13, als=0.1180):
    a = (als - 0.1184) / 0.0007
    lam = 0.12604 + 0.00206 * (Mh - 125.15) - 0.00004 * (Mt - 173.34)
    yt = 0.93690 + 0.00556 * (Mt - 173.34) - 0.00042 * a
    g3 = 1.1666 + 0.00314 * a - 0.00046 * (Mt - 173.34)
    g2 = 0.64779 + 0.00004 * (Mt - 173.34)
    gY = 0.35830 + 0.00011 * (Mt - 173.34)
    return np.array([gY, g2, g3, yt, lam])


def rge1(t, y):
    gp, g, g3, yt, lam = y[:5]
    d = np.empty(6)
    d[0] = 41 / 6 * gp ** 3 / K16
    d[1] = -19 / 6 * g ** 3 / K16
    d[2] = -7 * g3 ** 3 / K16
    d[3] = yt * (9 / 2 * yt ** 2 - 17 / 12 * gp ** 2 - 9 / 4 * g ** 2 - 8 * g3 ** 2) / K16
    d[4] = (24 * lam ** 2 - 6 * yt ** 4 + 3 / 8 * (2 * g ** 4 + (g ** 2 + gp ** 2) ** 2)
            + lam * (12 * yt ** 2 - 9 * g ** 2 - 3 * gp ** 2)) / K16
    d[5] = (12 * lam + 6 * yt ** 2 - 1.5 * gp ** 2 - 4.5 * g ** 2) / K16      # gamma_xi
    return d


def xi_factor_indep(Mt, mu):
    y0 = np.concatenate([buttazzo_match(Mt), [0.0]])
    sol = solve_ivp(rge1, (0, np.log(mu / Mt)), y0, rtol=1e-10, atol=1e-12)
    return float(np.exp(sol.y[5, -1])), sol.y[:5, -1].tolist()


def ds_var(m2, kmax=1.0, kmin=1e-30):
    """per-component super-horizon variance (units H^2) of a canonical BD field, k < kmax aH."""
    nu = np.sqrt(9 / 4 - m2)

    def P(lk):
        k = np.exp(lk)
        return k ** 3 * abs(hankel1(nu, k)) ** 2 / (8 * np.pi)
    # small-k asymptote P ~ k^{3-2nu}; integrate analytically below 1e-6
    lo = np.log(1e-6)
    I1 = quad(P, lo, np.log(kmax), limit=400)[0]
    Pk = P(lo)
    I0 = Pk / (3 - 2 * nu)                                # tail below 1e-6 as a power law
    return I1 + I0


def hm_tree(sm, H, xi0, mu0=None):
    """B = 8 pi^2 Delta V_b /(3 H^4), tree RG-improved, canonical phi = e^Gamma h0 with mu^2 = h0^2 + R."""
    R = 12 * H ** 2

    def V(lh):
        h = np.exp(lh)
        mu = np.sqrt(h * h + R)
        c = sm.cpl(mu)
        phi = np.exp(c['Gam']) * h
        xi = sm.xi(mu, xi0, mu0)
        return float(0.5 * xi * R * phi ** 2 + 0.25 * c['lam'] * phi ** 4)
    lhs = np.linspace(np.log(H), np.log(1e16), 800)
    Vs = np.array([V(x) for x in lhs])
    k = int(np.argmax(Vs))
    if k == 0 or k == len(lhs) - 1:
        return None, None
    r = minimize_scalar(lambda x: -V(x), bracket=(lhs[k - 1], lhs[k], lhs[k + 1]), tol=1e-10)
    return float(np.exp(r.x)), 8 * np.pi ** 2 * (-r.fun) / (3 * H ** 4)


def main():
    out = {}
    LN49 = np.log(4 * np.pi / 9)
    for Mt in (172.33, 172.60, 172.87):
        sm = S.SM(Mt)
        row = {}
        for H in (1.29e11, 1.51e11):
            F, cpl = xi_factor_indep(Mt, H)
            c, c0 = sm.cpl(H), sm.cpl(Mt)
            Ftrack = float(np.exp(c['Ixi'] - c0['Ixi']))
            xiH = 1 / 6 + (0.1 - 1 / 6) * F
            # curvature mass decomposition
            h = 1e-3 * H
            m2_loop = float(2 * sm.V_dS(np.array([h]), H, 0.1)[0] / h ** 2) / H ** 2
            m2_tree_h0 = float(2 * sm.V_dS(np.array([h]), H, 0.1, loop=False)[0] / h ** 2) / H ** 2
            cR = sm.cpl(np.sqrt(12) * H)
            m2_canon_tree = 12 * float(sm.xi(np.sqrt(12) * H, 0.1))
            e2G = float(np.exp(2 * cR['Gam']))
            m2_canon_loop = m2_loop / e2G
            # variances
            v_can = ds_var(m2_canon_tree)
            v_trk = ds_var(m2_loop)
            light_can = 3 / (8 * np.pi ** 2 * m2_canon_tree)
            # tree HM
            hb, B = hm_tree(sm, H, 0.1)
            target = 3 * 60 + LN49
            f = lambda x: (hm_tree(sm, H, x)[1] or 0.0) - target
            try:
                xmin = brentq(f, 0.0, 0.1, xtol=1e-5)
            except ValueError:
                xmin = None
            row[f'H={H:.3g}'] = dict(
                xi_factor_indep_1loop_Buttazzo=F, xi_factor_track=Ftrack, xi_H_indep=xiH,
                couplings_at_H_indep=dict(zip(['gp', 'g', 'g3', 'yt', 'lam'], cpl)),
                m2_track_one_loop=m2_loop, m2_tree_in_h0_units_12xi_e2Gamma=m2_tree_h0,
                m2_genuine_one_loop_shift=m2_loop - m2_tree_h0, m2_canonical_tree=m2_canon_tree,
                m2_canonical_with_loop=m2_canon_loop, e2Gamma_at_sqrtR=e2G,
                var_superhorizon_canonical_phi=v_can, var_superhorizon_in_h0_units=v_can / e2G,
                var_track_value_check_m2_0937=v_trk, light_formula_canonical=light_can,
                light_formula_h0_units=light_can / e2G,
                HM_tree_indep=dict(h_b_over_H=(hb / H if hb else None), B=B), xi_EW_min_tree_indep_N60=xmin)
            print(Mt, H, json.dumps(row[f'H={H:.3g}'], default=float), flush=True)
        out[f'Mt={Mt}'] = row
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'sk_infl.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
