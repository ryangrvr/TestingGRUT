"""Structural results for the non-supersymmetric hybrid with a heavy waterfall [D + C].

1. Mass relations (V0-dominated, lam v^4/4 = 3 M_P^2 H^2):  v = sqrt(24) M_P/m_chi,  sigma_c = 2 sqrt(3) M_P/m_sig,
   X = g^2 sigma_c (half the slope of M_chi^2 at sigma_c) = m_sig m_chi^2 / (4 sqrt(3) M_P)   [D].
2. One-loop flatness: with the full chi loop, the clock (tree c1, m2 retuned for A_s and n_s at each pivot) can supply
   the needed e-folds only while g^2 sigma_* is O(1) H. N_max(n_s) is the largest number of e-folds before sigma_c at
   which A_s and n_s can both be met with a monotonic potential [C]. Renormalisation scale mu^2 = kappa M_chi^2(2 sigma_c)
   (equivalently the tree quartic) is scanned.
3. Perturbative reheating bound: a field holding V0 with mass m oscillates with amplitude sqrt(6) H M_P/m; a channel
   with y A < m/2 has Gamma <= m^3/(64 pi A^2); Gamma >= H needs m >= (384 pi)^(1/5) (M_P/H)^(2/5) H   [D].
Usage: python3 -B structure.py -> ../results/structure.json
"""
import math, json, os
import numpy as np
import mfh_common as mc
import inflation as I

HERE = os.path.dirname(os.path.abspath(__file__))


def Nmax(p, mu2, eta_t):
    mP, V0 = p['mP'], p['V0']
    eps_t = V0 / (24 * math.pi ** 2 * mc.AS * mP ** 4)

    def N(dl):
        ss = p['sigma_c'] + dl * mP
        base = mc.Inflation(p, 0., 0., mu2=mu2)
        l1, l2, _ = base.derivs(ss)
        V1 = math.sqrt(2 * eps_t) * V0 / mP
        V2 = eta_t * V0 / mP ** 2
        m2 = l2 - V2
        c1 = V1 - l1 + m2 * (ss - p['sigma_c'])
        inf = mc.Inflation(p, c1, m2, mu2=mu2)
        xs = np.linspace(p['sigma_c'] * (1 + 1e-9), ss, 150)
        if min(inf.dV(z) for z in xs) <= 0:
            return 0.0
        return I.N_of_sigma(inf, ss)
    ds = np.geomspace(1e-3, 0.2, 12)
    Ns = [N(x) for x in ds]
    i = int(np.argmax(Ns))
    return float(Ns[i]), float(ds[i])


if __name__ == '__main__':
    H = 1.30e11
    mP = mc.MP_GEV / H
    rows = []
    for (mx, ms) in [(250, 600), (300, 600), (300, 750), (1000, 100), (4000, 3), (4000, 10), (4000, 30)]:
        p = mc.derived(H, float(mx), float(ms))
        r = dict(m_chi=mx, m_sig=ms, X=p['X'], sigma_c_MP=p['sigma_c_MP'], v_MP=p['v_MP'], lam=p['lam'], g=p['g'],
                 X_formula=ms * mx ** 2 / (4 * math.sqrt(3) * mP))
        for kap in (1.0, 10.0, 100.0):
            mu2 = kap * p['g'] ** 2 * (2 * p['sigma_c']) ** 2
            r['Nmax_kappa%g' % kap] = Nmax(p, mu2, 0.5 * (mc.NS_PLANCK - 1))
        rows.append(r)
        print(r, flush=True)
    mpert = (384 * math.pi) ** 0.2 * mP ** 0.4
    out = dict(label='D (relations) / C (N_max scan at n_s = 0.9649, H = 1.3e11 GeV)', note=__doc__, rows=rows,
               N_needed='N_* = 52.6 on the benchmarks (K = 2.05)',
               m_min_perturbative_fast_over_H=mpert,
               quiet_clock='for heavy chi (lam >~ 6e-3, m_chi >~ 3000 H) X <~ 1 gives m_sig = X/sqrt(lam) <~ 13 H: the clock is '
                           'quiet and holds V0 on the adiabatic valley (G4)')
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'structure.json'), 'w'), indent=1)
