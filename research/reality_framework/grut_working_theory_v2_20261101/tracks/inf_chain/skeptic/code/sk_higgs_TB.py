#!/usr/bin/env python3
"""inf_chain SKEPTIC (Claude, own code): Higgs survival DURING INFLATION at the quiet T-models' own Hubble scale
(TB1 H_i = 7.5e12, H_* = 9.0e12; TB2 H_i = 1.06e13, H_* = 1.6e13 GeV), which the track did not recompute ('Higgs
survival at the new H_i not recomputed').  inf_higgs computed it only at H_i = 1.29-1.50e11 GeV.

Criterion: inf_higgs's (Mantziris-Markkanen-Rajantie arXiv:2011.03763 form), re-implemented here:
    <N> = (4 pi/3) int dN_b e^{3(N_0 - N_b)} (H(N_b)/H_p)^3 e^{-B(N_b)},  B = 8 pi^2 DeltaV_b/(3 H^4),
    survival if <N> < 1; N_0 = N_obs + 1 (exit of today's horizon + 1);
    tree RG-improved barrier V = xi(mu) R phi^2/2 + lam(mu) phi^4/4, phi = e^{Gamma} h, mu^2 = h^2 + R,
    R = 6(2 - eps_H) H^2 along the benchmark's own H(N) (my sk_cmb_fluid background).
SM running: higgs_reheat/code/smpot.py (3-loop, imported read-only; PDG 2024 M_h = 125.20, alpha_s = 0.1180).
Regression: HT6t at H_i = 1.29e11, N_0 = 55.76 against inf_higgs inflation_HN.json (tree MD xi_min_HN).
Output ../results/sk_higgs_TB.json.  Labels: [C] own code on a model criterion [E].
"""
import json, math, sys, pathlib
import numpy as np
from scipy.optimize import brentq
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sk_cmb_fluid as sk
TRK = HERE.parent.parent.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(TRK / 'higgs_reheat' / 'code'))
import smpot as S           # read-only

MH, ALS = 125.20, 0.1180
H0C = 67.4 / 299792.458     # Mpc^-1


def history(m, f0, H_i_GeV, Nmax=4.0):
    sol, Nend = sk.background(m, f0)
    f, fN = sol.sol(Nend - 1.0); hS1, _ = sk.H_of(m, f, fN)
    sc = H_i_GeV / hS1
    Nb = np.unique(np.concatenate([np.linspace(0, 0.3, 16), np.linspace(0.3, 1.2, 13), np.linspace(1.2, Nmax, 10)]))
    rows = []
    for x in Nb:
        f, fN = sol.sol(Nend - x); h, eps = sk.H_of(m, f, fN)
        rows.append((x, h * sc, 6 * (2 - eps) * (h * sc) ** 2, eps))
    f, fN = sol.sol(Nend - 50.0); hp, _ = sk.H_of(m, f, fN)
    return dict(rows=np.array(rows), Hp=hp * sc, sol=sol, Nend=Nend, sc=sc, hS1=hS1)


def N_obs(hist, Nstar):
    """e-folds before eps_H = 1 at which k_0 = a_0 H_0 exited, given the pivot exit at Nstar."""
    sol, Nend = hist['sol'], hist['Nend']
    aH = lambda Nb: math.exp(-Nb) * sk.H_of(TMP['m'], *sol.sol(Nend - Nb))[0]
    k0 = aH(Nstar) * H0C / 0.05
    return brentq(lambda Nb: aH(Nb) - k0, Nstar, Nstar + 12)


TMP = {}


def barrier_B(sm, H, R, xi_EW, n=700):
    hs = np.exp(np.linspace(math.log(1e-2 * H), math.log(1e17), n))
    def V(h):
        mu = np.sqrt(h ** 2 + R)
        c = sm.cpl(mu); xi = sm.xi(mu, xi_EW)
        phi = np.exp(c['Gam']) * h
        return 0.5 * xi * R * phi ** 2 + 0.25 * c['lam'] * phi ** 4
    Vv = V(hs); k = int(np.argmax(Vv))
    if k == 0 or Vv[k] <= 0:
        return 0.0, None
    if k == n - 1:
        return np.inf, None
    hh = np.exp(np.linspace(math.log(hs[k - 1]), math.log(hs[k + 1]), 151))
    VV = V(hh); j = int(np.argmax(VV))
    return 8 * math.pi ** 2 * float(VV[j]) / (3 * H ** 4), float(hh[j])


def log10N(hist, sm, xi, N0):
    r = hist['rows']; Nb, H, R = r[:, 0], r[:, 1], r[:, 2]
    Bs = np.array([barrier_B(sm, H[i], R[i], xi)[0] for i in range(len(Nb))])
    Hp = hist['Hp']; Bp = barrier_B(sm, Hp, 12 * Hp ** 2, xi)[0]
    lnI = 3 * (N0 - Nb) + 3 * np.log(H / Hp) - Bs
    mx = np.max(lnI)
    part = mx + math.log(np.trapezoid(np.exp(lnI - mx), Nb))
    tail = 3 * (N0 - Nb[-1]) - Bp - math.log(3.0)
    tot = np.logaddexp(part, tail) + math.log(4 * math.pi / 3)
    return float(tot / math.log(10)), dict(B_end=float(Bs[0]), B_plateau=float(Bp))


def xi_min(f, lo=0.0, hi=0.6):
    xs = np.linspace(lo, hi, 13); fs = [f(x) for x in xs]
    for i in range(len(xs) - 1):
        if fs[i] > 0 >= fs[i + 1]:
            return float(brentq(f, xs[i], xs[i + 1], xtol=2e-5))
    return None


def main():
    out = dict(method=__doc__, rows={})
    cmb = json.load(open(HERE.parent / 'results' / 'sk_cmb_fluid.json'))
    cases = [('HT6t', sk.Hilltop(0.3469421539393091, 6, 1.5529238776548362e-05), 0.02, 1.29e11, None, 55.76),
             ('TB1', sk.TModel(0.3), None, cmb['TB1']['H_i_GeV'], cmb['TB1']['Nstar'], None),
             ('TB2', sk.TModel(1.0), None, cmb['TB2']['H_i_GeV'], cmb['TB2']['Nstar'], None)]
    sms = {Mt: S.SM(Mt, Mh=MH, als=ALS) for Mt in (172.28, 172.57, 172.86)}
    ih = json.load(open(TRK / 'inf_higgs' / 'results' / 'inflation_HN.json'))['HT6t|iH0']
    for key, m, f0, Hi, Nst, N0 in cases:
        TMP['m'] = m
        if f0 is None:
            f0 = m.phi_start(80.0)
        hist = history(m, f0, Hi)
        if N0 is None:
            N0 = N_obs(hist, Nst) + 1.0
        row = dict(H_i_GeV=Hi, H_plateau_GeV=hist['Hp'], N0=N0, by_Mt={})
        for Mt, sm in sms.items():
            xm = xi_min(lambda x: log10N(hist, sm, x, N0)[0])
            l01, info = log10N(hist, sm, 0.1, N0)
            c = sm.cpl(10 * hist['Hp'])
            row['by_Mt'][str(Mt)] = dict(xi_EW_min=xm, log10N_at_xi01=l01, lam_at_10Hp=float(c['lam']), **info)
            if key == 'HT6t':
                row['by_Mt'][str(Mt)]['inf_higgs_tree_MD_xi_min_HN'] = ih[str(Mt)]['tree MD xi_min_HN']
            print(key, Mt, row['by_Mt'][str(Mt)], flush=True)
        out['rows'][key] = row
    json.dump(out, open(HERE.parent / 'results' / 'sk_higgs_TB.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
