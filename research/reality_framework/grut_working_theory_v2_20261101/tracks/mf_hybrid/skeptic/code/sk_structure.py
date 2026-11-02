"""Independent check of the track's structural claims S1/S2 [D/C].

S1 claims: a heavy waterfall (m_chi >~ 3e3 H) forces m_sig <~ 13 H (X = g^2 sigma_c <~ 1) by one-loop flatness of the clock.
Here: my own e-fold capacity N_max(m_chi, m_sig, kappa) with the same definition (tree c1, m2 solved at each trial pivot
for A_s and n_s = 0.9649 at first order; N = int V/(mP^2 V'); the pivot must be reachable with V' > 0 throughout), but
with my own potential code, a log-spaced pivot scan up to 1 M_P, and a finer grid in m_sig, so that the largest m_sig
passing N_max >= 52.6 is found at m_chi = 3000 and 4000 H.
S2: perturbative bound re-derived for several decay types (scalar Yukawa Majorana/Dirac, trilinear to scalars, N_f species).
Usage: python3 -B sk_structure.py -> ../results/sk_structure.json
"""
import math
import numpy as np
import sk_common as C

mP = C.mP
AS = 2.1e-9


def params(m_chi, m_sig, kappa):
    lam = m_chi ** 4 / (48 * mP ** 2)
    v = m_chi / math.sqrt(2 * lam)
    g = m_sig / v
    sc = math.sqrt(lam) * v / g
    P = dict(lam=lam, v=v, g=g, sc=sc, cb=0.0, V0=0.25 * lam * v ** 4)
    P['mu2'] = kappa * g * g * (2 * sc) ** 2
    return P


def loop(P, s, n=0):
    """chi one-loop term and its first two sigma derivatives (analytic; independent re-derivation)."""
    m2 = C.M2(P, s)
    if m2 == 0.0:
        return 0.0
    L = math.log(abs(m2) / P['mu2'])
    a1, a2 = 2 * P['g'] ** 2 * s, 2 * P['g'] ** 2
    k = 64 * math.pi ** 2
    if n == 0:
        return m2 * m2 * (L - 1.5) / k
    if n == 1:
        return 2 * m2 * a1 * (L - 1.0) / k
    return 2 * (a1 * a1 * L + m2 * a2 * (L - 1.0)) / k


REGULAR = 0.0      # if > 0: require V' >= REGULAR * V'_* between sigma_c and the pivot


def nefolds(P, dl, ns=0.9649):
    eps_t = P['V0'] / (24 * math.pi ** 2 * AS * mP ** 4)
    eta_t = 0.5 * (ns - 1)
    ss = P['sc'] + dl * mP
    V1 = math.sqrt(2 * eps_t) * P['V0'] / mP
    V2 = eta_t * P['V0'] / mP ** 2
    m2 = loop(P, ss, 2) - V2
    c1 = V1 - loop(P, ss, 1) + m2 * (ss - P['sc'])
    grid = P['sc'] + (ss - P['sc']) * np.linspace(1e-4, 1.0, 1500) ** 2
    dV = np.array([c1 - m2 * (z - P['sc']) + loop(P, z, 1) for z in grid])
    if np.any(dV <= 0):
        return 0.0
    if REGULAR and dV.min() < REGULAR * V1:
        return 0.0          # a near-inflection between pivot and sigma_c (P_zeta spike > A_s/REGULAR^2): rejected
    Vv = np.array([P['V0'] + c1 * (z - P['sc']) - 0.5 * m2 * (z - P['sc']) ** 2 + loop(P, z) for z in grid])
    return float(np.trapezoid(Vv / (mP ** 2 * dV), grid))


def nmax(m_chi, m_sig, kappa, dlmax=0.2, n=100):
    """largest e-fold capacity on the small-field branch (pivot within dlmax M_P of sigma_c, as in the track's scan;
    beyond ~0.3 M_P the loop term itself becomes comparable to V0 and the V0-dominated construction is not valid)."""
    P = params(m_chi, m_sig, kappa)
    dls = np.geomspace(1e-4, dlmax, n)
    Ns = [nefolds(P, x) for x in dls]
    i = int(np.argmax(Ns))
    ss = P['sc'] + dls[i] * mP
    return float(Ns[i]), float(dls[i]), P, float(abs(loop(P, ss)) / P['V0'])


def kappa_scan(m_chi, m_sig, kappas=np.geomspace(0.1, 1e3, 21)):
    best = None
    rows = []
    for kap in kappas:
        N, dl, P, lf = nmax(m_chi, m_sig, kap)
        rows.append(dict(kappa=float(kap), Nmax=N, dl=dl, loop_over_V0=lf))
        if N >= 52.6 and lf < 0.05 and (best is None or N > best['Nmax']):
            best = rows[-1]
    npass = sum(1 for r in rows if r['Nmax'] >= 52.6 and r['loop_over_V0'] < 0.05)
    return dict(m_chi=m_chi, m_sig=m_sig, X=P['g'] ** 2 * P['sc'], sigma_c_MP=P['sc'] / mP,
                n_kappa_pass=npass, n_kappa=len(rows), best=best, rows=rows)


if __name__ == '__main__':
    rows = []
    for (mx, ms) in [(250, 600), (300, 750), (4000, 3), (4000, 10), (4000, 30), (4000, 90)]:
        r = dict(m_chi=mx, m_sig=ms)
        for kap in (1.0, 10.0, 100.0):
            N, dl, P, lf = nmax(mx, ms, kap)
            r['Nmax_k%g' % kap] = N
            r['dl_k%g' % kap] = dl
        r['X'] = P['g'] ** 2 * P['sc']
        r['sigma_c_MP'] = P['sc'] / mP
        r['passes_52.6_all_kappa'] = all(r['Nmax_k%g' % k] >= 52.6 for k in (1.0, 10.0, 100.0))
        rows.append(r)
        print(r, flush=True)
    ks = []
    REGULAR = 0.1
    for (mx, ms) in [(4000, 13), (4000, 30), (4000, 60), (4000, 90), (4000, 150), (4000, 300), (3000, 60), (3000, 150)]:
        ks.append(kappa_scan(mx, ms))
        ks[-1]['regularity'] = 'V prime >= 0.1 V prime_* between sigma_c and the pivot (P_zeta < 100 A_s)'
        print(mx, ms, 'X %.2f' % ks[-1]['X'], 'pass %d/%d' % (ks[-1]['n_kappa_pass'], ks[-1]['n_kappa']), ks[-1]['best'], flush=True)
    # S2 perturbative bound variants: amplitude A = sqrt(6) H mP/m, Gamma >= H
    b = {}
    b['Majorana_scalar_yukawa (track)'] = (384 * math.pi) ** 0.2 * mP ** 0.4
    b['Dirac_scalar_yukawa'] = (192 * math.pi) ** 0.2 * mP ** 0.4
    b['Majorana_Nf=10'] = (384 * math.pi / 10) ** 0.2 * mP ** 0.4
    b['trilinear_to_complex_scalar (A phi |h|^2, A A_amp < m^2/4)'] = (768 * math.pi / 4) ** 0.2 * mP ** 0.4
    C.dump('sk_structure.json', dict(label='C/D (skeptic)', note=__doc__, rows=rows, kappa_scans=ks, perturbative_bound_m_over_H=b))
