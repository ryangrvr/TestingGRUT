"""Axial fermion production: independent half-period check and the multi-period Pauli saturation [C].

Model (the track's): -(M/2) Nbar exp(i gamma5 phi/f) N, phi = S cos(m t) (no expansion, no back-reaction). After the
rotation each helicity h and momentum k is a two-level system H = [[p, M], [M, -p]], p = h k - mu5, mu5 = phidot/2f
= -A sin t (units m = 1, A = S/(2 f)).
Method (differs from the track's RK4): exact 2x2 exponential per step with the midpoint p (unitary to machine
precision), both helicities carried together through many half-periods, occupations read in the free basis at every
turning point (mu5 = 0).
Outputs: (1) the single half-period energy eps(A, M) against the track's axial_check.json and axial_tables.json;
(2) the cumulative energy fraction F_tot(n) = sum_h eps_h / (S^2/2) after n half-periods, i.e. what the track's fluid
assumes keeps growing as n F_1 (no Pauli blocking), for the benchmark channels at their slice amplitude.
Usage: python3 -B sk_axial.py -> ../results/sk_axial.json
       python3 -B sk_axial.py fsat -> ../results/sk_axial_Fsat_vs_A.json (saturation at A = 10, 20; 3, 6 for N1)
"""
import math, json, os, sys
import numpy as np
import sk_common as C


def vac(k, h, M):
    p = h * k
    E = np.sqrt(p * p + M * M)
    a, b = M * np.ones_like(k), -(E + p)
    n = np.sqrt(a * a + b * b)
    return np.stack([a / n, b / n]).astype(complex)


def occ(psi, k, h, M):
    p = h * k
    E = np.sqrt(p * p + M * M)
    a, b = M * np.ones_like(k), E - p
    n = np.sqrt(a * a + b * b)
    return np.abs((a * psi[0] + b * psi[1]) / n) ** 2


def evolve(psi, k, h, A, M, t0, t1, wdt=0.04):
    wmax = k.max() + A + M
    n = int(math.ceil((t1 - t0) * wmax / wdt))
    dt = (t1 - t0) / n
    for i in range(n):
        tm = t0 + (i + 0.5) * dt
        p = h * k + A * math.sin(tm)
        w = np.sqrt(p * p + M * M)
        c, s = np.cos(w * dt), np.sin(w * dt) / w
        a0 = c * psi[0] - 1j * s * (p * psi[0] + M * psi[1])
        a1 = c * psi[1] - 1j * s * (M * psi[0] - p * psi[1])
        psi = np.stack([a0, a1])
    return psi


def kgrid(A, M, nk):
    return np.linspace(1e-3, 1.25 * A + 3 * M + 2, nk)


def eps_of(n, k, M):
    return float(np.trapezoid(k ** 2 * n * np.sqrt(k * k + M * M), k) / (2 * math.pi ** 2))


def half_period(A, M, nk=1500, wdt=0.04):
    k = kgrid(A, M, nk)
    psi = evolve(vac(k, -1.0, M), k, -1.0, A, M, 0.0, math.pi, wdt)
    n = occ(psi, k, -1.0, M)
    return eps_of(n, k, M), float(n.max())


def multi(A, M, f_over_m, nhalf=16, nk=1500, wdt=0.04):
    """both helicities from the vacuum through nhalf half-periods; returns F_tot(n) and the band occupation."""
    k = kgrid(A, M, nk)
    S = 2 * A * f_over_m
    out = []
    psis = {h: vac(k, h, M) for h in (-1.0, 1.0)}
    for j in range(nhalf):
        t0, t1 = j * math.pi, (j + 1) * math.pi
        for h in (-1.0, 1.0):
            psis[h] = evolve(psis[h], k, h, A, M, t0, t1, wdt)
        ns = {h: occ(psis[h], k, h, M) for h in (-1.0, 1.0)}
        e = sum(eps_of(ns[h], k, M) for h in ns)
        band = k < A
        nbar = float(np.mean([ns[h][band].mean() for h in ns]))
        out.append(dict(nhalf=j + 1, eps_tot=e, F_tot=e / (0.5 * S * S), nbar_band=nbar))
    cap1 = 2 * A ** 4 / (8 * math.pi ** 2) / (0.5 * S * S)       # n = 1 in the band k < A, both helicities
    return dict(A=A, M_over_m=M, f_over_m=f_over_m, S=S, rows=out, F_band_capacity_n1=cap1,
                F_band_capacity_formula='2 (mu5max/4 pi f)^2 = %.4f' % (2 * (A / (4 * math.pi * f_over_m)) ** 2))


def fsat_scan():
    out = {}
    for nm, Mh, fm, As in [('N1', 250 / 600., 2., (3., 6.)), ('N2', 1.6, 8., (10., 20.)), ('N3', 0.8, 8., (10., 20.)),
                           ('N4', 0.48, 8., (10., 20.))]:
        for A in As:
            m = multi(A, Mh, fm, nhalf=16, nk=900)
            F = [r['F_tot'] for r in m['rows']]
            Fs = float(np.mean(F[6:]))
            out['%s_A%g' % (nm, A)] = dict(A=A, Mh=Mh, fm=fm, F1=F[0], Fsat=Fs, cap1=m['F_band_capacity_n1'],
                                           ratio=Fs / m['F_band_capacity_n1'])
            print(nm, A, out['%s_A%g' % (nm, A)], flush=True)
    json.dump(out, open(os.path.join(C.OUT, 'sk_axial_Fsat_vs_A.json'), 'w'), indent=1)


if __name__ == '__main__' and 'fsat' in sys.argv:
    fsat_scan()
elif __name__ == '__main__':
    res = dict(label='C (skeptic)', note=__doc__)
    tr = json.load(open(os.path.join(C.BENCH, 'axial_check.json')))['checks']
    chk = []
    for r in tr:
        e, nmx = half_period(r['A'], r['Mh'])
        e2, _ = half_period(r['A'], r['Mh'], nk=2500, wdt=0.02)
        chk.append(dict(A=r['A'], Mh=r['Mh'], eps_track=r['eps'], eps_mine=e, eps_mine_fine=e2, ratio=e / r['eps'],
                        nmax_track=r['nmax'], nmax_mine=nmx))
        print(chk[-1], flush=True)
    res['half_period_vs_track'] = chk
    # benchmark channels at the slice amplitude (track diagnostics.json: A_slice, M/m, f/m)
    bench = [('sigma->N1', 12.34, 250 / 600., 1200 / 600.), ('chi->N2', 42.3, 400 / 250., 2000 / 250.),
             ('chi->N3', 42.3, 200 / 250., 2000 / 250.), ('chi->N4', 42.3, 120 / 250., 2000 / 250.)]
    tabs = json.load(open(os.path.join(C.BENCH, 'axial_tables.json')))
    mres = {}
    for nm, A, Mh, fm in bench:
        m = multi(A, Mh, fm)
        # track's per-half-period fraction (its table) at this A
        key = min(tabs.keys(), key=lambda s: abs(float(s) - Mh))
        tA, te = np.array(tabs[key]['A']), np.array(tabs[key]['eps'])
        e1 = float(np.exp(np.interp(math.log(A), np.log(tA), np.log(te))))
        m['F1_track_table'] = e1 / (0.5 * m['S'] ** 2)
        m['track_unblocked_F_after_n'] = [r['nhalf'] * m['F1_track_table'] for r in m['rows']]
        mres[nm] = m
        print(nm, 'F1_track %.4f' % m['F1_track_table'], ['%.4f' % r['F_tot'] for r in m['rows']],
              'nbar', ['%.2f' % r['nbar_band'] for r in m['rows'][::3]], 'cap1 %.3f' % m['F_band_capacity_n1'], flush=True)
    res['multi_period'] = mres
    C.dump('sk_axial.json', res)
