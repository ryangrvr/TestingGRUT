#!/usr/bin/env python3
"""dq_heat SKEPTIC (Claude, own code): vector heat for explicit inflaton paths before the slice, mode by mode.

Each mode k < 1 left the horizon N(k) e-folds before the slice (k = e^{-N} H(N)/H_i).  Its frozen Goldstone power is
(H(N)/2 pi F(N))^2, so dq_heat's per-mode spectrum (static slice) is multiplied by
    f(N) = [ h_r(N) F(1,1) / F(s(N), h_r(N)) ]^2,
with F re-solved at the exit state (u = mu + s (10i - mu), R = 12 h_r^2; sk_exit.py).  Modes k >= 1 (release era) keep f = 1.
Paths [E] (M1r supplies none; m_I = 3 H_i is not a slow-roll mass, so these are embeddings, not M1r):
  static     : s = h_r = 1 (M1r as specified);
  fixedH_lin : s = 1 + 0.2 N, h_r = 1 (holder drifts, H fixed: dq_heat's family);
  phi2       : quadratic large field, slice at epsilon = 1: phi^2 = 2 + 4N -> s = h_r = (1 + 2N)^{1/2};
  plateau    : Starobinsky-type V0 (1 - e^{-b phi})^2, b = (2/3)^{1/2}, slice where m = 3H (e^{-b phi_s} = 1/3, epsilon = 1/3),
               slow-roll N(phi), s = phi/phi_s, h_r = (1 - e^{-b phi})/(2/3).
R_total(e) = band + SU + mediators + cascade (unchanged) + f-weighted piece A; S/Q from dq_entropy's closed form.
Output: ../results/sk_exit_path.json.   Labels: [C] given the path, [E] for the paths.
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_model_copy as fc

HERE = Path(__file__).resolve().parent
DQ = HERE.parents[1] / 'results'
RES = HERE.parent / 'results'
Q = np.array([1.0, 3.0, 6.0])
B = math.sqrt(2.0 / 3.0)
PHIS = math.log(3.0) / B


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def plateau_N(phi):
    g = lambda x: (math.exp(B * x) / B - x) / (2 * B)
    return g(phi) - g(PHIS)


def path(name, N):
    if name == 'static':
        return 1.0, 1.0
    if name == 'fixedH_lin':
        return 1.0 + 0.2 * N, 1.0
    if name == 'phi2':
        s = math.sqrt(1 + 2 * N); return s, s
    if name == 'plateau':
        phi = brentq(lambda x: plateau_N(x) - N, PHIS, 50.0) if N > 0 else PHIS
        return phi / PHIS, (1 - math.exp(-B * phi)) / (2.0 / 3.0)
    raise ValueError(name)


def main():
    p = fc.params(0.5)
    ref = fc.stationary(fc.U_PREP, p, R=12.0)
    F11 = float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(ref['z']) ** 2)))
    vm = json.load(open(DQ / 'vector_modes.json'))
    hb = json.load(open(DQ / 'heat_table.json'))['budgets']['portals_vanish_at_MY_benchmark']
    W0 = hb['W0']; rows_e = {r['e']: r for r in hb['rows']}
    Ng = np.linspace(0.0, 6.0, 49)
    out = dict(method=__doc__, F11=F11, paths={})
    for name in ('static', 'fixedH_lin', 'phi2', 'plateau'):
        fN = []; sN = []; hN = []
        guess = None
        for N in Ng:
            s, hr = path(name, N)
            u = p['mu'] + (fc.U_PREP - p['mu']) * s
            st = fc.stationary(u, p, R=12.0 * hr * hr, guess=guess)
            if st['grad'] > 1e-8:
                st = fc.stationary(u, p, R=12.0 * hr * hr, guess=None)
            guess = st['z']
            F = float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(st['z']) ** 2)))
            fN.append((hr * F11 / F) ** 2); sN.append(s); hN.append(hr)
        fN = np.array(fN); hN = np.array(hN)
        # exit e-fold of mode k: k = e^{-N} h_r(N)  (H_i units, a = 1 at the slice)
        lnk_of_N = -Ng + np.log(hN)

        def f_of_k(k):
            if k >= hN[0]:
                return 1.0
            N = float(np.interp(math.log(k), lnk_of_N[::-1], Ng[::-1]))
            return float(np.interp(N, Ng, fN))
        res = dict(N=Ng.tolist(), s=sN, h_r=hN.tolist(), f=fN.tolist(), by_e=[])
        for c in vm['cases']:
            if c['variant'] != 'full':
                continue
            e = c['e']; r = rows_e[e]
            ks = np.array([m['k'] for m in c['spectrum']]); dn = np.array([m['dN_dlnk'] for m in c['spectrum']])
            lo = ks <= 1.0
            fk = np.array([f_of_k(k) for k in ks])
            NA_lo_static = float(np.trapezoid(dn[lo], np.log(ks[lo])))
            NA_lo_path = float(np.trapezoid((dn * fk)[lo], np.log(ks[lo])))
            f_eff = NA_lo_path / NA_lo_static
            Rt = r['R_total'] - r['R_vec_k_le_1'] + f_eff * r['R_vec_k_le_1']
            res['by_e'].append(dict(e=e, f_eff=f_eff, R_total=Rt, S_per_Q=SQ(W0 * (1 + Rt))))
        es = sorted(x['e'] for x in res['by_e'] if x['e'] >= 1e-7)
        le = np.log(es); lW = [math.log(W0 * (1 + [x for x in res['by_e'] if x['e'] == e][0]['R_total'])) for e in es]
        fn = lambda x: np.interp(x, le, lW) - math.log(49.05678346694957)
        res['e_crit_bulk'] = math.exp(brentq(fn, le[0], le[-1])) if fn(le[0]) > 0 > fn(le[-1]) else None
        out['paths'][name] = res
        b6 = [x for x in res['by_e'] if x['e'] == 1e-6][0]
        print('%-10s s(N=1,2,3.8)=%s h_r=%s f(N=1,2,3.8)=%s | e=1e-6: f_eff=%.3f R=%.3f S/Q=%.3f -> %s ; e_crit(bulk)=%s' % (
            name, ['%.2f' % np.interp(x, Ng, sN) for x in (1, 2, 3.8)], ['%.2f' % np.interp(x, Ng, hN) for x in (1, 2, 3.8)],
            ['%.3f' % np.interp(x, Ng, fN) for x in (1, 2, 3.8)], b6['f_eff'], b6['R_total'], b6['S_per_Q'],
            'free X' if b6['S_per_Q'] > 3.4716191256166837 else 'Q-matter', '%.2e' % res['e_crit_bulk'] if res['e_crit_bulk'] else None))
    json.dump(out, open(RES / 'sk_exit_path.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
