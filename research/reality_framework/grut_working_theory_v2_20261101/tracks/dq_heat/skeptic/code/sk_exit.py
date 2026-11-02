#!/usr/bin/env python3
"""dq_heat SKEPTIC (Claude, own code): inflaton-exit sensitivity of the vector heat with H allowed to differ at exit.

dq_heat/finf_family.py moves only the holder (s = phi_I/phi_I(0)) at fixed H and then *estimates* [E] that a Hubble rate
growing with phi_I (large field, H proportional to phi_I) would nearly cancel the effect: s^2 (F(1)/F(s))^2 = 1.06-1.15.
That estimate keeps F(s) from the fixed-H family, i.e. it ignores that the prepared X condensate is held by the curvature
mass xi R |X|^2 = 12 xi H^2 |X|^2, which grows with H.  Here the conditional minimum is re-solved at the exit Hubble rate
(H_exit = h_r H_i, R = 12 h_r^2 in H_i units; physical M_Y, h, mu, kappa, lambda's unchanged; u = mu + s (10i - mu)), and
the frozen Goldstone angle gives the inflationary (piece-A) vector number scaling
    f(s, h_r) = [ (H_exit/F_exit) / (H_i/F_inf) ]^2 = [ h_r F(1,1) / F(s, h_r) ]^2.
R_vec,A(e) -> f R_vec,A(e); the release-era band (k >= 1) and the mediators are unchanged.  S/Q from dq_entropy's closed form
S/Q = 0.190255 (W - 0.97494)^0.75013 - 0.00428, W = |eps|^{-4/3} (1 + R) [Q dq_entropy, verified against dq_heat's column].
copied model: dq_model_copy.py (unchanged copy of dq_heat/code/dq_model.py).   Labels: [C] on the static family.
Output: ../results/sk_exit.json
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


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def main():
    p = fc.params(0.5)
    hb = json.load(open(DQ / 'heat_table.json'))['budgets']['portals_vanish_at_MY_benchmark']
    W0 = hb['W0']
    rows_e = {r['e']: r for r in hb['rows']}
    r6 = rows_e[1e-6]
    other = r6['R_total'] - r6['R_vec_k_le_1']          # band + SU + mediators + cascade (unchanged by the exit state)
    # f needed for S/Q = 3.4716 at e = 1e-6
    Wthr = 49.05678346694957
    f_crit = ((Wthr / W0 - 1) - other) / r6['R_vec_k_le_1']
    ref = fc.stationary(fc.U_PREP, p, R=12.0)
    F11 = float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(ref['z']) ** 2)))
    out = dict(method=__doc__, F_inf_ref=F11, f_crit_e1e6=f_crit, rows=[])
    print('F(1,1) = %.4f  f_crit(e=1e-6) = %.4f  (H_exit/F_exit must stay >= %.4f of the slice value)' % (F11, f_crit, math.sqrt(f_crit)))
    grid = []
    for s in (1.0, 1.25, 1.5, 2.0, 2.5):
        for hr in (1.0, 1.25, 1.5, 2.0, 2.5):
            grid.append((s, hr))
    # embeddings: H proportional to phi_I (large field) continued, and a plateau-like path h_r = 1 + 0.4 (s - 1)
    grid += [(2.25, 2.25), (2.4, 2.4), (3.0, 3.0), (4.0, 4.0), (5.0, 5.0)]
    grid += [(2.0, 1.4), (2.6, 1.64), (3.0, 1.8)]
    guess = {}
    for s, hr in grid:
        u = p['mu'] + (fc.U_PREP - p['mu']) * s
        st = None
        for g in (guess.get(s), guess.get(round(s - 0.25, 2)), None):
            try:
                st = fc.stationary(u, p, R=12.0 * hr * hr, guess=g)
                if st['grad'] < 1e-8 and st['eig'][1] > -1e-6:
                    break
            except Exception:
                st = None
        guess[s] = st['z']
        F = float(np.sqrt(2 * np.sum(Q ** 2 * np.abs(st['z']) ** 2)))
        f = (hr * F11 / F) ** 2
        Rt = other + f * r6['R_vec_k_le_1']
        sq = SQ(W0 * (1 + Rt))
        # bulk e_crit with the same scaling at every e (R_vec,A ~ e^-1.03 interpolated from the table)
        es = sorted(rows_e); le = np.log(es)
        lW = [math.log(W0 * (1 + rows_e[e]['R_total'] - rows_e[e]['R_vec_k_le_1'] + f * rows_e[e]['R_vec_k_le_1'])) for e in es]
        fn = lambda x: np.interp(x, le, lW) - math.log(Wthr)
        ec = math.exp(brentq(fn, le[0], le[-1])) if fn(le[0]) > 0 > fn(le[-1]) else None
        row = dict(s=s, h_r=hr, amps=np.abs(st['z']).tolist(), F_exit_over_Hi=F, F_exit_over_Hexit=F / hr,
                   lightest_m2_over_Hexit2=float(sorted(st['eig'])[1] / hr ** 2), grad=st['grad'], f=f, R_total_e1e6=Rt, SQ_e1e6=sq,
                   branch_e1e6='free X' if sq > 3.4716191256166837 else 'Q-matter', e_crit_bulk=ec)
        out['rows'].append(row)
        print('s=%.2f h_r=%.2f |X|,|Y|,|S| = %7.2f %6.2f %6.2f  F/H_exit=%7.2f  m2min/H^2=%6.2f  f=%.3f  R(1e-6)=%.3f  S/Q=%.3f  %-8s e_crit=%s' % (
            s, hr, *np.abs(st['z']), F / hr, row['lightest_m2_over_Hexit2'], f, Rt, sq, row['branch_e1e6'],
            '%.2e' % ec if ec else 'outside 1e-7..1e-5'))
    # dq_heat's [E] estimate for H proportional to phi_I: s^2 (F(1)/F(s, fixed H))^2
    fixH = {r['s']: r['F_exit_over_Hi'] for r in out['rows'] if r['h_r'] == 1.0}
    out['dq_heat_estimate_vs_resolved'] = [dict(s=s, estimate=s * s * (F11 / fixH[s]) ** 2,
                                                resolved=[r['f'] for r in out['rows'] if r['s'] == s and abs(r['h_r'] - s) < 1e-9][0])
                                           for s in (1.25, 1.5, 2.0, 2.5)]
    for d in out['dq_heat_estimate_vs_resolved']:
        print('H prop phi_I, s = %.2f: dq_heat estimate s^2 (F1/Fs)^2 = %.3f ; re-solved at H_exit = s H_i: f = %.3f' % (d['s'], d['estimate'], d['resolved']))
    json.dump(out, open(RES / 'sk_exit.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
