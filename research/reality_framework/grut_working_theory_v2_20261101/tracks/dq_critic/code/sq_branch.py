#!/usr/bin/env python3
"""dq_critic item 1 (Claude, own code, 2 Nov 2026): combine dq_heat's computed heat R(e) with dq_entropy's S/Q(R).

Inputs (all read-only, from reviewed tracks of this pass):
  dq_heat/results/heat_table.json         R(e) components on the M1r release (both portal conditions)   [C dq_heat]
  dq_heat/skeptic/results/sk_exit_path.json  f_eff(e) for three rolling-inflaton embeddings             [C on E paths]
  dq_entropy/skeptic/results/sk_late.json frozen release N_X, E_X (a <= 2582)                          [C]
  dq_entropy/results/summary.json         W_crit (Hartree and with item-3 entropy), S/Q(W) closed form [C]
Own steps:
  * frame: dq_heat's R is per its 400-1000 window E_X (188706.29); the relaxed plasma is set by the frozen release:
        W = (E_X,frozen + R E_X,win) / E_c(N_X,frozen),  E_c(N) = (N / c)^{4/3},  c = 2 sqrt(2 lX)(3 lX)^{-3/4};
  * S/Q(W) from the closed form, checked against a direct solve with an unchanged copy of dq_entropy's Hartree solver
    (hartree_copy.py) at three e;
  * entropy after relaxation: + 0.004 ... + 0.007 (dq_entropy item 3, confirmed by its skeptic);
  * error band on R: (i) computed multi-field check of the release-era band, |dR| <= 0.003 (1e-6/e)^{1/2}
    [dq_heat SKEPTIC]; (ii) dq_heat's outer full-vs-noX band +-0.05 (1e-6/e)^{1/2} [A]; (iii) omitted channels
    -1.7e-3 ... +2.6e-4 (dq_heat SKEPTIC N4) [E];
  * e_crit by root finding in ln e on log-log interpolated components (piece A ~ e^-1, band/cascade ~ e^-1/2);
  * variants: induced portals vanishing at M_P (dq_heat budget, with dq_entropy-skeptic's self-consistent-L shift
    dS/Q ~ -0.01 [E]); the rolling-inflaton embeddings (piece A x f_eff(e)).
Output: ../results/sq_branch.json.  Labels: [C] on the stated inputs.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
TR = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)

LX = 1e-3
C_NC = 2 * math.sqrt(2 * LX) * (3 * LX) ** (-0.75)        # 6.977568, N_c = c E^{3/4}
SIG = dict(bulk=3.4716191256166837, Q3e7=3.3839721242149268, Q2e5=2.9283778774213727, bind=2.5530624453282833)
INC = (0.0041, 0.0068)                                    # item-3 entropy after relaxation (thermal ... T=0 + surface)


def load():
    ht = json.load(open(os.path.join(TR, 'dq_heat/results/heat_table.json')))['budgets']
    ex = json.load(open(os.path.join(TR, 'dq_heat/skeptic/results/sk_exit_path.json')))['paths']
    late = json.load(open(os.path.join(TR, 'dq_entropy/skeptic/results/sk_late.json')))
    summ = json.load(open(os.path.join(TR, 'dq_entropy/results/summary.json')))
    return ht, ex, late, summ


def SQ_closed(W, cf):
    return cf['A'] * (W - cf['B']) ** cf['p'] + cf['C']


class Budget:
    """log-log interpolation of the R(e) components of one dq_heat budget."""
    def __init__(self, b, Nf=None, EXf=None):
        rows = sorted(b['rows'], key=lambda r: r['e'])
        self.le = np.log([r['e'] for r in rows])
        self.A = np.log([r['R_vec_k_le_1'] for r in rows])                       # piece A (frozen angle)
        self.B = np.log([r['R_vec_k_ge_1_release_era'] + r['R_vec_pieceB_SU'] for r in rows])
        self.C = np.log([r['R_cascade'] for r in rows])
        self.Y = rows[0]['R_Y_hom'] + rows[0]['R_Y_k']
        self.EXw = b['EX']; self.Nw = b['NX']
        self.Nf = Nf or b['NX']; self.EXf = EXf or b['EX']
        self.Ec = (self.Nf / C_NC) ** (4.0 / 3.0)
        self.rows = rows

    def _i(self, arr, e):
        x = math.log(e)
        if x < self.le[0]:
            s = (arr[1] - arr[0]) / (self.le[1] - self.le[0]); return math.exp(arr[0] + s * (x - self.le[0]))
        if x > self.le[-1]:
            s = (arr[-1] - arr[-2]) / (self.le[-1] - self.le[-2]); return math.exp(arr[-1] + s * (x - self.le[-1]))
        return math.exp(np.interp(x, self.le, arr))

    def parts(self, e):
        return dict(A=self._i(self.A, e), B=self._i(self.B, e), C=self._i(self.C, e), Y=self.Y)

    def R(self, e, fA=1.0, dR=0.0):
        p = self.parts(e)
        return fA * p['A'] + p['B'] + p['C'] + p['Y'] + dR

    def W(self, R):
        return (self.EXf + R * self.EXw) / self.Ec


def main():
    ht, ex, late, summ = load()
    cf = summ['for_critic']['SQ_of_W']
    fits = late['long_run']['fits']
    Nf = abs([f for f in fits['NX'] if f['model'] == 'a-2' and f['a_min'] > 800][0]['N_inf'])
    EXf = 188737.6639                                      # dq_entropy SKEPTIC N1 (exact conformal energy)
    bY = Budget(ht['portals_vanish_at_MY_benchmark'], Nf, EXf)
    bP = Budget(ht['portals_vanish_at_MP'])                 # M_P release: window values (frame shift negligible)
    out = dict(method=__doc__, frozen=dict(N_X=Nf, E_X=EXf, E_c=bY.Ec, W0=EXf / bY.Ec), rows=[], e_crit={}, checks={})

    def sq(b, e, fA=1.0, dR=0.0, inc=0.0, shift=0.0):
        return SQ_closed(b.W(b.R(e, fA, dR)), cf) + inc + shift

    def band(b, e, fA=1.0, shift=0.0, outer=False):
        x = (1e-6 / e) ** 0.5
        dlo = -(0.05 if outer else 0.003) * x - 1.7e-3
        dhi = +(0.05 if outer else 0.003) * x + 1.8e-4 / x + 8e-5
        return (sq(b, e, fA, dlo, INC[0], shift), sq(b, e, fA, 0.0, 0.5 * sum(INC), shift), sq(b, e, fA, dhi, INC[1], shift))

    paths = {k: v for k, v in ex.items() if k != 'static'}
    fpath = {k: (np.log([r['e'] for r in v['by_e']]), np.array([r['f_eff'] for r in v['by_e']])) for k, v in paths.items()}

    def fA_of(path, e):
        if path is None: return 1.0
        le, f = fpath[path]
        return float(np.interp(math.log(e), le, f))

    egrid = [1e-7, 3e-7, 1e-6, 1.5e-6, 2e-6, 2.5e-6, 3e-6, 4e-6, 5e-6, 1e-5, 3e-5, 1e-4]
    for e in egrid:
        p = bY.parts(e); R = bY.R(e)
        lo, mid, hi = band(bY, e); lo_o, _, hi_o = band(bY, e, outer=True)
        row = dict(e=e, R_parts=p, R_total=R, W=bY.W(R), SQ_hartree=SQ_closed(bY.W(R), cf),
                   SQ_final=mid, SQ_band_computed=[lo, hi], SQ_band_outer=[lo_o, hi_o],
                   branch='free X' if lo > SIG['bulk'] else ('Q-matter (bulk)' if hi < SIG['bulk'] else 'marginal'),
                   Q2e5_survive=hi < SIG['Q2e5'],
                   SQ_MP=sq(bP, e, inc=0.5 * sum(INC), shift=-0.01),
                   SQ_paths={k: sq(bY, e, fA_of(k, e), inc=0.5 * sum(INC)) for k in paths})
        out['rows'].append(row)
        print('e=%.2e R=%.4f W=%.3f S/Q=%.4f [%.4f, %.4f] outer [%.4f, %.4f] %-16s M_P %.3f  paths %s' % (
            e, R, row['W'], mid, lo, hi, lo_o, hi_o, row['branch'], row['SQ_MP'],
            ' '.join('%s %.3f' % (k, v) for k, v in row['SQ_paths'].items())))

    def ecrit(fun, target):
        g = lambda x: fun(math.exp(x)) - target
        return math.exp(brentq(g, math.log(1e-8), math.log(1e-3), xtol=1e-10))

    for th in ('bulk', 'Q2e5', 'bind'):
        t = SIG[th]
        d = dict(
            central=ecrit(lambda e: sq(bY, e, inc=0.5 * sum(INC)), t),
            hartree_only=ecrit(lambda e: sq(bY, e), t),
            computed_band=[ecrit(lambda e: band(bY, e)[0], t), ecrit(lambda e: band(bY, e)[2], t)],
            outer_band=[ecrit(lambda e: band(bY, e, outer=True)[0], t), ecrit(lambda e: band(bY, e, outer=True)[2], t)],
            portals_MP=ecrit(lambda e: sq(bP, e, inc=0.5 * sum(INC), shift=-0.01), t),
            paths={k: ecrit(lambda e, k=k: sq(bY, e, fA_of(k, e), inc=0.5 * sum(INC)), t) for k in paths})
        out['e_crit'][th] = d
        print('e_crit(%s, S/Q=%.4f): central %.4e  hartree %.4e  band [%.4e, %.4e]  outer [%.4e, %.4e]  M_P %.4e  paths %s' % (
            th, t, d['central'], d['hartree_only'], *d['computed_band'], *d['outer_band'], d['portals_MP'],
            ' '.join('%s %.3e' % (k, v) for k, v in d['paths'].items())))

    # f_eff needed for free X at the declared e = 1e-6 (static release, central increments)
    out['f_crit_e1e6'] = brentq(lambda f: sq(bY, 1e-6, f, inc=0.5 * sum(INC)) - SIG['bulk'], 0.05, 1.5)
    # margins at e = 2e-6 (the contradiction between dq_entropy reading B and the direct heat)
    r2 = [r for r in out['rows'] if r['e'] == 2e-6][0]
    out['checks']['e2e6'] = dict(SQ=r2['SQ_final'], band=r2['SQ_band_computed'], outer=r2['SQ_band_outer'],
                                 margin_central=r2['SQ_final'] / SIG['bulk'] - 1,
                                 reading_B_band_dq_entropy=3.370, note='direct heat supersedes the transfer')
    # direct Hartree solve with the copied solver at three e
    import hartree_copy as HC
    HC.set_release(Nf, EXf)
    for e in (1e-6, 2e-6, 3e-6):
        R = bY.R(e); Reff = R * bY.EXw / EXf
        st, _ = HC.solve(Reff)
        out['checks']['direct_hartree_e%.0e' % e] = dict(R=R, R_frozen_units=Reff, SQ_direct=st['S'] / Nf,
                                                        SQ_closed=SQ_closed(bY.W(R), cf), resid=st['resid'])
        print('direct Hartree e=%.0e: S/Q = %.5f (closed form %.5f)' % (e, st['S'] / Nf, SQ_closed(bY.W(R), cf)))
    print('f_eff needed for free X at e = 1e-6: %.4f' % out['f_crit_e1e6'])
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'sq_branch.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
