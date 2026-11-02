#!/usr/bin/env python3
"""dq_entropy item 4 + summary (Claude, own code).

S/Q as a function of R for the M1r release (item 2 scan + item 3 increments), and S/Q(e) for two readings of the heat:
  (A) 'provisional, literal': R(e) = 0.8056 + 4.0538 x + 0.10716 x^0.5, x = 1e-6/e (fx_chain skeptic; spot-release units)
      plugged in as R of the M1r release, as the brief asks.  This keeps the spot's Y_L debris (0.8056) and its vector heat
      per E_X, so it is an upper reading for M1r.
  (B) 'transferred': M1r's own mediator ledger (item 1: R_Y + (R_casc + R_Sigma) x^0.5) plus the spot vector heat moved in
      E_c units with wt_critic's nu = (H_i/H_spot)^2 (Y_spot/Y)^{4/3} [E]:  R = R_med(e) + nu (E/E_c)_spot/(E/E_c)_M1r R_vec(e),
      band (4.0538 x) and broad (4.7905 x); nu and H_i iterated with the branch's mass per charge.
Branch: free X if S/Q > 3.4716 (bulk Q-matter cannot survive), Q-matter if below; Q = 2e5 objects need S/Q < 2.9284.
p(S/Q) from wt_critic's mapping of fx_chain's gas-only network (results/chain_consistency.json, read-only) [E].
Usage: summary.py ../results  (reads freeze.json, sq_scan.json, entropy_increments.json; writes summary.json)
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

RES = Path(sys.argv[1])
WT = Path(__file__).resolve().parents[2]
FR = json.load(open(RES / 'freeze.json'))
SQ = json.load(open(RES / 'sq_scan.json'))
EI = json.load(open(RES / 'entropy_increments.json'))
CRIT = json.load(open(WT / 'wt_critic' / 'results' / 'chain_consistency.json'))

H_M1, YQ_REF, D, MX = 1.342842735e11, 5.07941423e-9, 4.4e-10, 0.0924089637375362
H_SPOT, Y_SPOT, W0_SPOT = 1.3334961939833566e11, 4.8294009076988175e-09, CRIT['spot']['E_over_Ec']
SIG = SQ['sigma_min']
RVEC = dict(band=4.0538, broad=4.7905)
pq = CRIT['p_scan_SQ']
LSQ, LP = np.log([x['S_per_Q'] for x in pq]), np.log([x['p'] for x in pq])


def p_of_SQ(s):
    """log-log interpolation; linear extrapolation in log below S/Q = 3.64 (p < 1e-4, irrelevant for M)."""
    x = math.log(s)
    if x <= LSQ[0]:
        sl = (LP[1] - LP[0]) / (LSQ[1] - LSQ[0]); return float(math.exp(LP[0] + sl * (x - LSQ[0])))
    if x >= LSQ[-1]:
        sl = (LP[-1] - LP[-2]) / (LSQ[-1] - LSQ[-2]); return float(math.exp(LP[-1] + sl * (x - LSQ[-1])))
    return float(math.exp(np.interp(x, LSQ, LP)))


class SofR:
    def __init__(self, key):
        r = SQ['releases'][key]
        self.w0 = r['E_over_Ec']
        self.R = np.array([x['R'] for x in r['scan']]); self.S = np.array([x['S_per_Q'] for x in r['scan']])
        self.T = np.array([x['T_h_over_mX_on'] for x in r['scan']]); self.n = np.array([x['nQ_on_mX3'] for x in r['scan']])
        self.Xb = np.array([x['Xbar_per_Q'] for x in r['scan']])
        pts = EI['points']
        Rd = np.array(sorted(p['R'] for p in pts.values()))
        order = sorted(pts.values(), key=lambda p: p['R'])
        self.dR = Rd; self.dT0 = np.array([p['T0']['dSQ'] for p in order]); self.dTh = np.array([p['thermal']['dSQ'] for p in order])
        self.dsurf = max(max(p['surface']['dSQ_if_released_at_T_0p8']) for p in order)
        # cubic spline-free: monotone interpolation in R on the 0.25 grid (S/Q is smooth; check vs direct solve below)

    def hartree(self, R):
        if R > self.R[-1] or R < 0:
            raise ValueError('R outside [0, 8]: %g' % R)
        return float(np.interp(R, self.R, self.S)) if False else float(self._interp(R))

    def _interp(self, R):
        # 4-point Lagrange on the uniform grid (S/Q is smooth in R)
        h = self.R[1] - self.R[0]; i = int(min(max(math.floor(R / h) - 1, 0), len(self.R) - 4))
        xs, ys = self.R[i:i + 4], self.S[i:i + 4]
        return sum(ys[j] * np.prod([(R - xs[k]) / (xs[j] - xs[k]) for k in range(4) if k != j]) for j in range(4))

    def final(self, R, mode='T0'):
        """S/Q after the item-3 increments: decomposition latent heat (T0 [C] or thermal [E]) + surface upper bound."""
        d = np.interp(R, self.dR, self.dT0 if mode == 'T0' else self.dTh)
        return self.hartree(R) + float(d) + (self.dsurf if mode == 'T0' else 0.0)


def H_for(yield_ratio, M):
    return H_M1 * (0.9374 / (M * yield_ratio)) ** (2.0 / 3.0)


def nu_of(H, M):
    Y = D / (M * MX)
    return (H / H_SPOT) ** 2 * (Y_SPOT / Y) ** (4.0 / 3.0)


def R_prov(e):
    x = 1e-6 / e
    return 0.8056 + 4.0538 * x + 0.10716 * x ** 0.5


def main():
    out = dict(sigma_min=SIG)
    rel = dict(benchmark=('benchmark_cMY0', 'base_1e5'), planck_bc=('planck_bc', 'planck_1e5'))
    es = [1e-6, 2e-6, 3e-6, 5e-6]
    for name, (sk, fk) in rel.items():
        f = FR[fk]; S = SofR(sk)
        y = f['yield_over_YQ_frozen']
        led = f['ledger']['1e-06']
        RY, RCS = led['R_Y'], led['R_cascade'] + led['R_Sigma']
        Rmed = lambda e: RY + RCS * (1e-6 / e) ** 0.5
        o = dict(E_over_Ec=S.w0, N_X=f['frozen']['NX']['value'], E_X=f['E_X_late'], yield_over_YQ=y, R_Y=RY, R_casc_plus_Sigma_e1e6=RCS,
                 R_crit={k: v['R_crit'] for k, v in SQ['releases'][sk]['R_crit'].items()})
        # R_crit including the item-3 increments
        o['R_crit_final'] = {}
        for k, s in SIG.items():
            row = {}
            for mode in ('T0', 'thermal'):
                try:
                    row[mode] = brentq(lambda R: S.final(R, mode) - s, 0.0, 8.0, xtol=1e-9)
                except ValueError:
                    row[mode] = None
            o['R_crit_final'][k] = row
        # (A) literal provisional R(e)
        A = {}
        for e in es:
            R = R_prov(e); s_h = S.hartree(R); s_f = S.final(R)
            A['%g' % e] = dict(R=R, S_per_Q_hartree=s_h, S_per_Q_final=s_f, branch='free X' if s_f > SIG['bulk_Q_inf'] else 'Q-matter',
                              Q2e5_objects_survive=bool(s_f < SIG['Q_2e5']))
        eA = brentq(lambda e: S.final(R_prov(e)) - SIG['bulk_Q_inf'], 1e-6, 2e-5)
        eA2 = brentq(lambda e: S.final(R_prov(e)) - SIG['Q_2e5'], 1e-6, 1e-4)
        o['reading_A_literal'] = dict(points=A, e_crit_bulk=eA, e_crit_Q2e5=eA2)
        # (B) transferred, nu and H_i iterated with the branch
        B = {}
        for src in ('band', 'broad'):
            def state(e):
                M = 1.0011; H = H_for(y, M)
                for _ in range(8):
                    nu = nu_of(H, M)
                    R = Rmed(e) + nu * W0_SPOT / S.w0 * RVEC[src] * (1e-6 / e)
                    s_f = S.final(R)
                    if s_f > SIG['bulk_Q_inf']:
                        M = 1 + 2 * p_of_SQ(s_f)
                    else:
                        M = 0.79                  # Q-matter relic, eta ~ 0.89-0.97 (fx_formation [E]): M = 0.782-0.802
                    H = H_for(y, M)
                return dict(R=R, nu=nu, H_i_GeV=H, M_per_charge=M, S_per_Q_hartree=S.hartree(R), S_per_Q_final=s_f,
                            branch='free X' if s_f > SIG['bulk_Q_inf'] else 'Q-matter', Q2e5_objects_survive=bool(s_f < SIG['Q_2e5']),
                            p_gas=p_of_SQ(s_f) if s_f > SIG['bulk_Q_inf'] else None)
            pts = {('%g' % e): state(e) for e in es}
            eb = brentq(lambda e: state(e)['S_per_Q_final'] - SIG['bulk_Q_inf'], 8e-7, 2e-5)
            eq = brentq(lambda e: state(e)['S_per_Q_final'] - SIG['Q_2e5'], 8e-7, 1e-4)
            B[src] = dict(points=pts, e_crit_bulk=eb, e_crit_Q2e5=eq)
        o['reading_B_transferred'] = B
        # compact S/Q(R) table for the critic (Hartree; + increments)
        o['SQ_of_R_table'] = [dict(R=float(R), S_per_Q_hartree=float(s), S_per_Q_final_T0=S.final(float(R)), S_per_Q_final_thermal=S.final(float(R), 'thermal'),
                                   T_h_over_mX_onset=float(T), nQ_onset_mX3=float(n), Xbar_per_Q=float(xb))
                              for R, s, T, n, xb in zip(S.R, S.S, S.T, S.n, S.Xb)]
        # cubic fit of the Hartree S/Q(R) on [0, 8]
        c = np.polyfit(S.R, S.S, 4)
        o['SQ_fit_quartic_in_R'] = dict(coef_highest_first=c.tolist(), max_abs_resid=float(np.max(np.abs(np.polyval(c, S.R) - S.S))))
        out[name] = o
        print('== %s: E/E_c = %.4f, Y/Y_Q = %.6f, R_Y = %.4f, R_casc+Sig(1e-6) = %.4f' % (name, S.w0, y, RY, RCS))
        print('   R_crit (Hartree):', {k: round(v, 4) for k, v in o['R_crit'].items() if v})
        print('   R_crit (final T0 / thermal):', {k: (round(v['T0'], 4), round(v['thermal'], 4)) for k, v in o['R_crit_final'].items()})
        for e in es:
            a = A['%g' % e]; b = B['band']['points']['%g' % e]; bb = B['broad']['points']['%g' % e]
            print('   e=%g  (A) R=%.3f S/Q=%.3f %-8s | (B band) R=%.3f S/Q=%.3f %-8s H=%.3e nu=%.3f | (B broad) R=%.3f S/Q=%.3f %s' % (
                e, a['R'], a['S_per_Q_final'], a['branch'], b['R'], b['S_per_Q_final'], b['branch'], b['H_i_GeV'], b['nu'], bb['R'], bb['S_per_Q_final'], bb['branch']))
        print('   e_crit bulk: (A) %.3e  (B band) %.3e  (B broad) %.3e ; Q=2e5: (A) %.3e (B band) %.3e (B broad) %.3e' % (
            eA, B['band']['e_crit_bulk'], B['broad']['e_crit_bulk'], eA2, B['band']['e_crit_Q2e5'], B['broad']['e_crit_Q2e5']))
    # H_i ledger at the benchmark
    yb = FR['base_1e5']['yield_over_YQ_frozen']; yp = FR['planck_1e5']['yield_over_YQ_frozen']
    pb = out['benchmark']['reading_B_transferred']['band']['points']['1e-06']
    out['H_ledger'] = dict(freeX_benchmark=pb['H_i_GeV'], freeX_benchmark_M=pb['M_per_charge'],
                           Q2e5_relic_benchmark=H_for(yb, 0.9374), Qmatter_benchmark=[H_for(yb, 0.802), H_for(yb, 0.782)],
                           freeX_planck=out['planck_bc']['reading_B_transferred']['band']['points']['1e-06']['H_i_GeV'],
                           Qmatter_planck=[H_for(yp, 0.802), H_for(yp, 0.782)],
                           yield_shift_vs_400_1000=FR['base_1e5']['frozen_over_400_1000'])
    print('H ledger:', json.dumps(out['H_ledger'], indent=0))
    # induced-portal sensitivity per unit ln Q at the benchmark (cport = +-1/(4 pi^2) runs, item 1), e = 1e-6, reading B:
    # W = (E/E_c)(1 + R_med) + nu W_vec,spot with nu ∝ y^{-4/3}; S/Q from the benchmark S(R) at the same W (exact invariance)
    Sb = SofR('benchmark_cMY0'); nu_b = nu_of(H_for(yb, 1.0011), 1.0011)
    der = {}
    for tag in ('cp_plus1_2e4', 'cp_minus1_2e4', 'base_1e5', 'planck_1e5'):
        if tag not in FR:
            continue
        f = FR[tag]; w0 = f['E_over_Ec']; led = f['ledger']['1e-06']
        nu = nu_b * (yb / f['yield_over_YQ_frozen']) ** (4.0 / 3.0)
        row = {}
        for src in ('band', 'broad'):
            W = w0 * (1 + led['R_central']) + nu * W0_SPOT * RVEC[src]
            Reff = W / Sb.w0 - 1
            row[src] = dict(W=W, S_per_Q_final=Sb.final(Reff), nu=nu)
        der[tag] = dict(E_over_Ec=w0, R_med=led['R_central'], yield_over_YQ=f['yield_over_YQ_frozen'], **row)
    if 'cp_plus1_2e4' in der:
        der['dSQ_per_unit_lnQ_band'] = (der['cp_plus1_2e4']['band']['S_per_Q_final'] - der['cp_minus1_2e4']['band']['S_per_Q_final']) / 2
        der['dlnH_per_unit_lnQ'] = -(2.0 / 3.0) * (math.log(der['cp_plus1_2e4']['yield_over_YQ']) - math.log(der['cp_minus1_2e4']['yield_over_YQ'])) / 2
    out['portal_sensitivity_e1e6'] = der
    # release-independent thresholds: S/Q depends only on W = (E/E_c)(1 + R) (item 2 invariance), so for any release
    # R_crit = W_crit/(E/E_c) - 1
    rb = out['benchmark']
    out['W_crit'] = {k: dict(hartree=rb['E_over_Ec'] * (1 + v) if v is not None else None,
                             final_T0=rb['E_over_Ec'] * (1 + rb['R_crit_final'][k]['T0']),
                             final_thermal=rb['E_over_Ec'] * (1 + rb['R_crit_final'][k]['thermal']))
                     for k, v in rb['R_crit'].items()}
    out['for_critic'] = dict(
        definition='R = injected conformal heat / E_X of the frozen M1r release (N_X = 8047.84, E_X = 188737.7, a = 1 at the slice, H_i units); '
                   'X-sector energy before onset (1 + R) E_X at fixed N_X',
        SQ_of_R='benchmark.SQ_of_R_table (0.25 grid; Hartree and final = Hartree + item-3 increments) or benchmark.SQ_fit_quartic_in_R',
        R_crit_benchmark=rb['R_crit'], R_crit_final_benchmark=rb['R_crit_final'], W_crit=out['W_crit'],
        rule='branch = free X if S/Q(R(e)) > 3.4716, Q-matter if below; Q = 2e5 objects need S/Q < 2.9284')
    # release-independent closed form S/Q(W), fitted on the benchmark scan and tested on the Planck-BC scan (different E/E_c)
    from scipy.optimize import curve_fit
    fW = lambda W, A, B, p, C: A * np.maximum(W - B, 1e-9) ** p + C
    Rg = np.array([x['R'] for x in rb['SQ_of_R_table']]); Sg = np.array([x['S_per_Q_hartree'] for x in rb['SQ_of_R_table']])
    cW, _ = curve_fit(fW, rb['E_over_Ec'] * (1 + Rg), Sg, p0=[0.3, 10, 0.75, 0.0], maxfev=20000)
    rp = SQ['releases']['planck_bc']
    Wp = np.array([rp['E_over_Ec'] * (1 + x['R']) for x in rp['scan']]); Sp = np.array([x['S_per_Q'] for x in rp['scan']])
    out['SQ_of_W_closed_form'] = dict(formula='S/Q = A (W - B)^p + C, W = (E/E_c)(1 + R) (Hartree, before item-3 increments)',
                                      A=cW[0], B=cW[1], p=cW[2], C=cW[3], W_range_fit=[float(rb['E_over_Ec']), float(rb['E_over_Ec'] * 9)],
                                      max_resid_fit=float(np.max(np.abs(fW(rb['E_over_Ec'] * (1 + Rg), *cW) - Sg))),
                                      max_resid_planck_scan=float(np.max(np.abs(fW(Wp, *cW) - Sp))), W_range_planck=[float(Wp[0]), float(Wp[-1])])
    out['for_critic']['SQ_of_W'] = out['SQ_of_W_closed_form']
    print('S/Q(W) closed form:', out['SQ_of_W_closed_form'])
    print('W_crit:', json.dumps(out['W_crit'], indent=0))
    print('quartic fit max resid:', rb['SQ_fit_quartic_in_R']['max_abs_resid'], rb['SQ_fit_quartic_in_R']['coef_highest_first'])
    print('portal sensitivity (e = 1e-6, B):', json.dumps({k: (v if not isinstance(v, dict) else {kk: (vv['S_per_Q_final'] if isinstance(vv, dict) else vv) for kk, vv in v.items()}) for k, v in der.items()}, indent=0))
    json.dump(out, open(RES / 'summary.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
