#!/usr/bin/env python3
"""inf_critic (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): the quiet-inflaton / CMB / dark-matter
triangle for single-field alpha-attractors.

Question: can one single-field plateau inflaton be (a) quiet after inflation (m_I ~ 2-4 H_i at the slice, the only window
where the M1r release is M1r-like, inf_potential skeptic MODE-0 scan), (b) CMB-normalised (A_s) and (c) dark-matter
consistent with FAST reheating (K ~ 2, the only history that keeps M1r's e-controlled branch point)?

Method (all model choices):
  1. T- and E-model alpha-attractors on an alpha grid.  V0 fixed by A_s = 2.1e-9 at N_* = 52.5 (exact Mukhanov-Sasaki,
     inf_potential inf_common, READ-ONLY via inf_chain ch_common, READ-ONLY).  Slice S1 = 1 e-fold before eps_H = 1, as in
     inf_chain.  Output: H_i = H(S1), m_I/H_i, r, n_s.
  2. The dark-matter requirement H_i,req(Gamma_I) of a quiet release, from inf_chain's computed Gamma scans of TB1
     (alpha = 0.3) and TB2 (alpha = 1) (results/chain_TB*_g*.json, read-only).  The two curves agree to <= 6 %, so they are
     used as an approximately universal quiet-release curve [E] for the alpha range where m_I/H_i <= 4.
  3. For each quiet alpha: the width Gamma_cal that closes the calibration, K there, and S/Q(1e-6) by log interpolation.
     For fast reheating: the alpha whose H_i equals H_i,req(Gamma -> 1e-3 ... 0.5) and its m_I/H_i.
Labels: step 1 [C]; steps 2-3 [E] (interpolation / extrapolation of computed scans).
Output: ../results/quiet_scan.json
"""
import json, math, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]
sys.path.insert(0, str(TRK / 'inf_chain' / 'code'))
import ch_common as cc          # noqa: E402  read-only (imports inf_common, calibrate, dq_model_copy read-only)

ic = cc.ic
RES = HERE.parent / 'results'


def tmodel_row(kind, alpha, Ns=52.5):
    m = ic.Model(kind, alpha=alpha)
    bg = ic.background(m, 75.0)
    sp = ic.spectrum_at_pivot(m, bg, Ns)
    V0 = ic.AS_TARGET / sp['Pz']
    b = cc.Bench('%s%g' % (kind, alpha), m, V0, Nslice=1.0, gam=0.0)
    # H at the pivot
    f, fN = ic.bg_at(bg, Ns)
    hs, _ = ic.h_of(m, f, fN)
    Hstar = math.sqrt(V0) * hs * ic.MPL_GEV
    return dict(kind=kind, alpha=alpha, N_star=Ns, H_i_GeV=b.Hi, H_star_GeV=Hstar, m_I_over_Hi=b.mI, r=sp['r'],
                n_s=sp['ns'], alpha_s=sp['alpha_s'], eps_H_S1=b.eps_s)


def chain_curves():
    out = {}
    for key in ('TB1', 'TB2'):
        rows = []
        for f in sorted((TRK / 'inf_chain' / 'results').glob('chain_%s_g1.000e-*.json' % key)):
            d = json.load(open(f))
            rows.append(dict(gam=d['gam'], K=d['K'], SQ=d['e1e6']['S_per_Q'] if isinstance(d['e1e6'], dict) else None,
                             H_req=d['relic_e1e6']['H_i_req_GeV'], M=d['relic_e1e6']['M'], SQ_inf=d['S_per_Q_e_inf'],
                             e_crit=d['e_crit_bulk']))
        rows.sort(key=lambda r: r['gam'])
        out[key] = rows
    return out


def loginterp(x, xs, ys):
    return float(np.exp(np.interp(np.log(x), np.log(xs), np.log(ys))))


def main():
    grid = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0]
    rows = []
    for kind in ('T', 'E'):
        for a in grid:
            try:
                r = tmodel_row(kind, a)
            except Exception as exc:          # noqa: BLE001
                r = dict(kind=kind, alpha=a, error=repr(exc))
            rows.append(r)
            print(json.dumps(r), flush=True)
    cur = chain_curves()
    # universality check of the quiet-release requirement curve
    g1 = [r['gam'] for r in cur['TB1']]; h1 = [r['H_req'] for r in cur['TB1']]
    h2 = [r['H_req'] for r in cur['TB2']]
    univ = [h2[i] / h1[i] for i in range(len(h1))]
    # log-log slope at the fast end, extrapolation to Gamma = 0.1, 0.5 [E]
    s_fast = math.log(h1[-1] / h1[-2]) / math.log(g1[-1] / g1[-2])
    ext = {g: h1[-1] * (g / g1[-1]) ** s_fast for g in (0.01, 0.1, 0.5)}
    # quiet alphas: closing width, K, S/Q
    close = []
    gs = np.array(g1); Ks = np.array([r['K'] for r in cur['TB1']]); SQs = np.array([r['SQ'] for r in cur['TB1']])
    Hs = np.array(h1)
    for r in rows:
        if 'error' in r or r['m_I_over_Hi'] > 4.5:
            continue
        H = r['H_i_GeV']
        if H > Hs.max() or H < Hs.min():
            close.append(dict(kind=r['kind'], alpha=r['alpha'], H_i_GeV=H, note='outside scan'))
            continue
        # H_req is decreasing in gam: invert in log space
        lg = float(np.interp(math.log(H), np.log(Hs[::-1]), np.log(gs[::-1])))
        gc = math.exp(lg)
        close.append(dict(kind=r['kind'], alpha=r['alpha'], m_I_over_Hi=r['m_I_over_Hi'], H_i_GeV=H, r=r['r'], n_s=r['n_s'],
                          Gamma_cal_over_Hi=gc, K=loginterp(gc, gs, Ks), S_per_Q_1e6=loginterp(gc, gs, SQs)))
    # fast reheating: which alpha has H_i = H_req(Gamma) ?
    Tr = [r for r in rows if r['kind'] == 'T' and 'error' not in r]
    aT = np.array([r['alpha'] for r in Tr]); HT = np.array([r['H_i_GeV'] for r in Tr]); mT = np.array([r['m_I_over_Hi'] for r in Tr])
    rT = np.array([r['r'] for r in Tr])
    fast = []
    for g, H in list(zip(g1[-3:], h1[-3:])) + list(ext.items()):
        la = float(np.interp(math.log(H), np.log(HT), np.log(aT)))
        fast.append(dict(Gamma_over_Hi=g, H_req_GeV=H, alpha_T=math.exp(la),
                         m_I_over_Hi=float(np.exp(np.interp(la, np.log(aT), np.log(mT)))),
                         r=float(np.exp(np.interp(la, np.log(aT), np.log(rT)))),
                         extrapolated=g > g1[-1]))
    # quiet threshold alpha (m_I/H_i = 4 and 2) in the T-model, with r there
    thr = {}
    for target in (4.0, 3.0, 2.0):
        la = float(np.interp(math.log(target), np.log(mT[::-1]), np.log(aT[::-1])))
        thr['m_I=%g' % target] = dict(alpha=math.exp(la), r=float(np.exp(np.interp(la, np.log(aT), np.log(rT)))),
                                     H_i_GeV=float(np.exp(np.interp(la, np.log(aT), np.log(HT)))))
    out = dict(method=__doc__, rows=rows, quiet_curve_TB1=cur['TB1'], quiet_curve_TB2=cur['TB2'],
               universality_TB2_over_TB1=univ, fast_slope=s_fast, H_req_extrapolated=ext, quiet_closure=close,
               fast_closure=fast, quiet_threshold_T=thr)
    RES.mkdir(exist_ok=True)
    json.dump(out, open(RES / 'quiet_scan.json', 'w'), indent=1)
    print('universality', [round(u, 3) for u in univ])
    print('fast slope', s_fast, 'extrap', ext)
    print('thresholds', json.dumps(thr, indent=1))
    print('quiet closure', json.dumps(close, indent=1))
    print('fast closure', json.dumps(fast, indent=1))


if __name__ == '__main__':
    main()
