#!/usr/bin/env python3
"""dq_critic item 3 (Claude, own code, 2 Nov 2026): recalibrate H_i per branch with the frozen N_X and the corrected
mass per charge, and give r.

  H_i = H_M1 (0.9374 / (M Y/Y_Q))^{2/3}     [D: Y_Q ~ H_i^{3/2} at a fixed dimensionless history, exponent 1.5000000002]
  Y/Y_Q = 0.9910212 at the M1 H_i (frozen release, dq_entropy SKEPTIC N1; dq_entropy 0.991023).
  Free X: M = 1 + 2p, p(S/Q) from wt_critic's gas-only map (log-log in S/Q; [E]); at e = 1e-7 p ~ 0.7 from
          REVIEW_KINETIC (f = 1.49 at R = 41.7, i.e. W = 557 at the spot; here W = 589) [E].
  Q-matter: M/m_X from qm_micro SKEPTIC (energy-conserving f_v, Bose + heat <Q>), perfect ... Beliaev trapping,
          linear in S/Q between its points (2.5, 3.0, 3.2; 1.5 for extrapolation below 2.5).
  Portals vanishing at M_P: Y/Y_Q(L) = 1.305763 (1 + 0.0246 (L - 13.31)), L = ln(M_P/(30 H_i)) solved self-consistently
          (dq_entropy SKEPTIC N2).
  Systematic of the X <-> Sigma dephasing at fragmentation: N_X -2.2% ... +0.64% -> H_i +1.5% ... -0.4% [C on a toy].
  r = 2 H^2/(pi^2 A_s M_P^2), A_s = 2.0989e-9 (Planck 2018, ln(1e10 A_s) = 3.044) [R], M_P = 2.435e18 GeV, with
          H = H_i: the value if H at CMB exit equals H_i; a rolling history has H_CMB >= H_i, so it is a floor [A].
Output: ../results/hi_calibration.json
"""
import json, math, os
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
TR = os.path.abspath(os.path.join(HERE, '..', '..'))
HM1 = 1.342842735e11
YYQ = 0.9910212
MPR = 2.435e18
AS = math.exp(3.044) * 1e-10


def p_of_SQ(SQ):
    cc = json.load(open(os.path.join(TR, 'wt_critic/results/chain_consistency.json')))['p_scan_SQ']
    x = np.log([r['S_per_Q'] for r in cc]); y = np.log([r['p'] for r in cc])
    lx = math.log(SQ)
    if lx < x[0]:
        s = (y[1] - y[0]) / (x[1] - x[0]); return math.exp(y[0] + s * (lx - x[0]))
    if lx > x[-1]:
        # [E] log-log bridge from the map's last point to REVIEW_KINETIC's case-L anchor p = 1.49/2.02855 at
        # W = 13.0412 x 42.707 (S/Q = 21.80 by the closed form); extrapolated beyond it with the same slope
        xa, ya = math.log(21.80), math.log(1.49 / 2.02855)
        s = (ya - y[-1]) / (xa - x[-1])
        return math.exp(y[-1] + s * (lx - x[-1]))
    return math.exp(np.interp(lx, x, y))


def M_qmatter(SQ):
    S = [1.5, 2.5, 3.0, 3.2]
    lo = [0.7883, 0.7974, 0.8023, 0.8048]; hi = [0.7944, 0.8076, 0.8152, 0.8182]
    return float(np.interp(SQ, S, lo)), float(np.interp(SQ, S, hi))


def H_of(M, yyq=YYQ):
    return HM1 * (0.9374 / (M * yyq)) ** (2.0 / 3.0)


def r_of(H):
    return 2 * H * H / (math.pi ** 2 * AS * MPR ** 2)


def H_MP(M):
    f = lambda H: H - H_of(M, 1.305763 * (1 + 0.0246 * (math.log(MPR / (30 * H)) - 13.31)))
    H = brentq(f, 5e10, 3e11)
    return H, math.log(MPR / (30 * H))


def main():
    sb = json.load(open(os.path.join(HERE, '..', 'results', 'sq_branch.json')))
    out = dict(method=__doc__, constants=dict(H_M1=HM1, Y_over_YQ=YYQ, A_s=AS, M_P=MPR), rows=[])
    sig = 3.4716191256166837
    for r in sb['rows']:
        e, SQ = r['e'], r['SQ_final']
        row = dict(e=e, SQ=SQ)
        if SQ > sig:
            p = p_of_SQ(SQ)
            if SQ > 6.742:
                row['p_note'] = '[E] beyond the wt_critic map: bridged to REVIEW_KINETIC case L (p = 0.73 at S/Q = 21.8)'
            # the p map is uncertain by x2-2.5 at fixed S/Q (REVIEW_KINETIC case L gives p = 2.3e-3 at S/Q 4.94
            # against the map's 9.1e-4); quoted as p_alt
            row['p_alt_caseL'] = p * 2.5
            M = 1 + 2 * p
            Ma = 1 + 2 * min(row['p_alt_caseL'], 0.9)
            row.update(branch='free X', p=p, M=[M, Ma], H_i=[H_of(Ma), H_of(M)])
        else:
            Ml, Mh = M_qmatter(SQ)
            row.update(branch='Q-matter', M=[Ml, Mh], H_i=[H_of(Mh), H_of(Ml)])
        row['H_i_dephasing'] = [row['H_i'][0] * (1.0064) ** (-2 / 3), row['H_i'][1] * (0.978) ** (-2 / 3)]
        row['r'] = [r_of(h) for h in row['H_i']]
        # M_P portals: branch from the M_P S/Q
        SQp = r['SQ_MP']
        if SQp > sig:
            pp = p_of_SQ(SQp) or 0.5
            Hp, L = H_MP(1 + 2 * pp); row['MP'] = dict(SQ=SQp, branch='free X', H_i=[Hp, Hp], L=L)
        else:
            Ml, Mh = M_qmatter(SQp)
            H1, L1 = H_MP(Mh); H2, L2 = H_MP(Ml)
            row['MP'] = dict(SQ=SQp, branch='Q-matter', H_i=[H1, H2], L=[L1, L2])
        row['MP']['r'] = [r_of(h) for h in row['MP']['H_i']]
        out['rows'].append(row)
        print('e=%.1e S/Q=%.3f %-8s M=%.4f-%.4f H_i=%.4e-%.4e (deph %.4e-%.4e) r=%.2e-%.2e | M_P: S/Q=%.3f %-8s H_i=%.4e-%.4e r=%.2e' % (
            e, SQ, row['branch'], *row['M'], *row['H_i'], *row['H_i_dephasing'], *row['r'], SQp, row['MP']['branch'],
            *row['MP']['H_i'], row['MP']['r'][0]))
    # rolling-inflaton embeddings at the declared e = 1e-6: Q-matter
    r6 = [r for r in sb['rows'] if r['e'] == 1e-6][0]
    out['e1e6_rolling_paths'] = {k: dict(SQ=v, M=M_qmatter(v), H_i=[H_of(M_qmatter(v)[1]), H_of(M_qmatter(v)[0])])
                                 for k, v in r6['SQ_paths'].items()}
    print('e = 1e-6 on rolling paths:', out['e1e6_rolling_paths'])
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'hi_calibration.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
