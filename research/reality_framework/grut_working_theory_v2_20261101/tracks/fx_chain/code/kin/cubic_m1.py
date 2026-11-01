"""fx_chain item 4 (Claude): the M1 cubic in the dominant pair-changing kernel.

The network's channel budget (m1_kin.py --budget) shows that ss -> X Xbar carries > 99.9 % of the pair-number change
at every temperature down to freeze-out.  Its tree amplitude (tracks/kin_rates/code/xsec22.py, read-only import)
    M = -2 g^2 - lam_s3 lam3/(s - m_s^2) - lam3^2 [1/(t - m_X^2) + 1/(u - m_X^2)]
contains the scalar cubic lam_s3 = V'''(v).  Untilted: 6 lP v.  M1: V''' (v) = 6 (lP - eps) v (-6.67 %), with V''(v),
m_s, lam3 = 2 g^2 v and the quartics unchanged [D, REVIEW_M1].  This script recomputes the MB rate density
gamma_ss(T) with both values on the kin_rates temperature grid and stores the ratio, which m1_kin.py applies as a
T-dependent factor on the ss_XXb channel.  The cubic also enters the s self-reactions (F, N, L, G2, ...), which carry
< 1e-3 of the pair change (budget) and are left unscaled.  [C]
"""
import sys, os, json
import numpy as np
TR = '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks'
sys.path.insert(0, TR + '/kin_rates/code')
import model as M
import xsec22 as X

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'results')
RAT = json.load(open(TR + '/kin_rates/results/rates_interface.json'))
EPS = 0.002


def main():
    P0 = X.Params()
    P1 = X.Params(lam_s3=6 * (M.LP - EPS) * M.V)
    tab = RAT['channels_22']['ss_XXb']['thermal_MB']
    Ts = np.array(tab['T_over_mX']); gtab = np.array(tab['gamma_sub'])
    rows = []
    for t, gt in zip(Ts, gtab):
        T = t * M.MX
        g0 = X.thermal_rate('ss_XXb', T, P0); g1 = X.thermal_rate('ss_XXb', T, P1)
        rows.append(dict(T_over_mX=float(t), gamma_table=float(gt), gamma_untilted=float(g0), gamma_M1=float(g1),
                         ratio=float(g1 / g0), table_check=float(g0 / gt) if gt > 0 else None))
        print('t=%-8.4g table=%.5e own=%.5e (x%.6f)  M1=%.5e  ratio=%.5f' % (t, gt, g0, g0 / gt if gt > 0 else 0, g1, g1 / g0), flush=True)
    # threshold amplitude decomposition (s = 4 m_s^2, t = u = m_X^2 - m_s^2 ... evaluated at the kinematic threshold)
    ms2, mx2 = P0.ms2, P0.mX2
    s = 4 * ms2; t = mx2 + ms2 - 2 * (np.sqrt(ms2) * np.sqrt(ms2) - 0.0)    # threshold: k_i = 0, E1 = m_s, E3 = m_s
    t = ms2 + mx2 - 2 * np.sqrt(ms2) * np.sqrt(ms2)
    parts = dict(contact=-2 * P0.g2, s_channel_untilted=-P0.lam_s3 * P0.lam3 / (s - ms2),
                 s_channel_M1=-P1.lam_s3 * P1.lam3 / (s - ms2), tu_exchange=-P0.lam3 ** 2 * 2 / (t - mx2))
    parts['total_untilted'] = parts['contact'] + parts['s_channel_untilted'] + parts['tu_exchange']
    parts['total_M1'] = parts['contact'] + parts['s_channel_M1'] + parts['tu_exchange']
    parts['ratio_M2_threshold'] = (parts['total_M1'] / parts['total_untilted']) ** 2
    print('threshold amplitude parts', parts)
    json.dump(dict(rows=rows, threshold_amplitude=parts, lam_s3_untilted=P0.lam_s3, lam_s3_M1=P1.lam_s3),
              open(os.path.join(OUT, 'cubic_m1.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
