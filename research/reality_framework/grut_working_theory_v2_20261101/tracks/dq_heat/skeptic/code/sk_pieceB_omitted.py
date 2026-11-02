#!/usr/bin/env python3
"""dq_heat SKEPTIC (Claude, own code): (1) piece B against the multi-field linear modes; (2) omitted vector/heat channels.

(1) dq_heat's vector_band.py uses the single-Higgs-direction Proca mode (F^2 = 2(|X|^2 + 9|Y|^2 + 36|Sigma|^2)).  dq_heat's
    own ymodes runs (rpmodes.c, global-U(1) limit e -> 0, exact de Sitter BD start for every eigenmode) contain the
    Sigma-phase eigenmode (m^2 ~ 0, S = 1), i.e. the eaten Goldstone = A_L for k/a >> m_A (Goldstone equivalence).  Its
    number over the same k window is compared with the Proca band (XYS <-> 'full', YS <-> 'noX'); the k = 15-80 part is
    compared with N_YH (cascade Y_H -> Y_L + G products, counted separately in R_cascade).
(2) Order-of-magnitude estimates [E] of channels not in dq_heat's ledger, at e = 1e-6 (deposit 1.7805 E_X per vector):
    a. inverse-decay freeze-in X Xbar -> A while T_h >> m_A:  N = 9 e^2 T_c^5/(192 pi^2 K m_A)  (Maxwell-Boltzmann,
       g_A = 3, int x^3 K_1 = 3 pi/2), T_c = (15 E/pi^2)^{1/4} with E = E_X (1 + R_Y);
    b. Y_L -> 3X + G through a virtual Y_H (unsuppressed vertex sqrt2 mu G Y_L Y_H):  BR ~ (2 mu^2/(m_H^2 - m_L^2)^2)
       Phi_4/Phi_3, Phi_4/Phi_3 = m_L^2/(96 pi^2) (massless phase space; propagator taken at E_G = 0, an upper value);
    c. condensate-induced two-body Y_L -> X X (vertex 3 kappa <X>):  int Gamma dt = 9 kappa^2 chi^2 ln(a_f/a_i)/(8 pi m_L K)
       (|X| = chi/a, chi^2 = late conformal mean) -- these Y_L deposit earlier, i.e. R_Y decreases by about that fraction.
Output: ../results/sk_pieceB_omitted.json.   Labels: [C] (1), [E] (2).
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DQ = HERE.parents[1] / 'results'
RES = HERE.parent / 'results'


def goldstone(fname, tau, windows):
    d = json.load(open(DQ / fname)); K = np.array(d['kgrid']); lk = np.log(K)
    T = [t for t in d['times'] if abs(t['tau'] - tau) < 1e-9][0]
    i = [j for j, c in enumerate(T['eig']) if c['S'] > 0.5 and c['m2'] < 1][0]
    n = np.array(T['n_k'][i]); integ = K ** 3 * n / (2 * np.pi ** 2)
    out = {}
    for lo, hi in windows:
        sel = (K >= lo) & (K <= hi)
        out['%g-%g' % (lo, hi)] = float(np.trapezoid(integ[sel], lk[sel]))
    spec = {('%.3f' % k): float(x) for k, x in zip(K, n) if 1 <= k <= 10}
    return out, spec


def main():
    W = [(1, 15), (1, 30), (15, 80)]
    vb = json.load(open(DQ / 'vector_band.json'))
    band = {}
    for tag in ('full_1e-06', 'noX_1e-06'):
        S = np.array(vb['sets'][tag]['spectrum']); k = S[:, 0]; n = S[:, 1]
        band[tag] = {'%g-%g' % w: float(np.trapezoid((k ** 3 * n / (2 * np.pi ** 2))[(k >= w[0]) & (k <= w[1])],
                                                       np.log(k[(k >= w[0]) & (k <= w[1])]))) for w in W[:2]}
        band[tag]['n_at_4.45_peak'] = float(np.interp(4.45, k, n)); band[tag]['n_at_4.713'] = float(np.interp(4.713, k, n))
    xys, xys_spec = goldstone('ymodes_XYS_r11.json', 40.0, W)
    ys40, _ = goldstone('ymodes_YS_r11.json', 40.0, W)
    ys100, _ = goldstone('ymodes_YS_r11.json', 100.0, W)
    bg = json.load(open(DQ / 'background.json'))
    NYH = bg['validation']['NYH']
    p1 = dict(proca_band=band, multifield_XYS_tau40=xys, multifield_YS_tau40=ys40, multifield_YS_tau100=ys100,
              XYS_goldstone_n_k_1_to_10=xys_spec, N_YH=NYH,
              ratio_XYS_over_full_1_30=xys['1-30'] / band['full_1e-06']['1-30'], ratio_YS100_over_noX_1_30=ys100['1-30'] / band['noX_1e-06']['1-30'])
    print('piece B, k in [1,30]: Proca full %.4f | multi-field XYS(tau 40) %.4f (ratio %.3f) ; Proca noX %.4f | YS(tau 40/100) %.4f/%.4f' % (
        band['full_1e-06']['1-30'], xys['1-30'], p1['ratio_XYS_over_full_1_30'], band['noX_1e-06']['1-30'], ys40['1-30'], ys100['1-30']))
    print('  Proca full n(4.45) = %.3e ; XYS Goldstone n(4.713) = %.3e (Proca full n(4.713) = %.3e) ; XYS k 15-80: %.4f vs N_YH %.4f' % (
        band['full_1e-06']['n_at_4.45_peak'], xys_spec.get('4.713', float('nan')), band['full_1e-06']['n_at_4.713'], xys['15-80'], NYH))

    # (2) omitted channels
    K = bg['late']['K_end']; EX = bg['validation']['EX']; chi2 = bg['late']['chi2_mean_late']
    kap, MY, hh, mu = 0.0034586847211042523, 30.0, 40.0, 7.159456159513793
    mL = math.sqrt(MY ** 2 - math.sqrt(2) * mu * hh); mH = math.sqrt(MY ** 2 + math.sqrt(2) * mu * hh)
    NL = bg['validation']['NYL'] + NYH
    est = {}
    for e in (1e-6, 1e-5):
        mA = 240 * e; G = e * e * mA / (48 * math.pi); dep = mA * math.sqrt(math.pi * K / (2 * G)) / EX
        Tc = (15 * EX * (1 + 0.3788) / math.pi ** 2) ** 0.25
        N_id = 9 * e * e * Tc ** 5 / (192 * math.pi ** 2 * K * mA)
        est['inverse_decay_freeze_in_e%.0e' % e] = dict(T_conf=Tc, N_A=N_id, R=N_id * dep)
    BR4 = 2 * mu ** 2 / (mH ** 2 - mL ** 2) ** 2 * mL ** 2 / (96 * math.pi ** 2)
    est['YL_to_3X_G'] = dict(BR_upper=BR4, N_A=NL * BR4, R_e1e6=NL * BR4 * 1.7805)
    frac = 9 * kap ** 2 * chi2 * math.log(5569 / 3.0) / (8 * math.pi * mL * K)
    est['YL_condensate_two_body'] = dict(fraction_of_YL=frac, dR_Y_upper=-frac * 0.3788)
    for k_, v in est.items():
        print(k_, {a: ('%.3g' % b) for a, b in v.items()})
    json.dump(dict(method=__doc__, pieceB=p1, omitted=est), open(RES / 'sk_pieceB_omitted.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
