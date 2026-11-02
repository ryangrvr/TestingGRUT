#!/usr/bin/env python3
"""inf_chain (Claude, own code): dark-matter calibration of the quiet T-model members (TB1: alpha = 0.3, TB2: alpha = 1;
S1 slice and holder) by the inflaton width.  For fixed alpha the potential fixes H_i = H(S1) through A_s, so the
calibration fixes Gamma_I instead of H_i:
  1. the chain (chain.analyse) on the fluid backgrounds at Gamma_I = 1e-10 ... 1e-3 H_i gives H_i,req(Gamma) at e = 1e-6
     (free X: M = 1 + 2p(S/Q)); log-log interpolation solves H_i,req(Gamma_cal) = H_i(potential);
  2. one consistency step: T_RH = (90/pi^2 g*)^{1/4} (Gamma M_P)^{1/2}, N_* with a w = 0 era from eps_H = 1 to T_RH
     (inf_potential calibrate.nstar, read-only), V0 re-normalised to A_s at that N_* (exact Mukhanov-Sasaki), new
     H_i = H(S1), Gamma_cal re-solved (the dimensionless release does not depend on V0);
  3. n_s, r, alpha_s at the pivot, H_* and the channel couplings Gamma = y^2 m_I/16 pi (Majorana, S1/S2 of inf_reheat).
Output ../results/calib_<KEY>.json.  Labels: [C] on the stated model choices; p(S/Q) [C on a conditional recipe / E].
"""
import json, math, sys
import numpy as np
import ch_common as ch
import chain as cn
from benchmarks_reg import get_bench

GAMS = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]


def solve_gamma(rows, Hi):
    lg = np.log([r['gam'] for r in rows]); lH = np.log([r['H_req'] for r in rows])
    # H_req decreases with Gamma
    if not (lH.min() <= math.log(Hi) <= lH.max()):
        return None
    o = np.argsort(lH)
    return float(math.exp(np.interp(math.log(Hi), lH[o], lg[o])))


def main(key):
    alpha = {'TB1': 0.3, 'TB2': 1.0}[key]
    rows = []
    for g in GAMS:
        o = cn.analyse(key, g)
        rows.append(dict(gam=g, H_req=o['relic_e1e6']['H_i_req_GeV'], S_per_Q=o['e1e6']['S_per_Q'], p=o['relic_e1e6'].get('p'),
                         K=o['K'], Y_over_YQ=o['Y_over_YQ_at_HM1'], R_vec=o['e1e6']['R_vec'], R_Y=o['e1e6']['R_Y'],
                         SQ_inf=o['S_per_Q_e_inf'], e_crit=o['e_crit_bulk']))
    b0 = get_bench(key)
    out = dict(method=__doc__, key=key, alpha=alpha, scan=rows, steps=[])
    Hi = b0.Hi; Nst = b0.Nstar
    for step in range(2):
        g = solve_gamma(rows, Hi)
        gG = g * Hi                                         # GeV
        TRH = (90 / (math.pi ** 2 * ch.GSTAR)) ** 0.25 * math.sqrt(gG * ch.MP)
        m = b0.m; bg = b0.bg
        f, fN = ch.ic.bg_at(bg, 1.0); hs, _ = ch.ic.h_of(m, f, fN)
        V0 = (Hi / ch.MP / hs) ** 2
        Nn, info = ch.cal.nstar(m, bg, V0, ('wT', 0.0, TRH), hs)
        sp = ch.ic.spectrum_at_pivot(m, bg, Nn)
        V0n = ch.ic.AS_TARGET / sp['Pz']
        Hin = math.sqrt(V0n) * hs * ch.MP
        fs, fNs = ch.ic.bg_at(bg, Nn); hst, _ = ch.ic.h_of(m, fs, fNs)
        Hstar = math.sqrt(V0n) * hst * ch.MP
        mI = math.sqrt(m.m2_min() * V0n) * ch.MP          # GeV
        # interpolated chain quantities at g
        lg = np.log([r['gam'] for r in rows])
        at = lambda k: float(math.exp(np.interp(math.log(g), lg, np.log([r[k] for r in rows]))))
        y = math.sqrt(16 * math.pi * g * Hi / mI)
        st = dict(step=step, H_i_target_GeV=Hi, Gamma_cal_over_Hi=g, Gamma_over_mI=g * Hi / mI, T_RH_GeV=TRH, Nstar_new=Nn,
                  V0_new=V0n, H_i_new_GeV=Hin, H_star_GeV=Hstar, m_I_GeV=mI, m_I_over_Hi=mI / Hin, ns=sp['ns'], r=sp['r'],
                  alpha_s=sp['alpha_s'], N_matter_efolds=math.log(info['rho_end'] / (math.pi ** 2 / 30 * ch.GSTAR * TRH ** 4)) / 3,
                  K=at('K'), S_per_Q_e1e6=at('S_per_Q'), p=at('p'), M=1 + 2 * at('p'), S_per_Q_e_inf=at('SQ_inf'),
                  R_vec_e1e6=at('R_vec'), y_majorana_N1=y)
        out['steps'].append(st)
        print('%s step %d: target H_i=%.4e -> Gamma_cal=%.3e H_i (%.2e m_I), T_RH=%.3e GeV, N*=%.3f, H_i(new)=%.4e, H*=%.4e, n_s=%.4f r=%.3e K=%.1f S/Q=%.1f p=%.2f y=%.2e N_MD=%.1f' % (
            key, step, Hi, g, st['Gamma_over_mI'], TRH, Nn, Hin, Hstar, sp['ns'], sp['r'], st['K'], st['S_per_Q_e1e6'], st['p'], y, st['N_matter_efolds']), flush=True)
        Hi = Hin
    json.dump(out, open(ch.RES / ('calib_%s.json' % key), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
