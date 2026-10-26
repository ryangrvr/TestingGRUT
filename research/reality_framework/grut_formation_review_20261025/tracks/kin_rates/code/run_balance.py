#!/usr/bin/env python3
"""kin_rates: table-level checks of the n<->2 kernels (own code).  usage: python3 -B run_balance.py <kernels_n2.json> <xsec_2n.json|none> <kernels_22.json> <out.json>
1. MB detailed balance from the tables: gamma_{n->2}(T) from Gt against gamma_{2->n}(T) from sigma for every pair,
   for the regulated-full and the subtracted kernels (same |M|^2 -> must agree up to MC/interpolation error).
2. Late non-relativistic limit against the old threshold constants A_a, A_b of gp_leftovers/amp42.py and the
   threshold sigma v of the old s route (SV_SS of lsrelic.py): A_old = S_in gamma(T)/n_eq^4 as T -> 0.
3. The on-shell (sequential) content removed by the subtraction at the regulator width, per channel and T, and
   the regulator spread (0.02 vs 0.04) of the remainder.
4. The old s route gamma_old = SV_SS n_s^2/2 against the energy-dependent ss -> X Xb rate.
"""
import json
import sys
import numpy as np
import model as M
import tables as TB

PAIRS = {'A_XXXbXb_XXb': 'XXb_XXbXXb', 'B_XXXXb_XX': 'XX_XXXXb', 'C_XXXbXb_ss': 'ss_XXXbXb', 'H_XXbss_XXb': 'XXb_XXbss',
         'J_XXXbs_Xs': 'Xs_XXXbs', 'K_XXbss_ss': 'ss_XXbss', 'N_ssss_ss': 'ss_ssss', 'D_XXXb_Xs': 'Xs_XXXb',
         'E_XXbs_XXb': 'XXb_XXbs', 'E2_XXs_XX': 'XX_XXs', 'F_sss_ss': 'ss_sss', 'G_XXbs_ss': 'ss_XXbs',
         'G2_Xss_Xs': 'Xs_Xss', 'L_sss_XXb': 'XXb_sss'}
A_OLD = dict(A_XXXbXb_XXb=0.0842146549366979, B_XXXXb_XX=0.002361976125302127)   # [Q own gp_leftovers] GeV^-8
SV_SS_OLD = 0.0001260281438027877                                                 # [Q own ls_threshold] GeV^-2
TS = (0.1, 0.3, 1.0, 3.0, 7.344, 12.0)


def main(n2path, twonpath, k22path, out):
    n2 = json.load(open(n2path))['channels']
    twon = json.load(open(twonpath))['channels'] if twonpath != 'none' else {}
    k22 = json.load(open(k22path))
    res = dict(detailed_balance=[], nr_limits={}, ris_content=[], old_s_route=[], regulator_spread=[])
    for ch, c in n2.items():
        tab = TB.Table(c, 'n2')
        tab04 = TB.Table(c, 'n2', freg=0.04) if any(abs(r['freg'] - 0.04) < 1e-12 for r in c['rows']) else None
        rev = PAIRS.get(ch)
        rtab = TB.Table(twon[rev], '2n') if rev in twon else None
        for t in TS:
            T = t * M.MX
            row = dict(channel=ch, reverse=rev, T_over_mX=t)
            for mode in ('full', 'sub'):
                gf = tab.thermal(T, mode)
                row[f'gamma_n2_{mode}'] = gf
                if mode == 'sub':
                    row['gamma_n2_sub_err'] = abs(tab.thermal(T, 'sub', shift=1.0) - gf)
                if rtab is not None:
                    gr = rtab.thermal(T, mode)
                    row[f'gamma_2n_{mode}'] = gr
                    row[f'ratio_{mode}'] = gf / gr if gr else None
                    if mode == 'sub':
                        row['gamma_2n_sub_err'] = abs(rtab.thermal(T, 'sub', shift=1.0) - gr)
            if tab04 is not None:
                g04 = tab04.thermal(T, 'sub')
                row['gamma_n2_sub_reg0.04'] = g04
                row['regulator_spread_rel'] = (g04 - row['gamma_n2_sub']) / (abs(row['gamma_n2_sub']) + row['gamma_n2_sub_err'] + 1e-300)
            row['ris_fraction_at_regulator'] = 1 - row['gamma_n2_sub'] / row['gamma_n2_full'] if row['gamma_n2_full'] else None
            res['detailed_balance'].append(row)
        res['ris_content'].append(dict(channel=ch, chain_threshold_sqrt_s_over_mX=tab.chain_threshold(c),
                                       threshold_sqrt_s_over_mX=tab.sq_th / M.MX, threshold_power_law_exponent=tab.k_th))
        if ch in A_OLD:
            rows = []
            for t in (0.03, 0.1, 0.3, 1.0, 3.0, 7.344, 12.0):
                T = t * M.MX
                n = TB.n_eq(M.MX, T)
                g = tab.thermal(T, 'sub'); gfull = tab.thermal(T, 'full')
                rows.append(dict(T_over_mX=t, A_eff_sub_GeV_minus8=tab.Sin * g / n**4, A_eff_fullreg_GeV_minus8=tab.Sin * gfull / n**4,
                                 ratio_sub_to_old=tab.Sin * g / n**4 / A_OLD[ch], ratio_fullreg_to_old=tab.Sin * gfull / n**4 / A_OLD[ch]))
            # threshold kernel value against |M_thr|^2 Phi_2 = A_old (2m)^4 S_out
            first = sorted([r for r in c['rows'] if abs(r['freg'] - 0.02) < 1e-12], key=lambda r: r['sqrt_s_over_mX'])[0]
            res['nr_limits'][ch] = dict(A_old=A_OLD[ch], rows=rows, first_point_sqrt_s_over_mX=first['sqrt_s_over_mX'],
                                        note='A_eff = S_in gamma(T)/n_eq^4 -> A_old as T -> 0 (old rate density = A_old n^4/S_in)')
    # old s route vs ss -> XXb
    th = k22['thermal']['ss_XXb']['rows']
    for r in th:
        t = r['T_over_mX']; T = t * M.MX
        ns = TB.n_eq(M.MS, T)
        g_old = SV_SS_OLD * ns * ns / 2
        res['old_s_route'].append(dict(T_over_mX=t, gamma_old_threshold_const=g_old, gamma_ss_XXb_sub=r['gamma_sub'], gamma_ss_XXb_nores=r['gamma_nores'],
                                       ratio_sub_to_old=r['gamma_sub'] / g_old, sigmav_sub=r['sigmav_sub'], sigmav_over_SV_SS=r['sigmav_sub'] / SV_SS_OLD))
    json.dump(res, open(out, 'w'), indent=1, default=float)
    for row in res['detailed_balance']:
        if row['T_over_mX'] in (0.3, 3.0, 12.0):
            print(f"{row['channel']:14s} T/m={row['T_over_mX']:5.1f} g_n2 sub={row['gamma_n2_sub']:+.3e}+-{row['gamma_n2_sub_err']:.1e} full={row['gamma_n2_full']:.3e} "
                  f"ratio_full={row.get('ratio_full')} ratio_sub={row.get('ratio_sub')} RIS={row['ris_fraction_at_regulator']:.3f} reg04={row.get('gamma_n2_sub_reg0.04')}")
    for ch, d in res['nr_limits'].items():
        print(ch, [(r['T_over_mX'], round(r['ratio_sub_to_old'], 4)) for r in d['rows']])
    print('old s route ratio:', [(r['T_over_mX'], round(r['ratio_sub_to_old'], 4)) for r in res['old_s_route'] if r['T_over_mX'] in (0.03, 0.3, 1.0, 3.0, 7.344, 12.0)])


if __name__ == '__main__':
    main(*sys.argv[1:5])
