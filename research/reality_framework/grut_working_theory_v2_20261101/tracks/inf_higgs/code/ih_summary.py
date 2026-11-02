#!/usr/bin/env python3
"""Track inf_higgs: collect the headline numbers into ../results/summary.json (and R < 0 window statistics)."""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np

import ih_common as C
import ih_run as IR

R_ = C.RES


def load(n):
    return json.load(open(R_ / n))


def gk(b, N):
    for k in b:
        if abs(float(k) - N) < 1e-6:
            return b[k]


def windows_table():
    out = {}
    xi = 0.0877
    bgs = {'M1r (3, 0.5)': C.Background(None, gam=0.5, b=3.0, tau_end=10.0)}
    for nm in C.BENCH:
        bgs[nm] = C.Background(C.load_bench(nm), gam=0.0, Nb0=2.0, tau_end=10.0)
    for k, bg in bgs.items():
        w = bg.windows(nmax=6, xi=xi)
        F = bg.fine()
        t, R = F['tau'], F['R']
        post = (t > bg.tau_endinf) & (t < bg.tau_endinf + 3.0)
        frac = float(np.mean(R[post] < 0))
        out[k] = dict(first_windows=w[:4], fraction_of_time_R_negative_first_3_Hinv=frac,
                      sum_exponent_first_3_Hinv=float(sum(x['exponent_xi'] for x in bg.windows(nmax=400, xi=xi)
                                                          if x['tau_start'] < bg.tau_endinf + 3.0)),
                      tau_endinf=bg.tau_endinf)
    return out


def main():
    S = {}
    B = load('bench_D.json')
    tab = {}
    for k, o in B.items():
        p = o['PDG']; N = o['N_obs']
        tab[k] = dict(N_obs=N, xi_H=o['xi_H'], H_GeV=o['H_GeV'],
                      D60={m: gk(p[m]['best'], 60)['D'] for m in p},
                      DNobs={m: gk(p[m]['best'], N)['D'] for m in p},
                      tau_c60={m: gk(p[m]['best'], 60).get('tau_c') for m in p},
                      D60_variants_at_172_86=p['172.86']['variants_D60'],
                      Mt_crit_60=o['Mt_crit_60'], Mt_crit_Nobs=o['Mt_crit_Nobs'], Mt_crit_60_criterion_band=o['Mt_crit_60_band'],
                      Mt_crit_60_Mh={m: v['Mt_crit_60'] for m, v in o['Mh_variants']['rows'].items()},
                      deadline60_tau_from_S1={m: v['deadline_60'] for m, v in o['Mt_grid'].items() if float(m) >= 172.85},
                      sigma_units_crit60=(o['Mt_crit_60'] - 172.57) / 0.29 if o['Mt_crit_60'] else None)
    S['gravitational_only_xiEW_0.1'] = tab
    xs = {}
    for f in sorted(R_.glob('xiscan_*.json')):
        d = json.load(open(f))
        xs[f.stem.replace('xiscan_', '')] = d['edges_xiEW']
    S['lower_edge_postinflation_xiEW'] = xs
    up = {}
    for f in sorted(R_.glob('upper_*.json')):
        d = json.load(open(f))
        up[f.stem.replace('upper_', '')] = {k: dict(xi_EW=v['xi_EW'], runaway=v['runaway'],
                                                     two_sigma_prod_over_hmax={m: v[m]['two_sigma_prod_over_hmax'] for m in ('172.28', '172.57', '172.86')},
                                                     D60={m: v[m]['D60'] for m in ('172.28', '172.57', '172.86')}) for k, v in d.items()}
    S['upper_edge'] = up
    inf = load('inflation_HN.json')
    S['inflation_xi_min'] = {k: {m: {kk: vv for kk, vv in o[m].items() if 'xi_min' in kk} for m in ('172.28', '172.57', '172.86')}
                             | dict(H_end_over_Hi=o['H_end_over_Hi'], R_over_H2_end=o['R_over_H2_end'],
                                    B_tree_MD_end_vs_const={m: [o[m]['tree MD']['B_end'], o[m]['tree MD']['B_constH']] for m in ('172.28', '172.57', '172.86')},
                                    log10N_tree_MD_HN={m: o[m]['tree MD']['log10N_HN'] for m in ('172.28', '172.57', '172.86')},
                                    log10N_loop_MD_HN={m: o[m]['loop MD']['log10N_HN'] for m in ('172.28', '172.57', '172.86')})
                             for k, o in inf.items()}
    S['inflation_gaussian_tail'] = {k: {m: v[m]['best_log10_patches'] for m in v} for k, v in load('inflation_gaussian_tail.json').items()}
    S['inflation_Aterm'] = load('inflation_Aterm.json')
    P = load('plasma.json')
    S['plasma'] = {k: dict(y_nat=v['y_nat'], lambda_IH=v['lambda_IH'],
                           rows={kk: {x: rr[x] for x in ('y', 'MN_GeV', 'Gamma_pert_open', 'yPhi1_over_half_mI', 'ncross',
                                                          'tau_m2_1em3', 'tau_m2_1em2', 'tau_thermal_pessimistic', 'tau_unblock',
                                                          'A_ind_GeV', 'A_max_inflation_GeV_Mt172_86', 'm2_at')}
                                 for kk, rr in v['rows'].items()}) for k, v in P.items()}
    S['D_with_plasma'] = load('final_plasmaD.json')
    S['portal'] = load('final_portal.json')
    S['controls'] = load('final_controls.json')
    S['R_negative_windows'] = windows_table()
    S['inputs'] = dict(PDG2024=C.PDG24, xi_EW=0.1, xi_H_at_Hi=tab['HT6t|MD|iH0']['xi_H'], T_RH_slow_for_Nobs=IR.T_RH_SLOW,
                       criterion='higgs_reheat scale-resolved D with the skeptic canonical normalisation; D < 0 survival')
    C.jdump(S, R_ / 'summary.json')
    # console table
    for k, v in tab.items():
        print(k, 'Nobs %.2f' % v['N_obs'], 'D60', {m: round(x, 1) for m, x in v['D60'].items()},
              'crit60 %.3f (%.2f sigma) Nobs %.3f band %s Mh %s' % (v['Mt_crit_60'], v['sigma_units_crit60'], v['Mt_crit_Nobs'],
                                                                   [round(x, 3) for x in v['Mt_crit_60_criterion_band']],
                                                                   {m: round(x, 3) for m, x in v['Mt_crit_60_Mh'].items()}))
    for k, v in S['R_negative_windows'].items():
        w = v['first_windows'][0]
        print(k, 'first R<0: R_min/H^2 %.2f duration %.3f/H exponent %.4f' % (w['R_min_over_H2'], w['duration_H'], w['exponent_xi']),
              'frac %.3f sum exponent (3/H) %.3f' % (v['fraction_of_time_R_negative_first_3_Hinv'], v['sum_exponent_first_3_Hinv']))


if __name__ == '__main__':
    main()
