#!/usr/bin/env python3
"""dq_heat (Claude, own code): collect the headline numbers into ../results/summary.json."""
import json
from pathlib import Path
import numpy as np

RES = Path(__file__).resolve().parent.parent / 'results'
L = lambda f: json.load(open(RES / f))
bg, bgp, ht, vm, vb, su, fk, fs = (L('background.json'), L('background_planck.json'), L('heat_table.json'), L('vector_modes.json'),
                                   L('vector_band.json'), L('sep_universe.json'), L('heat_table.json')['finite_k'], L('finf_sq.json'))
B = ht['budgets']['portals_vanish_at_MY_benchmark']; P = ht['budgets']['portals_vanish_at_MP']
cases = {(c['variant'], c['e']): c for c in vm['cases']}
out = dict(
    release=dict(NX=B['NX'], EX=B['EX'], eps=B['eps'], E_over_Ec=B['W0'], K=B['K'], Rr=bg['validation']['Rr_end'],
                 NYL=B['NYL_hom'], NYH=B['NYH_hom'], Nsig=B['Nsig_hom'], F_inf=bg['prep']['F_inf'], F_vac=240.0,
                 bitwise_vs_fx_holder=bg['validation']['rel_diff_vs_fx_holder']),
    release_planck_bc=dict(NX=P['NX'], EX=P['EX'], eps=P['eps'], E_over_Ec=P['W0'], NYL=P['NYL_hom'], NYH=P['NYH_hom'], Nsig=P['Nsig_hom'],
                           F_inf=bgp['prep']['F_inf'], NX_vs_skeptic=bgp['validation']['rel_diff_NX_vs_skeptic']),
    mediators=B['mediators'], deposit_per_YL_in_EX=B['deposit_per_YL_EX'],
    vector_scaling=vm['scaling_1e-6_to_1e-5'],
    vector_checks=dict(gate1_reproduction='inst variant with spot F_inf reproduces gate1 n_k to 1e-6 (5 modes)',
                       full_over_inst_e1e6=cases[('full', 1e-6)]['N_A'] / cases[('inst', 1e-6)]['N_A'], K_to_1p5=B['K'] ** 1.5,
                       full_over_noX_e1e6_k_le_1=cases[('full', 1e-6)]['N_A_k_le_1'] / cases[('noX', 1e-6)]['N_A_k_le_1'],
                       frozen_angle_ratio=(cases[('noX', 1e-6)]['N_A_k_le_1'] and (bg['prep']['F2_inf'] / (bg['prep']['F2_inf'] - 2 * bg['prep']['amps'][0] ** 2))),
                       band_k_ge_1_full=vb['sets']['full_1e-06']['N_A_k_ge_1'], band_k_ge_1_noX=vb['sets']['noX_1e-06']['N_A_k_ge_1'],
                       ratio_to_GMR_formula=[r['ratio_to_GMR_formula'] for r in B['rows']]),
    pieceB_separate_universe=[dict(e=r['e'], N_A_B=r['N_A_pieceB_SU'], ratio=r['ratio_B_over_A']) for r in su['per_e']],
    finite_k=dict(YS_N_YL_k_le_15=fk['YS_r11']['N_YL_k_le_15'], XYS_N_YL_k_le_15=fk['XYS_r11']['N_YL_k_le_15'],
                  N_YL_k_subtraction=B['N_YL_k_subtraction'], N_YL_all_k_le_80_plateau=B['N_YL_k80_plateau'],
                  cascade_linear_before_cut=B['cascade_linear_before_cut'], N_YH_k=B['N_YH_k'], N_sigma_k=B['N_sigma_k']),
    table={}, e_crit=dict(benchmark=B['e_crit'], planck_bc=P['e_crit']), inflaton_at_exit=fs)
for lab, b in (('benchmark_MY', B), ('planck_MP', P)):
    out['table'][lab] = [dict(e=r['e'], R_vec=r['R_vec'], R_vec_k_le_1=r['R_vec_k_le_1'], R_vec_k_ge_1=r['R_vec_k_ge_1_release_era'],
                              R_vec_SU=r['R_vec_pieceB_SU'], R_Y_hom=r['R_Y_hom'], R_Y_k=r['R_Y_k'], R_Y_k_low=r['R_Y_k_low'],
                              R_cascade=r['R_cascade'], R_total=r['R_total'], W=r['W'], S_per_Q=r['S_per_Q'], margin_bulk=r['margin_bulk'],
                              N_A=r['N_A'], m_A=r['m_A_over_H'], Gamma_A=r['Gamma_A_over_H'], a_mean_vec=r['a_mean_vec'], a_end_vec=r['a_end_vec'],
                              deposit_per_vector=r['deposit_per_vector_EX'],
                              timing={k: r[k] for k in r if k.startswith('timing_')}) for r in b['rows']]
json.dump(out, open(RES / 'summary.json', 'w'), indent=1, default=float)
for r in out['table']['benchmark_MY']:
    print('%.1e R_vec %.4f R_Yh %.4f R_Yk %.1e R_c %.4f R_tot %.4f S/Q %.3f' % (r['e'], r['R_vec'], r['R_Y_hom'], r['R_Y_k'], r['R_cascade'], r['R_total'], r['S_per_Q']))
print(json.dumps(out['vector_checks'], indent=1, default=float))
