#!/usr/bin/env python3
"""fx_holder track (Claude, own code): collects the headline numbers of tri_point.json, tri_runs.json, lam_scan.json and
controls.json into ../results/summary.json (no new computation)."""
import json
from pathlib import Path
R = Path(__file__).resolve().parent.parent / 'results'
tp = json.load(open(R / 'tri_point.json')); tr = json.load(open(R / 'tri_runs.json')); ls = json.load(open(R / 'lam_scan.json'))
ct = json.load(open(R / 'controls.json'))
S = {}
c = tp['coupling']['lam1_b3']
S['coupling'] = dict(g_I_b3=c['g_I'], g_I_b5=tp['coupling']['lam1_b5']['g_I'], alpha_I=c['alpha_I'], phi_I0_GeV_b3=c['phi_I0_GeV'], phi_I0_over_Mr_b3=c['I0_over_Mr'],
                     g_I_phi0_over_H=c['g_I_phi0_over_H'])
for lam in ('lam1', 'lam0.5'):
    sp = tp['stationary_points'][lam]
    S[f'preparation_{lam}'] = dict(amplitudes=sp['amplitudes'], delta=sp['delta'], physical_m2=sp['physical'], eig7_inflaton_shift=sp['eig7_shift_of_inflaton'],
                                   n_distinct_minima=len(sp['global_search_fixed_phiI0']))
ol = tp['one_loop']['lam1']
S['one_loop'] = dict(delta_mI2_MSbar_QMY=ol['delta_mI2_late_vacuum'], tadpole_shift=ol['inflaton_shift_from_tadpole'],
                     beta_lambda_IS=ol['induced_operators']['beta_lambda_IS'], beta_lambda_IY=ol['induced_operators']['beta_lambda_IY'],
                     beta_A_IS=ol['induced_operators']['beta_A_IS'], dmS2_phi0_per_lnQ=ol['induced_operators']['delta_mS2_at_phi0_per_unit_lnQ'],
                     dmY2_phi0_per_lnQ=ol['induced_operators']['delta_mY2_at_phi0_per_unit_lnQ'], two_loop_sunset_1e15=ol['two_loop_sunset_cutoff_1e15GeV'])
S['bounds_lam0.5'] = {k: v for k, v in tp['bounds']['lam0.5'].items() if k != 'late_vacuum_global_search'}
S['decays'] = tp['decays']
S['cp'] = dict(rephasing_residual_max=max(abs(x) for x in tp['cp_and_clock']['rephasing_invariant_residuals']),
               opposite_branch_prep=dict(delta=tp['cp_and_clock']['opposite_branch']['delta'], amps=tp['cp_and_clock']['opposite_branch']['amps'],
                                         lightest=tp['cp_and_clock']['opposite_branch']['lightest']))
runs = {}
for key, summ in tr['runs'].items():
    for k, v in summ.items():
        runs[f'{key}/{k}'] = dict(NX_400_1000=v['NX_400_1000'], NX_640_800=v['NX_640_800'], NX_1000_2000=v.get('NX_1000_2000'), run80=[v['NX_run80_min'], v['NX_run80_max']],
                                  yield_over_YQ=v['yield_over_YQ_400_1000'], H_new_GeV=v['H_new_GeV_400_1000'], H_new_range=v['H_new_GeV_run80_range'],
                                  NYL=v['NYL'], NYH=v['NYH'], R_Y=v['ledger_e1e-7']['R_Y'], R_central_e1e7=v['ledger_e1e-7']['R_central'],
                                  R_central_e1e6=v['ledger_e1e-6']['R_central'], checks=v['checks'])
S['runs'] = runs
S['comparisons'] = tr['comparisons']
S['lambda_scan'] = [dict(lam=r['lam'], prep_lightest=r['prep_lightest'], kappa_over_kappa_max=r['kappa_over_kappa_max'], tri_yield=r['trilinear']['yield_over_YQ'],
                         tri_H_new=r['trilinear']['H_new_GeV'], tri_NYL=r['trilinear']['NYL'], tri_R_Y=r['trilinear']['R_Y'], orig_yield=r['original']['yield_over_YQ'])
                    for r in ls['rows']]
S['controls'] = ct
json.dump(S, open(R / 'summary.json', 'w'), indent=1, default=float)
print(json.dumps(S['comparisons']['lam0.5_b3_g0.5'], indent=1, default=float)[:600])
