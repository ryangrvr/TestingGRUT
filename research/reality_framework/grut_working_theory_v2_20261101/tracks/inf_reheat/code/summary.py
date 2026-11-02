#!/usr/bin/env python3
"""inf_reheat (Claude, own code): headline numbers -> ../results/summary.json (reads the other results files)."""
import json, math, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rh_common as rc

L = lambda n: json.load(open(rc.RES / (n + '.json')))
eng, flu, cha, fm, bm = L('engine_runs'), L('fluid'), L('channels'), L('fermion_modes'), L('benchmarks')

out = dict(note='inf_reheat headline numbers; labels in NOTES.md. H_i units unless stated.')
out['engine_repro'] = dict(NX=eng['repro']['NX_400_1000'], target=eng['repro']['targets']['NX'], NYL=eng['repro']['NYL'],
                           EX=eng['repro']['EX'], D_K=eng['repro']['D_K'])
out['engine_gscan_b3'] = [dict(Gamma=r['Gamma'], NX=r['NX'], EX=r['EX'], NYL=r['NYL'], D=r['D'], K=r['K'],
                               Y_over_YQ_at_HM1=r['yield_over_YQ_at_HM1']) for r in eng['gscan']]
out['engine_slow_decay_vs_b'] = [dict(b=r['b'], NX=r['NX'], EX=r['EX'], NYL=r['NYL'], NYH=r['NYH'], eps=r['eps']) for r in eng['bscan']]
out['engine_Gamma_to_0_b3'] = [dict(Gamma=r['Gamma'], NX=r['NX'], EX=r['EX']) for r in eng['early']]
out['engine_plateau_at_rest_b3'] = [dict(wT=r['wT'], alpha=r['wT'] ** 2 / 6, NX=r['NX'], EX=r['EX'], NYL=r['NYL']) for r in eng['tmodel']]
out['K_law'] = dict(formula='K = c_K(b) (Gamma_I/H_i)^(-1/3), D = K^(3/2), Gamma_I <= 1e-6 (spread < 1e-6)',
                    c_K={k: v['c_K'] for k, v in flu['fits'].items()},
                    validation_vs_engine=[dict(Gamma=v['Gamma'], rel=v['rel']) for v in flu['validation_b3']],
                    N_eq_minus_two_thirds_ln=dict((g['b'], g['N_eq'] - (2 / 3) * math.log(1 / g['Gamma'])) for g in flu['grid']
                                                  if g['Gamma'] == 1e-8))
b3 = flu['fits']['3.0']['c_K']
NX0 = [r for r in eng['early'] if r['Gamma'] == 0.0][0]['NX']
out['yield_law_b3'] = dict(formula='Y/Y_Q(H_M1) = 0.991023 (|N_X|/8047.84) (2.79837/D)', N_X_slow=abs(NX0),
                           coefficient=0.991023 * abs(NX0) / 8047.84 * rc.K_REF ** 1.5 / b3 ** 1.5,
                           H_req_formula='H_i = 1.291e11 GeV x %.4f x (Gamma_I/H_i)^(-1/3)' % ((8047.84 / abs(NX0) * b3 ** 1.5 / rc.K_REF ** 1.5) ** (2 / 3)),
                           r_formula='r = 2.70e-7 (H_i/1.291e11 GeV)^2')
out['F0'] = bm['F0']
out['widths'] = cha['widths']['Nf_y2_for_Gamma_over_m_1_6']
out['optionB'] = {k: cha['optionB'][k] for k in ('A_for_Gamma_over_m_1_6', 'A_GeV', 'mH2_during_inflation_H2', 'tachyonic_growth_exponent_per_half_cycle',
                                                  'tree_quartic_shift', 'CW_over_V_inflation', 'CW_slope_over_tree_slope', 'A_max_inflation_branch_GeV',
                                                  'A_max_resonance_GeV', 'at_A_max_resonance', 'at_A_max_inflation')}
out['blocking_engine_amplitude'] = [{k: r[k] for k in ('y', 'Gamma0_over_m', 'yPhi_over_half_m', 'open_fraction', 'a_unblock',
                                                         'GK_energy_fraction_per_crossing')} for r in cha['blocking']['rows']]
out['GK_check'] = [{k: c[k] for k in ('yPhi', 'm', 'max_abs_dev_n1_vs_GK', 'N_ratio', 'E_ratio', 'max_n_after_2')} for c in fm['cases']]
out['seesaw'] = dict(rows=[{k: r[k] for k in ('M_over_Hi', 'ynu2', 'Gamma_N_over_Hi', 'tau_N_over_Hi')} for r in cha['seesaw']['rows']],
                     M_for_tau_8=cha['seesaw']['M_for_tau_N_8_over_Hi'])
out['hidden_leak'] = dict(thresholds=cha['hidden']['thresholds'], rows=[{k: r[k] for k in ('b', 'g_I', 'Gamma_hidden')} for r in cha['hidden']['rows']],
                          R_inf=cha['hidden']['R_inf_formula'])
tab = []
for b_ in bm['benchmarks']:
    s = b_['solution']
    if s is None:
        tab.append(dict(family=b_['family'], channel=b_['channel'], solution=None, reason=b_['reason']))
        continue
    row = dict(family=b_['family'], channel=b_['channel'], m_I_over_Hi=s['b'], H_i=s['H_i'], H_star=s['H_star'], r=s['r'], ns=s['ns'],
               N_star=s['N_star'], Gamma_over_Hi=s['Gamma_over_Hi'], Gamma_over_m=s['Gamma_over_m'], T_RH=s['T_RH_GeV'], K=s['K'], D=s['D'],
               N_matter_era=s['N_matter_era'], NYL=s['engine_at_b']['NYL'], hidden_leak_R=s['hidden_leak']['R_inf'],
               SQ_e1e6_E=s['SQ_e1e-6']['SQ'], R_vec_E=s['SQ_e1e-6']['R_vec'], R_Y_E=s['SQ_e1e-6']['R_Y'], e_crit_bulk_E=s['e_crit_bulk'],
               plasma_induced_tau=s['plasma']['tau_induced'], plasma_thermal_tau=s['plasma']['tau_thermal'], lambda_IH=s['lambda_IH'])
    if b_['channel'] == 'S1':
        row.update(y=s['y'], y_nat=s['y_nat'], tuning=s['flatness_tuning'], yPhi1_over_half_m=s['yPhi1_over_half_m'], tau_unblock=s['tau_unblock'],
                   plasma_no_IP_tau=s['plasma_no_instant_preheating']['tau_induced'])
    if b_['channel'] in ('S1', 'S2'):
        row.update(M_N=s['M_N_GeV'], tau_N_lab=s['tau_N_lab'], z_RH=s['z_RH'], thermal_lepto_eps_needed=s['thermal_lepto']['eps_needed'],
                   eps_DI=s['eps_DI'])
    if b_['channel'] == 'S2':
        row.update(y_axial=s['y_axial'], f_GeV=s['f_GeV'], Phi1_over_f=s['Phi1_over_f'], axial_NP=s['axial_nonperturbative'])
    if b_['channel'] == 'S3':
        row.update(g_inv_GeV=s['inv_g_GeV'], xi_A=s['xi_A'])
    tab.append(row)
out['benchmarks'] = tab
json.dump(out, open(rc.RES / 'summary.json', 'w'), indent=1, default=float)
print(json.dumps(out['yield_law_b3'], indent=1), json.dumps(out['K_law']['c_K']))
