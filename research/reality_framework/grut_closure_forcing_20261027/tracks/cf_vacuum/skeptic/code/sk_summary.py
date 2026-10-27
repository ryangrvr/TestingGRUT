"""Collect the skeptic's key numbers and verdict checks into ../results/sk_summary.json (run after sk_line.py and
sk_checks.py)."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
L = json.load(open(os.path.join(RES, 'sk_line.json')))
C = json.load(open(os.path.join(RES, 'sk_checks.json')))
lp = L['line_PDG2026']

# FCC-ee illustration with the track's recipe (precisions from the track: 40 MeV, 4 MeV, 0.00012) [E]
s_fcc = np.sqrt(0.04 ** 2 + (0.495 * 0.004) ** 2 + (0.352 * 0.00012 / 0.0009) ** 2 + 0.113 ** 2 + 0.032 ** 2)
fcc = {'sigma_tot_track_recipe': s_fcc, 'z_direct_central_kept': 1.603 / s_fcc,
       'z if a pole-mass renormalon ambiguity 0.11 GeV (1605.03609) is not cancelled': 1.603 / np.hypot(s_fcc, 0.11),
       'z with 0.25 GeV ambiguity [R, order of larger estimates]': 1.603 / np.hypot(s_fcc, 0.25),
       'z with xsec-like central and 0.11': 1.103 / np.hypot(s_fcc, 0.11)}

c0 = L['c6_operator_B15']['Mbar_Pl']['0']['Mt_crit']
summ = {
    'date': '2026-10-27',
    'by': 'Claude skeptic agent cf_vacuum_skeptic (closure-forcing pass, issue #13); git read-only, no commit',
    'line_PDG2026_Mt_crit': {'SMDR (track)': 170.997, 'B15 (Bednyakov 2015 matching + own 3-loop RGE)': lp['B15_loop3']['Mt_crit'],
                             'B13 (Buttazzo 2013 matching + own 3-loop RGE)': lp['B13_loop3']['Mt_crit'],
                             'B15 slopes (per GeV Mh, per 0.0009 alpha_s)': [lp['B15_dMt_dMh'], lp['B15_dMt_per_0.0009_as']],
                             'B15 mu_star': lp['B15_loop3']['mu_star'],
                             'loop2 - loop3 (B15)': lp['B15_loop2']['Mt_crit'] - lp['B15_loop3']['Mt_crit']},
    'z_direct': {'SMDR (track)': 3.45, 'B15': L['tension']['B15']['direct']['z'], 'B13': L['tension']['B13']['direct']['z'],
                 'B15 xsec': L['tension']['B15']['xsec']['z'], 'B15 FLAG alpha_s': L['tension']['B15']['direct, FLAG as']['z'],
                 'B15 EW-subfield alpha_s 0.1221(27)': L['tension']['B15']['direct, EW-fit as']['z'],
                 'B15 with Bednyakov upward th 0.17': L['tension']['B15_with_Bednyakov_upward_th_0.17']['direct']['z']},
    'Lambda_lambda_at_172.60_GeV': {'track (SMDR start)': 4.9e10, 'B15': L['Lambda_lambda']['B15']['172.6'],
                                    'B13': L['Lambda_lambda']['B13']['172.6']},
    'c6_shift_of_line_cut_Mbar_Pl': {k: (None if v is None else v['Mt_crit'] - c0)
                                     for k, v in L['c6_operator_B15']['Mbar_Pl'].items()},
    'c6_shift_of_line_cut_M_Pl': {k: (None if v is None else v['Mt_crit'] - c0)
                                  for k, v in L['c6_operator_B15']['M_Pl'].items()},
    'c6_needed_to_stabilise_PDG_central': L['c6_needed_to_stabilise_PDG_central'],
    'FCC_ee_illustrative': fcc,
    'hidden_sector_reproduced': {k: C['hidden_sector_P3a'][k] for k in
                                 ('lamS_min_lamSY0_lamY1', 'lamSY_needed_at_lamS=2.8e-05',
                                  'lamSY_needed_at_lamS=1.8e-10', 'window_point depth / rho_Lambda')},
    'CW_dlambda_Sigma': C['CW_Y_threshold']['numeric'],
    'portal_loop_needed': {'b': C['portal_loop_stabilisation']['b_needed'],
                           'lamHX': C['portal_loop_stabilisation']['lamHX_needed'],
                           'ratio_to_invisible_bound': [C['portal_loop_stabilisation']['ratio b_needed/b_max'],
                                                        C['portal_loop_stabilisation']['ratio lamHX_needed/lamHX_max']]},
    'hidden_landau_pole_GeV': {k: v for k, v in C['hidden_landau_pole'].items() if isinstance(v, dict)},
    'valley_identity_holds': C['valley_identity_holds'],
    'verdict_checks': {
        'P3a-SM': 'not refuted in substance; MIXED is outside the brief taxonomy -> STANDARD (relabelled published '
                  'absolute-stability hypothesis; testable; retrospectively disfavoured 1.6-3.7 sigma)',
        'P3a-GRUTref': 'verdict refuted: P3a + GRUT reference is internally consistent for any M_t <= M_t,crit; the '
                       'conflict is retrospective data tension, not branch inconsistency -> STANDARD (inherits P3a-SM)',
        'P3a-GRUT-Sigma-corner': 'UNTESTABLE (agree); all numbers reproduced',
        'P3d-MPP': 'MIXED -> STANDARD (relabelled Froggatt-Nielsen / Volovik MPP); equality UV-sensitive downward',
        'P3v-volume': 'STANDARD (agree)',
        'missed row, reading (b) with uncalibrated QV1': 'already fails retrospectively (Lambda_ours = 0, QV1 residual '
                                                         'problem); no separate SM content'}}
json.dump(summ, open(os.path.join(RES, 'sk_summary.json'), 'w'), indent=1, default=float)
print(json.dumps(summ, indent=1, default=float))
