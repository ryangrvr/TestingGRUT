#!/usr/bin/env python3
"""Collect the verifier's numbers into results/vh_summary.json (no new physics; interpolation of Langevin rows)."""
import json, os
import numpy as np

R = os.path.join(os.path.dirname(__file__), '..', 'results')
ld = lambda n: json.load(open(os.path.join(R, n)))
p1, p2, p3, p4, p5 = (ld('vh1_portals.json'), ld('vh2_rates.json'), ld('vh3_history.json'), ld('vh4_desitter.json'),
                      ld('vh5_reheating.json'))
pc = 0.3116077
rows = p4['langevin_fresh_start']
H = np.array([r['H'] for r in rows])


def cross(key):
    p = np.array([r[key] for r in rows])
    i = np.where(p >= pc)[0][0]
    f = (pc - p[i - 1]) / (p[i] - p[i - 1])
    return float(np.exp(np.log(H[i - 1]) + f * np.log(H[i] / H[i - 1])))


S = {
    'portals': {'max_mixed_beta_at_zero_portals': p1['zero_portal_max_mixed_beta'],
                'analytic_vs_tensor_max_rel': p1['analytic_vs_tensor_max_rel'],
                'lebedev_eq15_rel_diff': p1['lebedev_rel_diff']},
    'rates': {'beta_E': p2['beta_E_4dof'], 'pair_rate_per_dof': p2['pair_rate_per_dof']['total'],
              'kappa': p2['kappa'], 'central_charge_ratios': p2['central_charge_ratios'],
              'ours_over_clery_per_reaction': p2['ours_over_clery'],
              'ours_over_clery_per_particle_if_R_is_number_source': 2 * p2['ours_over_clery'],
              'ours_over_bernal_alpha': p2['ours_over_bernal'],
              'mean_T_over_E': p2['mc_mean_T_over_E'], 'mean_E_over_T': p2['mc_mean_E_over_T']},
    'history': {'TRH_max': p3['TRH_max'], 'TRH_max_band': p3['TRH_max_Tc_pm30'],
                'v_exponent': p3['scaling_exponent_v'], 'H_inst': p3['H_inst_reheat_at_TRHmax'],
                'TRH_max_shift_if_betaE/0.7242': p3['TRH_max_sensitivity_to_betaE_x1/0.7242'],
                'TRH_Pi': p3['TRH_Pi_threshold'], 'TRH_Pi_hubble': p3['TRH_Pi_hubble_damped'],
                'Gamma22_over_H': p3['Gamma22_over_H_at_Thmax'], 'n_h_max': p3['n_h_at_Thmax'],
                'Bh_max_1e10': p3['Bh_max']['1e+10'], 'Tmax_kinetic_selfconsistent': p3['Tmax_selfconsistent_kinetic'],
                'inflaton_dominated_Tmax_max': {str(r['He']): r['Tmax'] for r in p5['rows']},
                'inflaton_dominated_slope_He': p5['dlogTmax_dlogHe']},
    'desitter': {'restoration': p4['restoration'], 'thresholds': p4['thresholds'],
                 'finite_N_restoration_H_over_v': p4['finite_N_restoration_H_over_v'],
                 'langevin_rows': rows,
                 'percolation_pc_crossing_GeV': {k: cross(k) for k in
                                                 ('pminus_fluct_N50', 'pminus_fluct_N60', 'pminus_frozen_N50',
                                                  'pminus_frozen_N60')},
                 'langevin_note': 'fresh start phi=+v, X=D=0; 60000 paths; dN=0.01; Euler-Maruyama; MC error ~0.002'}
}
with open(os.path.join(R, 'vh_summary.json'), 'w') as f:
    json.dump(S, f, indent=1, default=float)
print(json.dumps(S['desitter']['percolation_pc_crossing_GeV'], indent=1))
print(json.dumps(S['history'], indent=1, default=float))
