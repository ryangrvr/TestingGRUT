#!/usr/bin/env python3
"""inf_reheat SKEPTIC (Claude): collect the headline recomputed numbers into ../results/sk_summary.json.
Run after sk_tmodel.py and sk_channels.py.  Labels as in SKEPTIC.md."""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / 'results'
t = json.load(open(RES / 'sk_tmodel.json'))
c = json.load(open(RES / 'sk_channels.json'))
out = dict(
    static_K_check=[dict(b=r['b'], Gamma=r['Gamma'], rel=r['rel']) for r in t['static_check']],
    universal_C=t['universal_C'],
    tmodel=[{k: r[k] for k in ('alpha', 'ns_track', 'ns_MS', 'r_track', 'r_MS', 'H_end_track', 'H_end_mine', 'N_star_track',
                              'N_star_mine_from_trackK', 'K_ratio', 'x0_rolling_over_He2', 'Gamma_shift', 'Phi1_track', 'Phi1_mine')}
            for r in t['rows']],
    cmb=c['cmb'],
    rolling_selfconsistent=c['rolling_selfconsistent'],
    sq_rolling=c['sq_rolling'],
    gk=c['gk'], width=c['width'], optionB=c['optionB'],
    s1_odd=[dict(alpha=o['alpha'], odd_slope_over_tree_muMN=o['mu_MN']['ratio_to_tree'],
                 odd_slope_over_tree_muyphi=o['mu_yphistar']['ratio_to_tree'], VCW_over_V=o['V_CW_even_over_Vtree_at_phistar_mu_yphi'],
                 A_ind_range=[o['A_ind']['A_no_log'], o['A_ind']['A_log_MP']], A_max_no_portal=o['A_max_inflation_branch_no_portal'],
                 lam_min_phi_slice=o['lam_min_phi_slice']) for o in c['s1_odd']],
    s1_induced_portal=c['s1_induced_portal'],
    induced_sigma_leak=c['induced_sigma_leak'],
    portal_annihilation=[dict(alpha=p['alpha'], rows=[{k: r[k] for k in ('lam_track', 'q', 'Gamma_eff_over_Hi', 'm2_induced', 'T_inst')}
                                                     for r in p['rows']]) for p in c['portal']],
    s2_eft=c['s2_eft'],
    pauli=c['pauli'],
    p_feedback=c['p_feedback'],
)
json.dump(out, open(RES / 'sk_summary.json', 'w'), indent=1, default=float)
print('wrote', RES / 'sk_summary.json')
