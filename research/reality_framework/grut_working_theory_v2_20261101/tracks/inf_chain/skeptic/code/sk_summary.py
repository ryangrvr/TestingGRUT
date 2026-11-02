#!/usr/bin/env python3
"""inf_chain SKEPTIC (Claude, own code): collect the skeptic's results into ../results/sk_summary.json."""
import json, pathlib
R = pathlib.Path(__file__).resolve().parent.parent / 'results'
L = lambda n: json.load(open(R / n))
cmb, vm, gs, hg, lg, es = (L('sk_cmb_fluid.json'), L('sk_vmodes.json'), L('sk_gamma_switch.json'), L('sk_higgs_TB.json'),
                           L('sk_ledger.json'), L('sk_estimates.json'))
out = dict(
    cmb_TB={k: dict(Nstar=cmb[k]['Nstar'], H_i=cmb[k]['H_i_GeV'], H_star=cmb[k]['H_star_GeV'], ns=cmb[k]['ns'], r=cmb[k]['r'],
                    alpha_s=cmb[k]['alpha_s'], x0=cmb[k]['x0'], C_K=cmb[k]['C_K'], K=cmb[k]['K'], T_RH=cmb[k]['T_RH_sudden_GeV'],
                    rel_to_track=cmb[k]['rel']) for k in ('TB1', 'TB2')},
    K_gamma05={k: cmb['K_' + k] for k in ('HT6t', 'T1')},
    vectors_TB1_static={str(r['e']): dict(mine=r['N_A_k_le_1_static'], track=r['track_N_A_k_le_1_static'], ratio=r['ratio_mine_over_track'])
                        for r in vm['rows']} | {'slope': vm['slope_dlnN_dlne']},
    gamma_switch={k: v.get('rel_change') for k, v in gs.items() if 'rel_change' in v} | {'HT6t_Send_NX': [gs['HT6t_Send']['gamma_from_S1']['NX'], gs['HT6t_Send']['gamma_from_end']['NX']]},
    higgs_inflation={k: {mt: v['xi_EW_min'] for mt, v in r['by_Mt'].items()} for k, r in hg['rows'].items()},
    closed_form_asymptote=lg['asymptote'],
    ledger_rel={k: v['rel'] for k, v in lg['ledger'].items()},
    R_Y_threshold={k: dict(R_Y=v['R_Y'], threshold=v['R_Y_threshold_bulk']) for k, v in lg['ledger'].items()},
    asym_share={k: v['asym_share_of_DM'] for k, v in lg['ledger'].items()},
    leak_MP_own={k: dict(factor=v['factor'], R_leak=v['R_leak_MP_own'], S_per_Q=v['S_per_Q_with_leak_MP_own']) for k, v in lg['leak'].items()},
    iteration=lg['iteration'],
    estimates={k: dict(np_per_crossing=v['frac_per_crossing_np'], pert_per_period=v['frac_per_period_pert'], F_np=v['F_np_until_unblock'],
                       N_unblock=v['N_unblock_after_S1'], R_grav_X=v['R_grav_X'], R_grav_pi=v['R_grav_pi_e1e6'], A_ind=v['A_ind_GeV'],
                       A_bound=v['A_bound_GeV']) for k, v in es.items()},
)
json.dump(out, open(R / 'sk_summary.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=0, default=float)[:3000])
