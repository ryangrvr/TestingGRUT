#!/usr/bin/env python3
"""mf_chain: assemble results/summary.json (headline numbers) from the run outputs.
Usage: python3 -B summary.py"""
import json, glob, os
import mfc_common as M
import table

if __name__ == '__main__':
    rows = table.main([])
    fl = json.load(open(M.RES / 'final_ledger.json'))['cases']
    its = {os.path.basename(f)[8:-5]: json.load(open(f)) for f in glob.glob(str(M.RES / 'iterate_*.json'))}
    cmb = json.load(open(M.RES / 'cmb_rescale.json'))
    kap = json.load(open(M.RES / 'cmb_kappa_scan.json'))
    grav = {os.path.basename(f)[5:-5]: {k: dict(N_total=v['N_total'], E_total_over_EX=v['E_total_over_EX'])
                                       for k, v in json.load(open(f))['fields'].items()} for f in glob.glob(str(M.RES / 'grav' / 'grav_*.json'))}
    out = dict(
        track='mf_chain (Claude, multi-field inflaton pass; GitHub issue #13 context)',
        labels='[C] given the model choices unless stated; hidden-heat leak and post-switch gravitational tail [E]; ledger, closed form, '
               'p(S/Q) map inherited [Q]. Every history, rescaling rule and criterion is a model choice; nothing here validates or excludes anything.',
        engine='mf_ideal rpmfi.c MODE 5 + vmodes.c + drivers, imported read-only and unchanged (mfc_common.py redirects outputs)',
        regression=dict(note='HFa as exported, raw R', S_per_Q_1e6_ours=[r for r in rows if r['tag'] == 'HFa_raw'][0]['SQ_1e6'],
                        S_per_Q_1e6_mf_ideal_smoke=14.693300033331639),
        history_builder_validation=dict(note='unblocked friction, t_sw = 0.6, design H_i: rebuilt history vs bench_HFa',
                                        K=[r for r in rows if r['tag'] == 'val_track'][0]['K'], K_bench=2.0481855289913504,
                                        S_per_Q_1e6=[r for r in rows if r['tag'] == 'val_track'][0]['SQ_1e6']),
        all_runs=rows, final_ledger=[dict(key=c['key'], label=c['label'], H_i_GeV=c['H_i_GeV'], K=c['K'], eps=c['eps'], E_over_Ec=c['E_over_Ec'],
                                          N_YL=c['N_YL'], N_YH=c['N_YH'], N_sigma=c['N_sigma'], S_per_Q_e_inf=c['S_per_Q_e_inf'],
                                          e_crit=c['e_crit_bulk'], e1e6={k: c['by_e']['e1e-6'][k] for k in ('S_per_Q', 'branch', 'M', 'p', 'H_i_req_GeV', 'R_vec', 'R_Y', 'R_casc', 'dR_grav')},
                                          e3e6={k: c['by_e']['e3e-6'][k] for k in ('S_per_Q', 'branch', 'M', 'H_i_req_GeV')},
                                          leak_e1e6={k: dict(R_leak=v['R_leak'], S_per_Q=v['S_per_Q'], M=v['M'], H_req_over_H_i=v['H_req_over_H_i'])
                                                     for k, v in c['by_e']['e1e-6']['with_hidden_heat'].items()},
                                          spectator=c['spectator'], isocurvature=c.get('isocurvature'), map=c['mf_ideal_map_same_K_ms3'])
                                     for c in fl],
        iterations={k: dict(converged=v.get('converged'), steps=[dict(H_i_GeV=s['H_i_GeV'], K=s['K'], S_per_Q_1e6=s.get('SQ_1e6', s.get('SQ_leak')),
                                                                      H_req_GeV=s['H_req_GeV'], ratio=s['ratio']) for s in v['steps']])
                    for k, v in its.items()},
        cmb_at_converged=cmb['rows'], cmb_regression=cmb['regression_design'], kappa_scan=kap,
        gravitational=grav)
    M.dump(M.RES / 'summary.json', out)
    print('wrote', M.RES / 'summary.json')
