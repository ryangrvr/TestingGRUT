#!/usr/bin/env python3
"""inf_potential SKEPTIC: collect the headline numbers of the skeptic runs into ../results/sk_summary.json."""
import json
from pathlib import Path
R = Path(__file__).resolve().parent.parent / 'results'
L = lambda f: json.load(open(R / f))
cmb, rl, mi, kin, born, flo, eng = L('sk_cmb.json'), L('sk_reledger.json'), L('sk_mi_scan.json'), L('sk_kinematics.json'), L('sk_yy_born.json'), L('sk_floquet.json'), L('sk_engine.json')
S = dict(note='inf_potential SKEPTIC headline numbers; labels and context in ../SKEPTIC.md', cmb={}, reledger={}, mI_scan=[], kinematics={}, born_HT6t=[], floquet_first={}, engine_late={})
for n, r in cmb['runs'].items():
    S['cmb'][n] = {k: r[k] for k in ('Nstar_skeptic', 'Nstar_track', 'As_at_trackcal', 'ns', 'alpha_s', 'r', 'n_t', 'ns_SR', 'r_SR', 'H_S1_for_As2p1', 'Hstar_over_HS1', 'm_I_over_Hi')}
    S['cmb'][n]['ns_pull_ACT_DR6_PACTLB'] = (r['ns'] - 0.9743) / 0.0034
    S['cmb'][n]['ns_pull_PACT'] = (r['ns'] - 0.9709) / 0.0038
for n, r in rl['benchmarks'].items():
    S['reledger'][n] = dict(A_S1_e1e6=r['A_S1']['1e-06'], A_S1_ecrit=r['A_S1_ecrit'], B_S1_Kamp_e1e6=r['B_S1_Kamp']['1e-06'], B_ecrit=r['B_S1_Kamp_ecrit'],
                            A_Send_e1e6=r['A_Send']['1e-06'], C_Send_plusYY_e1e6=r['C_Send_plusYY']['1e-06'], D_S1_plusYY_e1e6=r['D_S1_plusYY']['1e-06'],
                            F_track_rough=r['F_track_rough_e1e-6'], E_full={e: v['SQ'] for e, v in r['E_full'].items()}, E_full_e_inf=r['E_full_e_inf'],
                            E0_noYY={e: v['SQ'] for e, v in r['E0_noYY'].items()}, E0_noYY_e_inf=r['E0_noYY_e_inf'],
                            E1_noYY_noYLdeposit={e: v['SQ'] for e, v in r['E1_noYY_noYLdeposit'].items()}, inputs=r['inputs'])
S['mI_scan'] = [{k: x[k] for k in ('m_I', 'K', 'NX', 'EX', 'eps', 'NYH', 'NYL', 'Nsig', 'SQ_e1e-6', 'SQ_e1e-5', 'SQ_e1e-6_noYLdeposit')} for x in mi['rows']]
S['kinematics'] = kin if isinstance(kin, dict) else {}
S['born_HT6t'] = [dict(t_on=x['t_on'], N_YL=x['N_YL'], N_YH=x['N_YH'], track=x['track']) for x in born['runs']]
S['floquet_first'] = {n: v[0] for n, v in flo['runs'].items()}
S['engine_late'] = {k: v['late_mean'] for k, v in eng['runs'].items()}
json.dump(S, open(R / 'sk_summary.json', 'w'), indent=1)
print(json.dumps(S['cmb'], indent=0)[:1500])
