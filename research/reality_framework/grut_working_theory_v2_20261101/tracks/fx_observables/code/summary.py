"""fx_observables: collect the headline numbers into results/summary.json (no new physics)."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results')
vt = json.load(open(os.path.join(RES, 'vapour_transfer.json')))
rows = []
for r in vt['B']:
    if 'transfer' not in r:
        continue
    k = np.array([t['k'] for t in r['transfer']]); P = np.array([t['P_ratio'] for t in r['transfer']])
    def first_below(x):
        j = np.where(P < x)[0]
        if len(j) == 0:
            return None
        j = j[0]
        if j == 0:
            return float(k[0])
        # log interpolation
        return float(math.exp(np.interp(x, [P[j], P[j - 1]], [math.log(k[j]), math.log(k[j - 1])])))
    rows.append(dict(Q=r['Q'], S_over_Q=r['S_over_Q'], sticking=r['sticking'], zeta=r['sticking'] * r['Q']**(-1 / 3),
                     Tvis_f_keV=r['Tvis_f_GeV'] * 1e6, one_minus_eta_f=r['one_minus_eta_f'], f_v_mass=r['f_v_mass'],
                     v_rms_at_decoupling=r['v_rms_at_m'], m_WDM_equiv_keV=r['m_WDM_equiv_keV'], r_s_kpc=r['r_s_Mpc'] * 1e3,
                     P10=float(P[k == 10][0]), P20=float(P[k == 20][0]), P50=float(P[k == 50][0]), P100=float(P[k == 100][0]),
                     k_at_P_0p99=first_below(0.99), k_at_P_0p95=first_below(0.95), k_at_P_0p25=first_below(0.25)))
rows.sort(key=lambda z: (z['zeta']))
out = dict(note='P = P_M1/P_CDM at z = 0 (sub-horizon Newtonian toy, see vapour_transfer.py); zeta = s Q^-1/3 controls the capture freeze-out',
           scenario_A_codex_gas=vt['A'], scenario_B=rows)
json.dump(out, open(os.path.join(RES, 'summary_transfer.json'), 'w'), indent=1, default=float)
print('%-8s %-5s %-6s %-9s %-8s %-6s %-6s %-9s %-8s %-7s %-7s %-7s %-7s %-7s %-7s' % ('Q', 'S/Q', 's', 'zeta', 'Tf_keV', '1-eta', 'f_v', 'mWDMkeV', 'r_s_kpc', 'P10', 'P20', 'P50', 'P100', 'k99', 'k95'))
for z in rows:
    print('%-8.0e %-5.2f %-6g %-9.2e %-8.3g %-6.3f %-6.3f %-9.3g %-8.3g %-7.4f %-7.4f %-7.4f %-7.4f %-7s %-7s' % (
        z['Q'], z['S_over_Q'], z['sticking'], z['zeta'], z['Tvis_f_keV'], z['one_minus_eta_f'], z['f_v_mass'], z['m_WDM_equiv_keV'],
        z['r_s_kpc'], z['P10'], z['P20'], z['P50'], z['P100'], '%.3g' % z['k_at_P_0p99'] if z['k_at_P_0p99'] else '-', '%.3g' % z['k_at_P_0p95'] if z['k_at_P_0p95'] else '-'))
