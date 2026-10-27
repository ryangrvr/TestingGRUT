"""en_particles SKEPTIC: collect the headline numbers from skeptic/results/*.json and Codex's FIELD_CURVES.npz
into skeptic/results/sk_summary.json. Usage: python3 -I sk_summarize.py <skeptic/results> <FIELD_CURVES.npz>"""
import sys, json
import numpy as np
R = sys.argv[1]; J = lambda n: json.load(open(f'{R}/{n}'))
out = dict(source='en_particles_skeptic; own code on reproduction-type histories (Codex engine, unmodified) and on Codex published arrays')
snap = J('sk_snapshots.json')['rows']
tau = lambda s: 20+8.869218738806*s+1.9857962750634526*s*s/2
out['snapshots_tau700_800'] = [dict(tau=round(tau(r['s']), 3), massless=r['massless']['nminus_over_Q'], hartree=r['hartree']['nminus_over_Q'],
                                    ratio=r['massless']['nminus_over_Q']/r['hartree']['nminus_over_Q'], mH2=r['mH2'], Q0_over_Q=r['Q0_over_Q'],
                                    outside_mean=r['outside_mean_fraction']) for r in snap if 699 < tau(r['s']) < 801][:6]
d = np.load(sys.argv[2]); cur = {}
for tag in ['N32', 'N48', 'N64', 'halfstep48', 'seed2_48', 'cut10_48', 'box16_64']:
    t = d[tag+'_tau']; m = (t >= 700) & (t <= 800.01); row = {}
    for key in ['energy_outside_spatial_mean_fraction', 'spatial_mean_charge_fraction', 'phi_plus_interaction_energy_fraction']:
        v = d[tag+'_'+key][m]; row[key] = dict(min=float(v.min()), max=float(v.max()), mean=float(v.mean()), at800=float(d[tag+'_'+key][-1]))
    cur[tag] = row
out['codex_field_curves_tau700_800'] = cur
w = J('sk_W_late_full.json')['rows'][0]
out['late_window_estimators'] = {k: dict(nonzero=w[k]['nminus_nonzero'], zero=w[k]['nminus_zero'], total=w[k]['nminus_all']) for k in ['hann', 'dpss', 'rect', 'virial']}
out['late_window_dpss_groups'] = {k: w['dpss'][k] for k in ['nminus_k_le_1p6', 'nminus_k_1p6_4', 'nminus_k_ge_4']}
out['late_window_centroid_m2_eff'] = [[r['k'], r['nminus'], r['m2_from_centroid_neg']] for r in w['dpss_shells'][:14]]
ms = J('sk_mass_scan_late.json'); out['mass_scan_late'] = ms
subs = []
for nm in ['W_first', 'W_early', 'W_late', 'W_ext']:
    for r in J(f'sk_sub100_{nm}.json')['rows']:
        subs.append(dict(window=nm, s=r['s'], hann=r['hann']['nminus_all'], virial=r['virial']['nminus_all'], dpss=r['dpss']['nminus_all'],
                         hann_nonzero=r['hann']['nminus_nonzero'], virial_nonzero=r['virial']['nminus_nonzero'],
                         lowest_shell_dpss=r['dpss']['nminus_shell1']))
out['subwindow_totals'] = dict(rows=subs, hann_range=[min(x['hann'] for x in subs), max(x['hann'] for x in subs)],
                               virial_range=[min(x['virial'] for x in subs), max(x['virial'] for x in subs)],
                               dpss_range=[min(x['dpss'] for x in subs), max(x['dpss'] for x in subs)])
out['orbit_tau20_exact'] = J('sk_orbit_tau20.json')
out['box_seed_tau547_800_virial'] = {n: dict(nonzero=J(f)['nminus_nonzero_virial'], zero=J(f)['nminus_zero_virial'],
                                            total=J(f)['nminus_nonzero_virial']+J(f)['nminus_zero_virial'],
                                            frac_k_le_1p6=J(f)['frac_k_le_1p6'], lowest_shell_k=J(f)['lowest_shell_k'],
                                            frac_lowest_shell=J(f)['frac_lowest_shell'])
                                     for n, f in [('L12_seed1', 'sk_box_L12_late.json'), ('L12_seed2', 'sk_seed2_late.json'),
                                                  ('L16_N64', 'sk_box_L16_late.json'), ('L12_extension_tau1100_1460', 'sk_box_L12_ext.json')]}
b = J('sk_box16_tau800.json')['rows'][0]
out['L16_regenerated_tau800_vs_codex'] = dict(ours=dict(massless=b['massless']['nminus_over_Q'], hartree=b['hartree']['nminus_over_Q'], Q0_over_Q=b['Q0_over_Q'],
                                                       outside=b['outside_mean_fraction']),
                                              codex=dict(massless=2.412004466176687, hartree=1.265916637314318, Q0_over_Q=1.1468525754940613, outside=0.6493096376292042))
rp = J('sk_response_tau700.json')['shells']
out['response_tau700'] = dict(s0=J('sk_response_tau700.json')['s0'], Z_plus_range=[min(r['Z_plus_rho'] for r in rp), max(r['Z_plus_rho'] for r in rp)],
                              Z_minus_range=[min(r['Z_minus_rho'] for r in rp), max(r['Z_minus_rho'] for r in rp)],
                              rho_m2_kfd4_k4_7p5=[r['m2_rho_peak_kfd4'] for r in rp if 4 <= r['k'] <= 7.5],
                              pencil_m2_k3_8p5=[r['m2_pencil_pos'] for r in rp if 3 <= r['k'] <= 8.5],
                              pencil_resid_K4_low_k=[r['pencil'][4]['resid'] if 4 in r['pencil'] else r['pencil']['4']['resid'] for r in rp[:3]])
out['synthetic_response_bias'] = J('sk_synth_response.json')
out['export'] = J('sk_export.json'); out['phi_quanta'] = J('sk_phi_late.json')
json.dump(out, open(f'{R}/sk_summary.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ['snapshots_tau700_800', 'late_window_estimators', 'late_window_dpss_groups', 'box_seed_tau547_800_virial', 'L16_regenerated_tau800_vs_codex']}, indent=1))
print('subwindow ranges', out['subwindow_totals']['hann_range'], out['subwindow_totals']['virial_range'], out['subwindow_totals']['dpss_range'])
print('codex curves N48', cur['N48']); print('all-control envelope outside', min(c['energy_outside_spatial_mean_fraction']['min'] for c in cur.values()), max(c['energy_outside_spatial_mean_fraction']['max'] for c in cur.values()),
      'phi+int', min(c['phi_plus_interaction_energy_fraction']['min'] for c in cur.values()), max(c['phi_plus_interaction_energy_fraction']['max'] for c in cur.values()),
      'Q0/Q', min(c['spatial_mean_charge_fraction']['min'] for c in cur.values()), max(c['spatial_mean_charge_fraction']['max'] for c in cur.values()))
print('response', {k: out['response_tau700'][k] for k in ['Z_plus_range', 'Z_minus_range', 'rho_m2_kfd4_k4_7p5', 'pencil_m2_k3_8p5', 'pencil_resid_K4_low_k']})
