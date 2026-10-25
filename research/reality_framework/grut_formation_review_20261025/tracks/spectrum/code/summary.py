#!/usr/bin/env python3
"""Collect the spectrum-track numbers into results/spectrum_summary.json and print markdown tables."""
import json, os
import numpy as np
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
L = lambda n: json.load(open(os.path.join(RES, n)))
F1, F2, F3 = L('F1_analysis.json'), L('F2dt_analysis.json'), L('F3h_analysis.json')
F4 = L('F4long_analysis.json') if os.path.exists(os.path.join(RES, 'F4long_analysis.json')) else None
C = L('codex_trace_spectra_REPRODUCTION_DATA.json')
kb = F1['k_bins']; Q0, E0 = F1['Q0'], F1['E0']
out = dict(Q0=Q0, E0=E0, k_bins=kb)

# 1. numerical convergence at a common observable: R=240, window 12000
conv = {}
for n, d in [('F1', F1), ('F2dt', F2), ('F3h', F3)]:
    w = d['fixed_radius']['240.0']['windows']['12000']
    conv[n] = dict(source='independent', h=d['h'], dt=d['dt'], Q60_9600_12000=d['retention']['Q60_mean_9600_12000'], NX=w['NX'], NXbar=w['NXbar'],
                   EX=w['EX'], EXbar=w['EXbar'], Ns=w['Ns'], Es=w['Es'], NX_k_ge_005=sum(v for k, v in zip(kb[:-1], w['X_bins']) if k >= 0.05),
                   NXbar_spatial_T12000_Rw240=d['spatial']['T12000_Rw240']['NXbar'],
                   Q0_minus_Qin240_T12000=d['fixed_radius']['240.0']['time_domain'].get('Q0_minus_Qinside_T') if d['T'] == 12000 else None)
for n in ['P025_time', 'P0125', 'P00625']:
    r = C['runs'][n]; w = r['windows']['12000']
    conv[n] = dict(source='Codex arrays (reproduction data), this track\'s analysis', Q60_9600_12000=r['core60_late_mean_fraction_codex'],
                   NX=w['NX'], NXbar=w['NXbar'], EX=w['EX'], EXbar=w['EXbar'], Ns=w['Ns'], Es=w['Es'],
                   NX_k_ge_005=sum(v for k, v in zip(kb[:-1], w['X_bins']) if k >= 0.05), Qpass_240_T12000=r['Qpass_T_continuum_flux'])
out['convergence_R240_W12000'] = conv

# 2. certified antiparticle and neutral numbers (converged runs F2dt, F3h; window/radius spread from F1)
best = [conv['F2dt'], conv['F3h']]
xb = [c['NXbar'] for c in best] + [c['NXbar_spatial_T12000_Rw240'] for c in best]
fr1 = F1['fixed_radius']
xb_window = [fr1['240.0']['windows'][w]['NXbar'] for w in ['12000', '18000', '24000']] + [fr1['480.0']['windows']['24000']['NXbar'], fr1['960.0']['windows']['24000']['NXbar']] + [F1['spatial'][f'T24000_Rw{R}']['NXbar'] for R in [120, 240, 360]]
out['antiparticles'] = dict(NXbar_converged_mean=float(np.mean(xb)), NXbar_converged_spread=float(np.ptp(xb)),
                            NXbar_F1_window_radius_method_values=xb_window, NXbar_F1_rel_spread=float(np.ptp(xb_window) / np.mean(xb_window)),
                            mean_omega=float(np.mean([c['EXbar'] / c['NXbar'] for c in best])),
                            Xbar_bins_F3h_R240_W12000=F3['fixed_radius']['240.0']['windows']['12000']['Xbar_bins'],
                            fraction_of_Q0=float(np.mean(xb) / Q0))
out['neutral'] = dict(Ns=float(np.mean([c['Ns'] for c in best])), Es=float(np.mean([c['Es'] for c in best])),
                      Ns_spread_with_codex_fine=float(np.ptp([conv[n]['Ns'] for n in ['F2dt', 'F3h', 'P0125', 'P00625']])))

# 3. particle spectrum: window/radius/method comparison at T=24000 (F1) and T=48000 (F4long)
def cmp_table(d, T, Rs):
    fr, sp = d['fixed_radius'], d['spatial']
    cols = {f'fr{int(R)}_W{T}': fr[f'{R:.1f}']['windows'][str(T)]['X_bins'] for R in Rs if str(T) in fr[f'{R:.1f}']['windows']}
    for Rw in [120, 240, 360]:
        key = f'T{T}_Rw{Rw}'
        if key in sp: cols['sp_' + key] = sp[key]['X_bins']
    return cols
out['X_bins_F1_T24000'] = cmp_table(F1, 24000, [240., 480., 960.])
out['X_bins_F1_T12000'] = cmp_table(F1, 12000, [240., 480., 960.])
if F4:
    out['X_bins_F4_T48000'] = cmp_table(F4, 48000, [240., 480., 960.])
    out['X_bins_F4_T24000'] = cmp_table(F4, 24000, [240., 480., 960.])

def certify(cols, tol):
    arr = np.array(list(cols.values()))
    rel = np.ptp(arr, axis=0) / np.maximum(np.mean(arr, axis=0), 1e-30)
    return [float(x) for x in rel]
out['X_bins_rel_spread_F1_T24000'] = certify(out['X_bins_F1_T24000'], 0)
if F4: out['X_bins_rel_spread_F4_T48000'] = certify(out['X_bins_F4_T48000'], 0)

# 4. budget and uncertified component at the longest window
def budget(d, T):
    assert abs(d['T'] - T) < 1, 'budget uses end-of-run inside charges and fluxes; only valid at the run end'
    fr = d['fixed_radius']; w240 = fr['240.0']['windows'][str(T)]; w480 = fr['480.0']['windows'][str(T)]; w960 = fr['960.0']['windows'][str(T)]
    qin = d['Qinside_final']
    lowk = sum(v for k, v in zip(kb[:-1], w240['X_bins']) if k < 0.03)
    midk = sum(v for k, v in zip(kb[:-1], w240['X_bins']) if 0.03 <= k < 0.05)
    return dict(T=T, NX_R240=w240['NX'], NXbar_R240=w240['NXbar'], EX_R240=w240['EX'], EXbar_R240=w240['EXbar'], Ns_R240=w240['Ns'], Es_R240=w240['Es'],
                NX_R480=w480['NX'], NXbar_R480=w480['NXbar'],
                closure_R480_charge=(w480['NX'] - w480['NXbar']) / w480['tapered_time_domain_Qpass'] - 1,
                closure_R480_energy=(w480['EX'] + w480['EXbar'] + w480['Es']) / fr['480.0']['time_domain']['Epass_T'] - 1,
                incoming_fraction_R480=w480['NX_incoming'] / w480['NX'],
                closure_R960_charge=(w960['NX'] - w960['NXbar']) / w960['tapered_time_domain_Qpass'] - 1,
                closure_R960_energy=(w960['EX'] + w960['EXbar'] + w960['Es']) / fr['960.0']['time_domain']['Epass_T'] - 1,
                incoming_fraction_R960=w960['NX_incoming'] / w960['NX'],
                X_lowk_lt_0p02_R240=sum(v for k, v in zip(kb[:-1], w240['X_bins']) if k < 0.02),
                X_lowk_lt_0p03_R240=lowk, X_k_0p03_0p05_R240=midk,
                Q_core_r30=qin['30'], Q_core_r60=qin['60'], Q_inside_240=qin['240'],
                halo_60_240=qin['240'] - qin['60'], halo_30_240=qin['240'] - qin['30'])
out['budget_F1_T24000'] = budget(F1, 24000)
if F4: out['budget_F4_T48000'] = budget(F4, 48000)
out['nonlinearity_F1'] = {R: F1['fixed_radius'][R]['nonlinearity'] | {k: F1['fixed_radius'][R]['windows']['24000'][k] for k in ['k10_X', 'k50_X', 'fracX_k2_below_10dm2slow', 'phase_beyond_R_at_k10']} for R in ['240.0', '480.0', '960.0']}
out['retention'] = {n: d['retention'] for n, d in [('F1', F1), ('F2dt', F2), ('F3h', F3)] + ([('F4long', F4)] if F4 else [])}
out['core_frequency_F1'] = F1['core_frequency_from_phase_at_r5']
if F4: out['core_frequency_F4'] = F4['core_frequency_from_phase_at_r5']
json.dump(out, open(os.path.join(RES, 'spectrum_summary.json'), 'w'), indent=1)

# markdown
print('| run | source | Q60 9600-12000 | N_X | N_Xbar | E_X | E_Xbar | N_s | E_s | N_X(k>=0.05) |')
print('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|')
for n, c in conv.items():
    print(f"| {n} | {'indep' if c['source']=='independent' else 'Codex arrays'} | {100*c['Q60_9600_12000']:.4f}% | {c['NX']:.1f} | {c['NXbar']:.2f} | {c['EX']:.1f} | {c['EXbar']:.2f} | {c['Ns']:.2f} | {c['Es']:.1f} | {c['NX_k_ge_005']:.1f} |")
for key in ['X_bins_F1_T24000'] + (['X_bins_F4_T48000'] if F4 else []):
    cols = out[key]; names = list(cols)
    print(f'\n{key}\n| k bin | ' + ' | '.join(names) + ' | rel spread |'); print('|---' * (len(names) + 2) + '|')
    rel = certify(cols, 0)
    for i in range(len(kb) - 1):
        print(f'| {kb[i]:.3f}-{kb[i+1]:.3f} | ' + ' | '.join(f'{cols[n][i]:.1f}' for n in names) + f' | {rel[i]:.3f} |')
print('\nXbar bins F3h R240 W12000:', [round(x, 3) for x in out['antiparticles']['Xbar_bins_F3h_R240_W12000']])
print(json.dumps({k: out[k] for k in ['antiparticles', 'neutral', 'budget_F1_T24000'] + (['budget_F4_T48000'] if F4 else [])}, indent=1))
