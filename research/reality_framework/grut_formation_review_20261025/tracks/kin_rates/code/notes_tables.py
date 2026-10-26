#!/usr/bin/env python3
"""kin_rates: print the markdown tables of NOTES.md sections 3.3, 4.3 and 5 from results/rates_interface.json and results/balance.json."""
import json
import numpy as np
RES = '../results'
I = json.load(open(f'{RES}/rates_interface.json')); B = json.load(open(f'{RES}/balance.json'))
KT = (0.3, 1.0, 3.0, 7.344, 12.0)


def fmt(x, e=None):
    if x is None:
        return '–'
    if e is not None and e > 0:
        return f'{x:+.2e} ± {e:.0e}'
    return f'{x:.3e}'


print('### n->2 thermal table (MB, mu=0): gamma_sub [GeV^4] (resolved*) / gamma_full_reg')
print('| key | ' + ' | '.join(f'T = {t} m_X' for t in KT) + ' |')
print('|---|' + '---|' * len(KT))
for ch, c in I['channels_n2'].items():
    tb = c['thermal_MB']
    cells = []
    for t in KT:
        i = int(np.argmin(np.abs(np.array(tb['T_over_mX']) - t)))
        g, e, gf, res, bnd = tb['gamma_sub'][i], tb['gamma_sub_err'][i], tb['gamma_full_reg'][i], tb['remainder_resolved'][i], tb['remainder_bound'][i]
        cells.append(f"{g:+.2e}{'*' if res else ''} / {gf:.1e}" if abs(g) > 0 else f"0 / {gf:.1e}")
    print(f"| {ch.split('_')[0]} ({ch.split('_',1)[1]}) | " + ' | '.join(cells) + ' |')
print()
print('### 2->n thermal table (MB): gamma_sub / gamma_full_reg')
print('| channel | ' + ' | '.join(f'T = {t} m_X' for t in KT) + ' |')
print('|---|' + '---|' * len(KT))
for ch, c in I['channels_2n'].items():
    tb = c['thermal_MB']
    cells = []
    for t in KT:
        i = int(np.argmin(np.abs(np.array(tb['T_over_mX']) - t)))
        cells.append(f"{tb['gamma_sub'][i]:+.2e}{'*' if tb['remainder_resolved'][i] else ''} / {tb['gamma_full_reg'][i]:.1e}")
    print(f"| {ch} | " + ' | '.join(cells) + ' |')
print()
print('### NR limits: A_eff/A_old')
for ch, d in B['nr_limits'].items():
    print(ch, 'A_old', d['A_old'], ' '.join(f"T={r['T_over_mX']}: sub {r['ratio_sub_to_old']:.3f} (fullreg {r['ratio_fullreg_to_old']:.3g})" for r in d['rows']))
print()
print('### table-level detailed balance: ratio n2/2n (full_reg, sub) per pair at T')
worst_full = 0; worst_sub = []
for row in B['detailed_balance']:
    if row.get('ratio_full') is None:
        continue
    rf = row['ratio_full']; rs = row.get('ratio_sub')
    worst_full = max(worst_full, abs(rf - 1))
    # sub: significance using both errors
    g1, e1, g2, e2 = row['gamma_n2_sub'], row['gamma_n2_sub_err'], row['gamma_2n_sub'], row['gamma_2n_sub_err']
    sig = (g1 - g2) / np.hypot(e1, e2) if (e1 + e2) > 0 else 0
    worst_sub.append((abs(sig), row['channel'], row['T_over_mX'], rf, rs, g1, e1, g2, e2))
    if row['T_over_mX'] in (0.3, 3.0, 12.0):
        print(f"{row['channel']:14s} T={row['T_over_mX']:5.1f} full ratio {rf:.4f}  sub: {g1:+.2e}±{e1:.0e} vs {g2:+.2e}±{e2:.0e} (ratio {rs if rs is None else round(rs,3)}, {sig:+.1f}σ)")
print('worst |full ratio - 1|:', worst_full)
ws = sorted(worst_sub, reverse=True)[:5]
print('largest sub deviations:', [(round(w[0], 1), w[1], w[2]) for w in ws])
print()
print('### old s route')
for r in B['old_s_route']:
    if r['T_over_mX'] in (0.03, 0.3, 1.0, 3.0, 7.344, 12.0):
        print(f"T={r['T_over_mX']}: gamma_old {r['gamma_old_threshold_const']:.3e} gamma_sub {r['gamma_ss_XXb_sub']:.3e} ratio {r['ratio_sub_to_old']:.4f}")
print()
print('### regulator spread per channel at T (gamma_sub 0.02 vs 0.04)')
for row in B['detailed_balance']:
    if row['T_over_mX'] in (0.3, 3.0, 12.0):
        print(f"{row['channel']:14s} T={row['T_over_mX']:5.1f} sub02 {row['gamma_n2_sub']:+.2e}±{row['gamma_n2_sub_err']:.0e} sub04 {row.get('gamma_n2_sub_reg0.04')} RIS_frac {row['ris_fraction_at_regulator']}")
print()
print('### s self-reactions F, N and chain thresholds')
for ch in ('F_sss_ss', 'N_ssss_ss', 'L_sss_XXb'):
    tb = I['channels_n2'][ch]['thermal_MB']
    print(ch, [(t, f"{tb['gamma_sub'][int(np.argmin(np.abs(np.array(tb['T_over_mX'])-t)))]:.2e}") for t in KT])
print({ch: c['chain_threshold_sqrt_s_over_mX'] for ch, c in I['channels_n2'].items()})
