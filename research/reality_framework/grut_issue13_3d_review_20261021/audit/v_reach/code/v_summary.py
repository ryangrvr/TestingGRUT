"""Verifier v_reach: Richardson-combine the two-resolution Q-ball minimisations, locate the binding threshold,
check monotonicity of Omega(Q) and the archived values, and write e(Q) for v_cosmo.py."""
import json, sys
import numpy as np
R = sys.argv[1]
raw = json.load(open(f'{R}/v_qball_min_raw.json'))
lg = json.load(open(f'{R}/v_qball_large.json'))
rows = {}
for h in ('0.08', '0.04'):
    for br in ('down', 'up'):
        for o in raw[h][br]:
            rows.setdefault(round(o['Q_phys'], 1), {})[h] = o
for k, o in lg.items():
    h, Qb = k.split('_')
    rows.setdefault(round(100 * float(Qb), 1), {})[h] = dict(e=o['e'], Omega=o['Omega'], R90E=o['R90'])
tab = []
for Q in sorted(rows):
    d = rows[Q]
    if '0.08' in d and '0.04' in d:
        eR = d['0.04']['e'] + (d['0.04']['e'] - d['0.08']['e']) / 3
        OR = d['0.04']['Omega'] + (d['0.04']['Omega'] - d['0.08']['Omega']) / 3
        tab.append(dict(Q=Q, e_h008=d['0.08']['e'], e_h004=d['0.04']['e'], e_rich=eR, Omega_rich=OR, R90_h004=d['0.04']['R90E']))
Qs = np.array([t['Q'] for t in tab]); es = np.array([t['e_rich'] for t in tab]); Om = np.array([t['Omega_rich'] for t in tab])
i = np.where(np.diff(np.sign(es - 1.0)) != 0)[0][0]
Qb = float(np.exp(np.interp(0.0, [es[i + 1] - 1, es[i] - 1], [np.log(Qs[i + 1]), np.log(Qs[i])])))
arch = [t for t in tab if abs(t['Q'] - 39944) < 1][0]
# thin-wall estimate with sigma_Q = 0.68, n* = 1.778, m* = 0.56234
Rtw = 3 * 0.68 / (1.778 * (1 - 0.56234))
Qtw = 100 * 4 * np.pi / 3 * 1.778 * Rtw**3
mc = 5.042203313927454  # (E/Q) r [GeV] from v_cosmo
summ = dict(table=tab, binding_threshold_Q=Qb, thin_wall_threshold_Q=float(Qtw),
            Omega_monotone_decreasing_in_Q=bool(np.all(np.diff(Om) < 0)), Omega_range=[float(Om.min()), float(Om.max())],
            Q_range=[float(Qs.min()), float(Qs.max())],
            archived=dict(e_rich=arch['e_rich'], Omega_rich=arch['Omega_rich'], target_e=0.91108652),
            e_at={f'Q={t["Q"]:.4g}': t['e_rich'] for t in tab if t['Q'] in (39944.0, 100000.0, 618530.0, 1519910.0)},
            r_allDM_v1TeV={f'Q={t["Q"]:.4g}': mc / (100 * t['e_rich']) for t in tab})
json.dump(summ, open(f'{R}/v_qball_summary.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != 'table'}, indent=1))
for t in tab:
    print(round(t['Q']), round(t['e_h008'], 7), round(t['e_h004'], 7), round(t['e_rich'], 7), round(t['Omega_rich'], 5), round(t['R90_h004'], 3))
