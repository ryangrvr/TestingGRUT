#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code): compact summary of the sk_late.c runs -> results/sk_late.json.
Adopted frozen values: N_inf + c2 a^-2 over a >= a_end/3 of the long run (a_end = 2582); spread = max deviation of the other
fits (a^-2 over a >= a_end/6, a_end/2; a^-1 + a^-2 over a >= a_end/6; const over a >= a_end/2).
Usage: sk_late_summary.py LONG.json CONV.json OUT.json"""
import json, sys
L = json.load(open(sys.argv[1])); C = json.load(open(sys.argv[2]))
fz = {}
for nm in ('NX', 'NS', 'Econf', 'Eeng'):
    F = {(f['model'], round(f['a_min'], 1)): f for f in L['fits'][nm]}
    keys = sorted({k[1] for k in F})
    a6, a3, a2 = keys
    ad = F[('a-2', a3)]
    alts = [F[('a-2', a6)], F[('a-2', a2)], F[('a-1,a-2', a6)], F[('const', a2)]]
    fz[nm] = ad['N_inf']; fz[nm + '_err'] = (ad['err'] ** 2 + max(abs(f['N_inf'] - ad['N_inf']) for f in alts) ** 2) ** 0.5
fz['gauss_NX_plus_6NS'] = fz['NX'] + 6 * fz['NS']
wl = {(w['eta_lo'], w['eta_hi']): w for w in L['windows']}
conv = max(abs(w['NX'] - wl[(w['eta_lo'], w['eta_hi'])]['NX']) for w in C['windows'])
o = dict(method='sk_late.c: conformal fields, exact radiation a(eta), 4th-order Yoshida velocity Verlet, h = c_h/(40 a); start: rpdq state at H_i t = 1000',
         long_run=dict(info=L['info'], end=L['end'], c_h=L['c_h'], eta_end=L['eta_end'], fits=L['fits'], std_slope_vs_a=L['std_slope_vs_a'],
                       windows=L['windows']),
         step_convergence=dict(c_h=[C['c_h'], L['c_h']], max_window_NX_diff=conv, end=C['end']),
         vs_rpdq_same_windows=C['compare_rpdq'], vs_rpdq_matched_19p75_windows=C['compare_rpdq_matched_windows'],
         frozen=fz,
         dq_entropy=dict(NX=-8047.840302004927, NX_err=0.2194639130876233, NS=1341.3067184118302, EX=188737.70931829128, EX_err=1.10164))
json.dump(o, open(sys.argv[3], 'w'), indent=1)
print(json.dumps(fz, indent=1)); print('step convergence', conv)
print('matched-window max |dNX|', max(abs(m['d']) for m in C['compare_rpdq_matched_windows']))
