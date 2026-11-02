#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code): the induced-portal boundary condition 'vanish at M_P' -- which ln(M_P/M_Y)?
Summarises sk_state1000.py runs of rpdq.c (400-1000 windows) at c(M_Y) = -L/(4 pi^2) for
  L = 13.31 (dq_entropy), L = ln(2.435e18/(30 H_M1)) = 13.3121 (dq_heat, fx_holder skeptic),
  L = ln(2.435e18/(30 x 1.0764e11)) = 13.5332 (self-consistent with the free-X H_i of that release).
Usage: sk_portal_summary.py RUNDIR OUT.json"""
import json, math, sys
CC = 2 * math.sqrt(2e-3) * (3e-3) ** -0.75
H_M1 = 1.3428427349130815e11
RD = sys.argv[1]
rows = {}
for t in ('base_1000', 'planck1331', 'planck_dqheat', 'planck_sc'):
    d = json.load(open('%s/state_%s.json' % (RD, t)))
    N, E = abs(d['NX_400_1000']), d['EX_400_1000']
    rows[t] = dict(cport=d['cport'], L=-d['cport'] * 4 * math.pi ** 2, NX_400_1000=N, EX_400_1000=E, E_over_Ec=E / (N / CC) ** (4 / 3),
                   NYL=d['NYL_400_1000'], NYH=d['NYH_400_1000'], Nsig=d['Nsig_400_1000'], Y_over_YQ=d['yield_over_YQ_400_1000'],
                   R_Y_rel=(d['NYL_400_1000'] + d['NYH_400_1000']) / E, END=d['END'])
b, p, q, s = rows['base_1000'], rows['planck1331'], rows['planck_dqheat'], rows['planck_sc']
out = dict(rows=rows,
           dq_heat_vs_dq_entropy_NX=q['NX_400_1000'] / p['NX_400_1000'] - 1,
           dq_heat_reported_NX=10606.16025730126,
           selfconsistent_vs_1331=dict(dlnNX=math.log(s['NX_400_1000'] / p['NX_400_1000']), dL=s['L'] - p['L'],
                                       dlnNX_per_unit=math.log(s['NX_400_1000'] / p['NX_400_1000']) / (s['L'] - p['L']),
                                       dlnH=-(2 / 3) * math.log(s['Y_over_YQ'] / p['Y_over_YQ']),
                                       dE_over_Ec=s['E_over_Ec'] - p['E_over_Ec'], R_Y_ratio=s['R_Y_rel'] / p['R_Y_rel']))
H0 = 1.0764e11
for it in range(6):          # iterate L(H) with the local derivative
    L = math.log(2.435e18 / (30 * H0))
    dl = (L - p['L']) * out['selfconsistent_vs_1331']['dlnNX_per_unit']
    H0 = 1.0764e11 * math.exp(-(2 / 3) * dl)
out['H_freeX_planck_selfconsistent'] = H0
out['L_selfconsistent'] = math.log(2.435e18 / (30 * H0))
out['dq_entropy_claim'] = 'Self-consistent H_i would give 13.49, about 0.2% in H_i'
json.dump(out, open(sys.argv[2], 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'rows'}, indent=1))
for k, v in rows.items(): print(k, {kk: (round(vv, 6) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != 'END'})
