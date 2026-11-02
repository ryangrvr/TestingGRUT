#!/usr/bin/env python3
"""Track mf_higgs: assemble results/summary.json from the result files (and two cached mode runs for the band)."""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np

import mh_common as C
import mh_run as RN

R = C.RES
L = lambda f: json.load(open(R / f))
post = L('post_D.json')
infl = L('inflation.json')
sens = L('inflation_sens.json')
pl = L('plasma.json')
th = L('thermal_reach.json')
al = L('alphas.json')
out = dict(label='C on the stated model choices; criterion D and the hard-gas plasma are [E]',
           inputs=dict(PDG2024=C.PDG24, xi_EW=0.1, note='xi_H = 0.1 read as xi(M_t); xi(H_i) = 0.0877'))

# post-inflation, per history
hist = {}
for k, o in post.items():
    ext = o['Mt_ext']
    hist[k] = dict(K=o['K'], N_obs=o['N_obs'],
                   D_grav_PDG={m: ext[m]['grav']['D'] for m in ('172.28', '172.57', '172.86')},
                   D_grav_N60_PDG={m: ext[m]['grav_D60'] for m in ('172.28', '172.57', '172.86')},
                   D_plasma_ext={m: ext[m]['plasma']['D'] for m in ext},
                   plasma_note={m: ext[m]['plasma'].get('note') for m in ext},
                   plasma_x_at_best={m: ext[m]['plasma'].get('x') for m in ext},
                   Mt_crit_grav=o['Mt_crit_grav'], Mt_crit_grav_band=o['Mt_crit_grav_band'],
                   Mt_crit_grav_N60=o['Mt_crit_grav_N60'],
                   Mt_crit_variants={kk: v['Mt_crit'] for kk, v in o['variants'].items()},
                   deadlines={m: v['deadline'] for m, v in o['deadlines'].items()},
                   Mt_crit_plasma=o['Mt_crit_plasma'], Mt_crit_plasma_c1e_4=o['Mt_crit_plasma_c1e-4'])
out['post_inflation'] = hist

# criterion band at exactly 172.86 (EXP, PAULI0)
band = {}
for h in ('EXP', 'PAULI0'):
    bg = RN.bgc('HFa', h)
    Rr = RN.run_cached('HFa', h, RN.xi_at(bg.Hi_GeV))
    N = Rr['res']['N_obs']
    vals = [RN.D_of(Rr, 172.86, N, cH=c)['D'] for c in (0.5, 2.0)] + [RN.D_of(Rr, 172.86, N, kc_f=f)['D'] for f in (0.5, 2.0)]
    band[h] = dict(D=RN.D_of(Rr, 172.86, N)['D'], band=[min(vals), max(vals)])
out['criterion_band_172.86'] = band

# inflation
inf = {}
for name in ('HFa', 'HFb'):
    o = infl[name]
    inf[name] = dict(N0=o['N0'], caps=o['Mt_cap_xi0.1'],
                     xi_min={m: {kk: v for kk, v in r.items() if 'xi_min' in kk} for m, r in o['rows'].items()},
                     log10N={m: {kk: v for kk, v in r.items() if 'log10N' in kk} for m, r in o['rows'].items()})
out['inflation'] = inf
out['inflation_sens'] = sens
out['inflation_tail'] = {n: {m: v['best_log10_patches'] for m, v in d.items()} for n, d in L('inflation_tail.json').items()}
out['alphas'] = al

# plasma
out['plasma'] = {n: {h: dict(K=v['K'], t_fSM=v['t_fSM'], t_m2ind=v['t_m2ind'], m2ind_at=v['m2ind_at'],
                             t_thermal_pessimistic=v['t_thermal_pessimistic']) for h, v in d.items()} for n, d in pl.items()}
out['thermal_only_reach'] = {h: v['Mt_reach_thermal_only'] for h, v in th.items()}

# xi window
edges = {}
for h in ('EXP', 'PAULI0'):
    e = L(f'xiscan_HFa_{h}.json')['edges_xiEW']
    edges[h] = e
up = {}
for h in ('EXP', 'PAULI0'):
    u = L(f'upper_HFa_{h}.json')
    xs = sorted(float(x) for x in u)
    res = {}
    for Mt in ('172.28', '172.57', '172.86'):
        r = [u[f'{x}'][Mt]['two_sigma_prod_over_hmax'] for x in xs]
        cross = None
        for i in range(len(xs) - 1):
            if r[i] < 1 <= r[i + 1]:
                cross = math.exp(math.log(xs[i]) + (math.log(1 / r[i])) * (math.log(xs[i + 1]) - math.log(xs[i])) /
                                 (math.log(r[i + 1]) - math.log(r[i])))
        res[Mt] = dict(xi_H_meanfield_edge=cross,
                       xi_EW_meanfield_edge=(RN.smt(172.57).xi_EW(1.30002e11, cross) if cross else None),
                       max_D_grav=max(u[f'{x}'][Mt]['grav'] for x in xs),
                       any_runaway=any(u[f'{x}']['runaway'] for x in xs))
    up[h] = res
out['xi_window'] = dict(lower_post_inflation_xiEW=edges, upper_gravitational=up,
                        upper_with_plasma_barrier=L('upper_plasma.json'))

# portals and induced terms
out['g6_portal'] = {h: {k: dict(Lam=v['Lam'], lam_chiH=v['lam_chiH'], runaway=v['runaway'],
                                D_grav={m: v[m]['D_grav'] for m in ('172.57', '172.86', '173.45', '174.05', '174.6')},
                                inflation_log10N={m: v[m]['inflation_log10N_tree'] for m in ('172.86', '174.6')})
                        for k, v in L(f'g6_portal_{h}.json').items() if k != 'conversion'} for h in ('EXP', 'PAULI0')}
out['induced'] = L('induced.json')
out['windows'] = {h: {k: v for k, v in w.items() if k != 'first'} | {'first_two': w['first'][:2]}
                  for h, w in L('windows.json').items()}
out['controls'] = L('controls.json')
out['H_i'] = {h: {k: dict(Mt_crit_grav=v['Mt_crit_grav'], Mt_crit_plasma=v['Mt_crit_plasma'], tau_m2ind_1=v['tau_m2ind_1'])
                  for k, v in L(f'Hi_HFa_{h}.json').items()} for h in ('EXP', 'PAULI0')}
C.jdump(out, R / 'summary.json')
print(json.dumps(dict(band=band, upper=up, edges=edges, thermal=out['thermal_only_reach'], Hi=out['H_i']), indent=1,
                 default=float))
