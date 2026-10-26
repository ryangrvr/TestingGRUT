#!/usr/bin/env python3
"""hw_verify (Claude, own code): derived headline numbers from the result files (reproducible digest input).
Usage: python3 -B summary.py ../results     (reads ladder/phases/scaling/prep/release_*/rl_actual*.json, writes summary.json)
"""
import json, sys
from pathlib import Path
import numpy as np
R = Path(sys.argv[1])
out = {}

def cross(pts, target=0.01):
    """log-linear crossing of f = target in R from (R, f) points."""
    pts = sorted(pts)
    for (r0, f0), (r1, f1) in zip(pts, pts[1:]):
        if f0 <= target < f1:
            t = (np.log(target) - np.log(f0)) / (np.log(f1) - np.log(f0))
            return float(r0 + t * (r1 - r0)), [r0, r1]
    return None, None

# direct case-L threshold at the actual eps
rows = []
for fn in ('rl_actual.json', 'rl_actual_fine.json'):
    rows += [r for r in json.load(open(R / fn))['rows'] if 'f' in r]
rl = {}
for ym in sorted({r['ymult'] for r in rows}):
    pts = [(r['R'], r['f']) for r in rows if r['ymult'] == ym]
    v, br = cross(pts)
    rl['ymult_%g' % ym] = dict(R_L=v, bracket=br, points=sorted(pts))
out['R_L_direct_eps_0.1457'] = rl
lad = json.load(open(R / 'ladder.json'))['rows']
spot = next(r for r in lad if r['name'] == 'm30' and abs(r['delta'] - np.pi / 2) < 1e-6)
out['R_L_table_at_spot'] = spot['R_L']
out['R_L_table_over_direct'] = spot['R_L'] / rl['ymult_1']['R_L']

# ladder table and the debris edge
tab = []
for r in sorted([r for r in lad if r['name'].startswith('m') and abs(r['delta'] - np.pi / 2) < 1e-6], key=lambda r: r['couplings_over_H']['M']):
    tab.append(dict(m=r['couplings_over_H']['M'], eps=r['eps'], R_Y=r['R_Y'], R_Sigma=r['R_Sigma'], R_cascade=r['R_cascade'],
                    R_central=r['R_central'], R_L_table=r['R_L'], ratio=r['R_central'] / r['R_L'], yield_over_YQ=r['yield_over_YQ'],
                    c6=r['c6'], NYH=r['NYH'], NYL=r['NYL'], Nsig=r['Nsig']))
out['ladder_pi_over_2'] = tab
a, b = tab[0], tab[1]
slope = np.log(a['ratio'] / b['ratio']) / np.log(b['m'] / a['m'])
out['debris_edge_m'] = dict(table_RL=float(a['m'] * a['ratio']**(1 / slope)), slope=float(slope),
                            corrected_RL=float(a['m'] * (a['ratio'] / out['R_L_table_over_direct'])**(1 / slope)))
# allowances at the spot point (direct R_L)
RLd = rl['ymult_1']['R_L']
per_light = spot['deposit_per_Y'] / spot['EX']; per_heavy = per_light + spot['deposit_per_A'] / spot['EX']
out['allowance_at_spot'] = dict(R_per_light_Y=per_light, R_per_heavy_Y=per_heavy, R_per_vector=spot['deposit_per_A'] / spot['EX'],
                                heavy_quanta_allowed=(RLd - spot['R_Y'] - spot['R_Sigma']) / per_heavy,
                                light_quanta_allowed=(RLd - spot['R_cascade'] - spot['R_Sigma']) / per_light,
                                per_unit_charge_heavy=(RLd - spot['R_Y'] - spot['R_Sigma']) / per_heavy / abs(spot['NX']),
                                per_unit_charge_light=(RLd - spot['R_cascade'] - spot['R_Sigma']) / per_light / abs(spot['NX']))

# phases and edge fraction
ph = json.load(open(R / 'phases.json'))['rows']
rel = sorted([r for r in ph if r['name'] == 'phase'], key=lambda r: r['delta'])
out['phase_relaxed'] = dict(R_range=[min(r['R_central'] for r in rel), max(r['R_central'] for r in rel)],
                            best=max(rel, key=lambda r: abs(r['eps']))['delta'], best_eps=max(rel, key=lambda r: abs(r['eps']))['eps'],
                            band_abs_eps_ge_0p10_deg=[r['delta'] * 180 / np.pi for r in rel if abs(r['eps']) >= 0.10][::len([r for r in rel if abs(r['eps']) >= 0.10]) - 1],
                            band_abs_eps_ge_0p13_deg=[r['delta'] * 180 / np.pi for r in rel if abs(r['eps']) >= 0.13][::len([r for r in rel if abs(r['eps']) >= 0.13]) - 1],
                            max_ratio_where_defined=max(r['R_central'] / r['R_L'] for r in rel if r['R_L']))
cod = sorted([r for r in ph if r['name'] == 'phase_codex'], key=lambda r: r['delta'])
out['phase_codex'] = dict(R_range=[min(r['R_central'] for r in cod), max(r['R_central'] for r in cod)],
                          min_ratio_where_defined=min(r['R_central'] / r['R_L'] for r in cod if r['R_L']))
fr = sorted([r for r in ph if r['name'] == 'frac' and 'eps' in r], key=lambda r: r['frac'])
out['edge_fraction'] = [dict(frac=r['frac'], eps=r['eps'], R=r['R_central'], R_L=r['R_L'], yield_over_YQ=r['yield_over_YQ']) for r in fr]
out['edge_fraction_failed'] = [dict(frac=r['frac'], status=r['status']) for r in ph if 'status' in r]

# scaling
sc = json.load(open(R / 'scaling.json'))['checks']
out['scaling'] = dict(max_abs_d_eps=max(abs(c['d_eps']) for c in sc), max_rel_dR=max(abs(c['d_Rcentral_rel']) for c in sc),
                      max_rel_yield_dev=max(abs(c['yield_ratio'] / c['s_pow_1p5'] - 1) for c in sc),
                      max_rel_c6_dev=max(abs(c['c6_ratio'] / c['s_pow_m2'] - 1) for c in sc), s_values=[c['s'] for c in sc])

# preparation
pr = json.load(open(R / 'prep.json'))['results']
held = [(float(r['tag'].split('_N')[1]), r) for r in pr['releases'] if r['tag'].startswith('codex_held')]
out['prep'] = dict(modes_omega_over_H=[m['omega_over_H'] for m in pr['linearised']['held_modes']],
                   amplitude_decay_rate_over_H=1.5, NY_decay_fit=pr['NY_decay_fit'],
                   dY_codex=pr['displacement_codex_vs_relaxed']['ymodes_codex']['abs_dY'],
                   dY_relaxed=pr['displacement_codex_vs_relaxed']['ymodes_relaxed']['abs_dY'],
                   codex_held=[dict(N=N, eps=r['eps'], NY=r['NYH'] + r['NYL'], R=r['R_central'], R_L=r['R_L'], dY=r['ymodes_at_release']['abs_dY']) for N, r in held],
                   interp=[dict(f=float(r['tag'].split('_f')[1]), eps=r['eps'], NYH=r['NYH'], R=r['R_central'], R_L=r['R_L']) for r in pr['releases'] if r['tag'].startswith('interp')],
                   free=[dict(tag=r['tag'], eps=r['eps'], R=r['R_central'], argX=r['phases_at_release']['argX'], dS=r['phases_at_release']['argS_minus_pi'] - np.pi / 2) for r in pr['releases'] if 'free' in r['tag']],
                   phase_gradient=pr['linearised']['phase_gradient'], phase_mass2=pr['linearised']['phase_mass2_over_H2'])
# N at which R(held) first drops below R_L (interpolated in log R)
pts = [(N, r['R_central'], r['R_L']) for N, r in held if r['R_L']]
for (N0, R0, L0), (N1, R1, L1) in zip(pts, pts[1:]):
    if R0 > L0 and R1 <= L1:
        t = np.log(R0 / L0) / (np.log(R0 / L0) - np.log(R1 / L1)); out['prep']['N_cross_R_L'] = float(N0 + t * (N1 - N0)); break
ip = sorted([(d['f'], d['R'], d['R_L']) for d in out['prep']['interp']])
for (f0, R0, L0), (f1, R1, L1) in zip(ip, ip[1:]):
    if R0 > L0 and R1 <= L1:
        t = np.log(R0 / L0) / (np.log(R0 / L0) - np.log(R1 / L1)); out['prep']['f_cross_R_L'] = float(f0 + t * (f1 - f0)); break

# finite release
rel = {}
for g in ('heavy', 'stated_codex', 'stated_relaxed'):
    rows = json.load(open(R / ('release_%s.json' % g)))['rows']
    rows.sort(key=lambda r: -(r['gamma_over_H'] or 1e9))
    rel[g] = [dict(gamma_over_H=r['gamma_over_H'], eps=r['eps'], NX=r['NX'], NY=r['NYH'] + r['NYL'], Nsig=r['Nsig'], R=r['R_central'],
                   R_reheat=r['R_central_reheat'], R_L_table=r['R_L'], yield_over_YQ=r['yield_over_YQ'], yield_diluted=r['yield_diluted_over_YQ'],
                   Rr_end=r['Rr_end']) for r in rows]
    base = rel[g][0]
    for d in rel[g]:
        d['eps_over_instant'] = d['eps'] / base['eps']; d['R_over_instant'] = d['R'] / base['R']
out['finite_release'] = rel
(R / 'summary.json').write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(out, indent=1, default=float))
