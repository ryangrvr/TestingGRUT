"""Verifier v_reach: recompute Codex 479f111 production component statistics (read-only extract).
Usage: python3 -I v_codex479.py <dir with production_*.csv.components.csv> <out.json>
Qbar units (lattice), physical Q = 100 Qbar (lam_phi = 0.01). Box 64^3 ell^3 (L = 32 half-width, h = 0.5, N = 128).
"""
import sys, csv, json
import numpy as np
d, out = sys.argv[1], sys.argv[2]
QB_PHYS = 1.77e4      # reach claim; verifier's own value inserted below if supplied
res = {}
tot = 0.2200852108 * 64**3
for grid in ('coarse', 'medium', 'fine'):
    rows = list(csv.DictReader(open(f'{d}/production_{grid}.csv.components.csv')))
    for th in sorted(set(r['threshold'] for r in rows), key=float):
        Q = np.array([float(r['Q']) for r in rows if r['threshold'] == th])
        Q = np.sort(Q)[::-1]
        above = Q[Q * 100 > QB_PHYS]
        res[f'{grid}_{th}'] = dict(n=int(len(Q)), charge_fraction=float(Q.sum() / tot), top5_Qbar=[float(x) for x in Q[:5]],
                                   cw_mean_Qbar=float((Q**2).sum() / Q.sum()), n_above_Qb177=int(len(above)),
                                   frac_of_component_charge_above_Qb=float(above.sum() / Q.sum()),
                                   median_Qbar=float(np.median(Q)))
json.dump(dict(box_Qbar_total=tot, stats=res), open(out, 'w'), indent=1)
for k, v in res.items():
    print(k, v['n'], round(v['charge_fraction'], 3), [round(x, 1) for x in v['top5_Qbar']], round(v['cw_mean_Qbar'], 1), v['n_above_Qb177'], round(v['frac_of_component_charge_above_Qb'], 3))
