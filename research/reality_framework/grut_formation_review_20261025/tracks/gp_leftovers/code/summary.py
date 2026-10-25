#!/usr/bin/env python3
"""Collect the allowed-eps statements of track gp_leftovers into results/summary.json (own code)."""
import json
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parents[1] / 'results'
bud = json.loads((R / 'budget.json').read_text())
eq = json.loads((R / 'equilibrium.json').read_text())
cat = json.loads((R / 'catalysis.json').read_text())
ch = json.loads((R / 'channels.json').read_text())
vec = json.loads((R / 'vector_source.json').read_text())


def crossing(rows, level):
    pts = sorted([(r['eps'], r['f_extra_pairs']) for r in rows if 'f_extra_pairs' in r and r['f_extra_pairs'] > 0])
    for (e1, f1), (e2, f2) in zip(pts[:-1], pts[1:]):   # increasing eps, decreasing f
        if f1 >= level > f2:
            t = (np.log(level) - np.log(f1)) / (np.log(f2) - np.log(f1))
            return float(np.exp(np.log(e1) + t * (np.log(e2) - np.log(e1))))
    return None


H = [r for r in eq if r.get('case') == 'H' and r.get('energy_rule') == 'sudden']
L = [r for r in eq if r.get('case') == 'L' and r.get('energy_rule') == 'sudden']
out = dict(
    frozen_remainder=bud['remainder_laws']['thresholds'],
    max_equilibrium_homogeneous_condensate={str(l): crossing(H, l) for l in (0.01, 0.1, 1.0)},
    max_equilibrium_lumps_conversion=[dict(eps=r['eps'], p_final=r.get('p_final'), min_cond_frac=r.get('min_cond_frac'),
                                           final_cond_frac=r.get('final_cond_frac'), status=r.get('status')) for r in L],
    conversion_charge_budget=cat['thresholds'],
    no_free_gas_state_sudden_eps_ge=[r['eps'] for r in H if 'status' in r],
    Gamma_equals_H=ch['Gamma_equals_H_T_vis_GeV'],
    vector_source=[dict(H_inf=v['H_inf'], pairs_at_decay=v['pairs_per_net_charge_at_decay'], p_final=v['p_final'],
                        f_extra=v['f_extra']) for v in vec])
(R / 'summary.json').write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(out, indent=1, default=float))
