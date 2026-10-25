#!/usr/bin/env python3
"""Controls for equilibrium.py (case H, sudden energy): step doubling, 4->2 coefficients x10 and x0.1."""
import json
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import equilibrium as eq

OUT = Path(__file__).resolve().parents[1] / 'results'


def job(args):
    import warnings; warnings.filterwarnings('ignore')
    eps, tag = args
    A0 = dict(eq.A)
    if tag == 'A_x10':
        for k in eq.A: eq.A[k] = 10 * A0[k]
    if tag == 'A_x0.1':
        for k in eq.A: eq.A[k] = 0.1 * A0[k]
    if tag == 'steps_x2':

        orig = np.geomspace
        eq.np.geomspace = lambda a, b, n: orig(a, b, 2 * n if n == 240 else n)
    r = eq.run(eps, 'H', 'sudden')
    return dict(eps=eps, control=tag, p_final=r.get('p_final'), f_extra=r.get('f_extra_pairs'), T_vis_freeze=r.get('T_vis_freeze'))


if __name__ == '__main__':
    jobs = [(e, t) for e in (0.3, 0.2, 0.1, 0.03) for t in ('baseline', 'steps_x2', 'A_x10', 'A_x0.1')]
    with Pool(3, maxtasksperchild=1) as p:  # fresh interpreter state per control
        out = p.map(job, jobs, chunksize=1)
    for r in out: print(r)
    (OUT / 'equilibrium_controls.json').write_text(json.dumps(out, indent=1, default=float))
