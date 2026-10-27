"""Part 2e: sensitivity of the transition outcome to a 1e-12 relative change of the phi seed (chaos test),
and the outcome distribution versus initial orbital phase at fixed seed.  Independent code."""
import sys, os, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *
from part2_floquet import orbit_info
from part2_transition import transition_run, phase_state

def job(a):
    eps, ph, rel = a
    y4 = tau20_state(); info = orbit_info(y4); st = phase_state(y4, info['T_r'], ph)
    r = transition_run(st, info, eps, 1e-6 * (1 + rel), keep_trace=False)
    return dict(eps=eps, phase=ph, rel=rel, nminus=r['nminus'], Nphi=r['Nphi'])

if __name__ == '__main__':
    jobs = [(eps, ph, rel) for eps in [1e-2, 1e-3, 1e-4] for ph in [0.0, 0.25, 0.5, 0.75] for rel in [0.0, 1e-12]]
    with Pool(2) as p: rows = p.map(job, jobs)
    out = []
    for eps in [1e-2, 1e-3, 1e-4]:
        for ph in [0.0, 0.25, 0.5, 0.75]:
            a = [r for r in rows if r['eps'] == eps and r['phase'] == ph]
            d = dict(eps=eps, phase=ph, nminus=[x['nminus'] for x in a], Nphi=[x['Nphi'] for x in a], dn=abs(a[0]['nminus'] - a[1]['nminus']))
            out.append(d); print(d)
    save_json('transition_chaos.json', dict(rows=out, note='pairs differ only by a 1e-12 relative change of the phi seed'))
