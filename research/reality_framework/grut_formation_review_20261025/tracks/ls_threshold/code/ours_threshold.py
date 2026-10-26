#!/usr/bin/env python3
"""ls_threshold, reproduction with OUR earlier code (not Codex's): case-L threshold R_L by bisection.

Imports run_inj from tracks/sp_yledger/code/thermal_relic.py UNCHANGED (which uses gp_leftovers/equilibrium.py:
adaptive-quad Bose integrals, 260 log steps to T_vis = 1e-6 GeV, explicit endpoint exponential update,
review g*s table).  Nothing is written into those tracks (no bytecode).
Inputs [Q from b448f25 DECAY_RESULTS.json]: epsilon_net and raw ymult of candidates 10 and 17.
Usage: python3 -B ours_threshold.py <out.json>     (~40-60 min on 4 cores)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'sp_yledger' / 'code'))
import thermal_relic as TR  # noqa: E402

CANDS = {10: dict(eps=0.15741283451037602, ymult=5.69009662564405, R_Y=113.0689908402781, R_lower=199.4042010902021),
         17: dict(eps=0.1470461788202921, ymult=5.309416605350421, R_Y=122.80645189533328, R_lower=249.69026892092234)}


def bisect(job):
    cand, ym = job
    c = CANDS[cand]; lo, hi = 60.0, 140.0; hist = []
    for _ in range(7):
        mid = 0.5 * (lo + hi)
        r = TR.run_inj(c['eps'], 'L', mid, ym)
        hist.append(dict(R=mid, f=r['f_extra'], sigma=r['sigma'], max_cond=r['max_cond_frac']))
        if r['f_extra'] > 0.01: hi = mid
        else: lo = mid
    return dict(candidate=cand, ymult=ym, R_L_bracket=[lo, hi], history=hist)


def point(job):
    cand, name, ym = job
    c = CANDS[cand]; r = TR.run_inj(c['eps'], 'L', c[name], ym)
    return dict(candidate=cand, ledger=name, ymult=ym, **{k: r[k] for k in ('R', 'sigma', 'T_h0_over_m', 'p_final', 'f_extra', 'max_cond_frac', 'T_vis_freeze')})


if __name__ == '__main__':
    jobs = [(k, ym) for k in CANDS for ym in (1.0, CANDS[k]['ymult'])]
    pts = [(k, n, CANDS[k]['ymult']) for k in CANDS for n in ('R_Y', 'R_lower')]
    with Pool(4) as p:
        th = p.map(bisect, jobs)
        Path(sys.argv[1]).write_text(json.dumps(dict(thresholds=th), indent=1, default=float))
        pr = p.map(point, pts)
    out = dict(method=__doc__, thresholds=th, points=pr)
    Path(sys.argv[1]).write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))
