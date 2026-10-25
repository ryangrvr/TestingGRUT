#!/usr/bin/env python3
"""Track sp_critic (Claude, own driver; reuses our sp_yledger/gp_leftovers benchmark code).

Question: in the maximal-equilibrium benchmark, how much thermalised debris energy R = E_inj/E_X can a release
with charge fraction eps carry before the symmetric relic exceeds 1% of calibrated dark matter (f = 0.01)?
This turns the point-specific relic (sp_yledger) into a threshold R_max that can be compared with the Y excitation
of every release in Codex's published scan (critic_checks.py, part D).

Imports run_inj from tracks/sp_yledger/code/thermal_relic.py (which itself imports tracks/gp_leftovers/code).
No bytecode is written into those directories.
Usage: python3 -B rmax_relic.py <out.json> [eps,eps,...] [R,R,...] [H|L|HL]   (~10-20 min, 4 cores)
  default grid: eps 0.2402, 0.19; R 0-30; both cases.  Extension used for case L:
  python3 -B rmax_relic.py ../results/rmax_relic_L.json 0.240201575085664,0.19,0.3267 30,50,75,100,150,200 L
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.dont_write_bytecode = True
YL = Path(__file__).resolve().parents[2] / 'sp_yledger' / 'code'
sys.path.insert(0, str(YL))
import thermal_relic as TR  # noqa: E402  (our sp_yledger code)

EPS = [0.2402015750856640, 0.19]
RS = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]


def job(a):
    eps, case, R = a
    try:
        r = TR.run_inj(eps, case, R, 1.0)
        return {k: r[k] for k in ('eps', 'case', 'R', 'sigma', 'p_final', 'f_extra', 'cond_frac0',
                                  'max_cond_frac', 'final_cond_frac', 'T_h0_over_m')}
    except Exception as exc:  # solver failure is recorded, not hidden
        return dict(eps=eps, case=case, R=R, status='solver failure %r' % exc)


def rmax(rows, eps, case):
    """log-linear interpolation of the first crossing of f = 0.01 in R."""
    import math
    pts = sorted((r['R'], r['f_extra']) for r in rows
                 if r.get('eps') == eps and r.get('case') == case and 'f_extra' in r)
    for (r0, f0), (r1, f1) in zip(pts, pts[1:]):
        if f0 <= 0.01 < f1:
            if r0 == 0:
                return dict(bracket=[r0, r1], R_max=None)
            t = (math.log(0.01) - math.log(max(f0, 1e-300))) / (math.log(f1) - math.log(max(f0, 1e-300)))
            return dict(bracket=[r0, r1], R_max=math.exp(math.log(r0) + t * (math.log(r1) - math.log(r0))))
    return dict(bracket=None, R_max=None, note='no crossing in scanned R')


def main(out):
    global EPS, RS
    cases = ('H', 'L')
    if len(sys.argv) > 4:
        EPS = [float(x) for x in sys.argv[2].split(',')]
        RS = [float(x) for x in sys.argv[3].split(',')]
        cases = tuple(sys.argv[4])
    jobs = [(e, c, R) for e in EPS for c in cases for R in RS]
    with Pool(4) as p:
        rows = p.map(job, jobs)
    res = dict(method=__doc__, eps=EPS, R=RS, rows=rows,
               R_max={f'{e}_{c}': rmax(rows, e, c) for e in EPS for c in cases},
               scope='gp_leftovers maximal-equilibrium benchmark (full thermalisation, 4->2, case H homogeneous, '
                     'case L lumps + conversion). Benchmark, not kinetics. f = symmetric mass / calibrated DM.')
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    for r in rows:
        print({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()})
    print(json.dumps(res['R_max'], indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
