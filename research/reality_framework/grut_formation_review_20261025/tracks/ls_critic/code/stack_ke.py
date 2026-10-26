#!/usr/bin/env python3
"""ls_critic (Claude): the one corner where stacking favourable choices passes the benchmark proxy (candidate 10,
Y-relaxed preparation, conversion sigma v x72, on-shell s route, zero-cascade ledger R = 148.86) re-run in the
energy-conserving ke mode of tracks/ls_threshold/code/lsrelic.py (unchanged), mu_c = E_c = 0.861195 and 0.923316.
Usage: OMP_NUM_THREADS=1 python3 -B stack_ke.py <out.json>     (~5 min)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'ls_threshold' / 'code'))
import lsrelic as LR  # noqa: E402
from stack_check import REL  # noqa: E402


def job(a):
    uc, sv = a
    LR.UC = uc; LR.EC = uc; LR.SV *= sv
    r = LR.ke_run(REL[10]['eps'], REL[10]['R_lower'], REL[10]['ym'], N=600, antilumps=True, s_channel=True)
    r.pop('history')
    r.update(kind='ke_stack', candidate=10, mu_c=uc, sv_factor=sv)
    return r


if __name__ == '__main__':
    with Pool(2, maxtasksperchild=1) as p:
        rows = p.map(job, [(0.8611951308217225, 72.0), (0.9233161353726057, 72.0)])
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        print('ke stack mu_c=%.4f f=%.4g fallbacks=%d' % (r['mu_c'], r['f'], r['root_fallbacks']))
