#!/usr/bin/env python3
"""ls_critic (Claude): does stacking favourable modelling choices inside the case-L benchmark rescue a candidate?

Each track probed one choice at a time.  Here the most favourable combination found across tracks is applied
together, in the benchmark's own isentropic proxy (not ke, which is stricter):
  * preparation: Y-relaxed (ls_scaling; re-checked in relaxed_check.py): cand 10 eps 0.16487, raw yield 5.970,
    R_lower 148.86, R_central 178.95; cand 17 eps 0.15783, yield 5.701, R_lower 192.69, R_central 224.61;
  * conversion sigma v x72 (ls_threshold's hot-daughter value), which raised R_L at Codex's preparation 103 -> 148;
  * the on-shell s route (ls_threshold), which lowers f.
Thresholds by bisection in proxy mode (lsrelic.threshold, N = 400) and direct f at the ledger values (N = 800).
Imports tracks/ls_threshold/code/lsrelic.py unchanged.
Usage: OMP_NUM_THREADS=1 python3 -B stack_check.py <out.json>     (~30 min on 4 cores)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'ls_threshold' / 'code'))
import lsrelic as LR  # noqa: E402

REL = {10: dict(eps=0.164868911425625, ym=5.970174242153348, R_lower=148.8550348931443, R_central=178.95271220937167),
       17: dict(eps=0.157831918500502, ym=5.701037145166892, R_lower=192.68519438524606, R_central=224.61344302089202)}


def thr(a):
    cand, sv, s_ch, lo, hi = a
    LR.SV *= sv
    r = LR.threshold(REL[cand]['eps'], REL[cand]['ym'], lo, hi, 9, N=400, s_channel=s_ch)
    r.update(kind='threshold', candidate=cand, sv_factor=sv, s_channel=s_ch)
    return r


def pt(a):
    cand, sv, s_ch, which = a
    LR.SV *= sv
    r = LR.proxy_run(REL[cand]['eps'], REL[cand][which], REL[cand]['ym'], N=800, s_channel=s_ch)
    r.update(kind='point', candidate=cand, sv_factor=sv, s_channel=s_ch, ledger=which)
    return r


if __name__ == '__main__':
    jobs = [(thr, (10, 1.0, False, 60.0, 220.0)), (thr, (10, 72.0, False, 60.0, 300.0)),
            (thr, (10, 72.0, True, 60.0, 300.0)), (thr, (17, 72.0, True, 60.0, 300.0)),
            (pt, (10, 72.0, True, 'R_lower')), (pt, (10, 72.0, True, 'R_central')),
            (pt, (17, 72.0, True, 'R_lower')), (pt, (10, 1.0, False, 'R_lower'))]
    rows = []
    with Pool(4, maxtasksperchild=1) as p:
        res = [(f.__name__, a, p.apply_async(f, (a,))) for f, a in jobs]
        for name, a, r in res:
            try:
                rows.append(r.get())
            except Exception as exc:
                rows.append(dict(kind='failed', job=name, args=list(a), status=repr(exc)))
            Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, inputs=REL, rows=rows), indent=1, default=float))
            x = rows[-1]
            print({k: x.get(k) for k in ('kind', 'candidate', 'sv_factor', 's_channel', 'ledger', 'R', 'f', 'bracket')},
                  flush=True)
