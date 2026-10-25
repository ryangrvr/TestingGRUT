#!/usr/bin/env python3
"""sp_yledger (Claude, own code): T_R dependence of the inflaton-matter release at the leading point,
delta = pi/2, Hubble mass removed at a = 1 (bg 1).  Usage: python3 md_leading_TR.py <out.json>  (~2 min)"""
import json, sys
from multiprocessing import Pool
import numpy as np
from md_release import one, HLEAD, TR_STATED
from ymodel import minimum, entropy_inst

if __name__ == '__main__':
    m = minimum(HLEAD, 1.0); amp = (m['x'], m['y'], m['s'])
    TRs = [1e13, TR_STATED, 1e14, 3e14, 1e15]
    with Pool(4) as p:
        rows = p.map(one, [(HLEAD, np.pi / 2, TR, 1, 1.5e6, amp) for TR in TRs])
    json.dump(dict(T_RI_GeV=entropy_inst(HLEAD)[1], rows=rows), open(sys.argv[1], 'w'), indent=1, default=float)
    for r in rows:
        w = r['600-1000']
        print(f"TR={r['T_R_GeV']:.3g} TR/TRI={r['T_R_over_instant']:.4f} eps={w['eps']:+.4f} NX={w['NX']:.4g} Y/YQ={w['YX_over_YQ']:+.4g}")
