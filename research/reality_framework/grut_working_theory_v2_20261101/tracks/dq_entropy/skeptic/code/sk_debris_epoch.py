#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code): the free-Y and sigma quanta at the epochs where they decay (sigma: Gamma_sigma = H
at a ~ 2.2; Y_H: <a> = 9.4), against the late values (a >= 315) that dq_entropy's ledger uses.  rpdq.c (track source,
read-only build) to H_i t = 60 with dense samples; physical-time running means over one Y/sigma oscillation are not needed
because N_YH, N_YL, N_sigma are quadratic adiabatic numbers already (rpdq diagnostic).  Usage: sk_debris_epoch.py OUT.json"""
import json, os, sys
from pathlib import Path
import numpy as np
os.environ.setdefault('DQ_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_entropy_skeptic')
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'code'))
import dq_common as dq
p, st, G, du = dq.prepare(0.5, 3.0, 0.0)
rr = dq.run('debris60', p, 3.0, 0.5, st['z'], G=G, cport=0.0, end=60.0, rtol=1e-10, atol=1e-13, nsamp=6000, binw=0.5, cbin=0.0, save=False)
S = rr['S']; C = dq.COLS
t, a = S[:, 0], S[:, 1]
def at(name, aval):
    i = int(np.argmin(abs(a - aval))); j = slice(max(i - 50, 0), i + 50)
    return float(np.mean(S[j, C.index(name)]))
rows = []
for av in (1.5, 2.2, 3.0, 5.0, 7.5, 9.4, 12.0, 15.0):
    rows.append(dict(a=av, N_sigma=at('Nsig', av), N_YH=at('NYH', av), N_YL=at('NYL', av)))
    print(rows[-1])
late = dict(N_sigma=1.9879757446145123e-05, N_YH=0.010303911565136449, N_YL=0.5664202557423563)
json.dump(dict(rows=rows, late_dq=late), open(sys.argv[1], 'w'), indent=1)
