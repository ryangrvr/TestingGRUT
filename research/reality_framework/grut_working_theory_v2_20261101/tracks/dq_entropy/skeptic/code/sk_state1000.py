#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code): run dq_entropy's rpdq.c (built from the track's source into the skeptic scratch)
to H_i t = T_END for a given portal coefficient and save the full final state (the 'Y' line), the conformal time and the
400-1000 window, as the hand-off to the independent late-time integrator sk_late.c, and for the portal-coefficient
checks.  Usage: sk_state1000.py TAG CPORT T_END OUT.json
dq_common.py is imported read-only with DQ_SCRATCH pointing at the skeptic scratch."""
import json, os, sys
from pathlib import Path
os.environ.setdefault('DQ_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_entropy_skeptic')
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'code'))
import dq_common as dq
import numpy as np

tag, cport, end = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
p, st, G, du = dq.prepare(0.5, 3.0, cport)
rr = dq.run(tag, p, 3.0, 0.5, st['z'], G=G, cport=cport, end=end, rtol=1e-10, atol=1e-13, nsamp=int(end / 5), binw=20.0, cbin=0.25, save=False)
o = dict(tag=tag, cport=cport, end=end, Y=rr['Y'], eta_end=rr['eta_end'], END=rr['END'], params=p, G=[G.real, G.imag], prep=[[c.real, c.imag] for c in st['z']])
if end >= 1000:
    for k in ('NX', 'EX', 'NYL', 'NYH', 'Nsig', 'NS'):
        o[k + '_400_1000'] = dq.window_mean(rr['bins'], k, 400.0, 1000.0)
    o['yield_over_YQ_400_1000'] = dq.fc.yield_over_target(dict(S=rr['S']), abs(o['NX_400_1000']))[1]
last = dict(zip(dq.COLS, rr['S'][-1].tolist()))
o['last'] = last
cb = rr['cbins']
o['cbins_tail'] = cb[-400:].tolist()
json.dump(o, open(sys.argv[4], 'w'), indent=1, default=float)
print(tag, 'NX(400-1000)=%.9f' % o.get('NX_400_1000', float('nan')), 'eta_end', o['eta_end'], 'a', last['a'], 'H', last['H'], o['END'])
