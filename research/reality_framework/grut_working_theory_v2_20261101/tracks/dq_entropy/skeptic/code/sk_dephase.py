#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code): a computed toy of dq_entropy's fragmentation-dephasing systematic [E].
dq_entropy argues that if fragmentation (H_i t ~ 60-70, M1 value) dephases the X <-> Sigma exchange, N_X keeps the dressing
it had then, and estimates -0.5..+0.2% (4-orbit running means), up to -1.5..+0.7% (2-orbit).  Toy: from the rpdq state at
H_i t_f (sk_state1000.py), switch the only X-number-violating coupling kappa off with a cos^2 ramp over deta_r (sudden,
1, 2, 4, 8 X orbits; one orbit ~ 1.7 in eta) in sk_late.c, and record the frozen N_X relative to the full-coupling frozen
value -8047.8259 (sk_late long run).  Usage: sk_dephase.py OUT.json T_f [T_f ...]"""
import json, math, os, subprocess, sys
from pathlib import Path
import numpy as np
SP = Path(os.environ.get('SK_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_entropy_skeptic'))
HERE = Path(__file__).resolve().parent
BIN = SP / 'build' / 'sk_late'
NX_FROZEN = -8047.825918932946
ORBIT = 1.7
out = dict(method=__doc__, NX_frozen_full=NX_FROZEN, rows=[])
for tf in sys.argv[2:]:
    tfv = float(tf)
    stf = SP / 'runs' / ('state_t%g.json' % tfv)
    if not stf.exists():
        subprocess.run([sys.executable, '-B', str(HERE / 'sk_state1000.py'), 't%g' % tfv, '0', '%g' % tfv, str(stf)], check=True,
                       env=dict(os.environ, DQ_SCRATCH=str(SP), OMP_NUM_THREADS='1'))
    st = json.load(open(stf))
    Y = st['Y']; p = st['params']; last = st['last']
    a0 = math.exp(Y[0]); K = a0 * a0 * last['H']; eta0 = st['eta_end']
    line = [a0, eta0, K, p['lx'], p['ly'], p['ls'], p['kap'], p['M'], p['h'], p['mu']] + Y[3:9] + Y[9:15]
    inp = ' '.join('%.17g' % v for v in line)
    for norb in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        dr = max(norb * ORBIT, 1e-9)
        eend = eta0 + dr + 1.0
        o = subprocess.run([str(BIN), '%.17g' % eend, '0.2', '1', '0.25', '%.17g' % eta0, '%.17g' % dr], input=inp, capture_output=True,
                           text=True, check=True).stdout
        F = [list(map(float, l.split()[1:])) for l in o.split('\n') if l.startswith('F ')]
        # bins fully after the ramp: N_X is exactly conserved there (kappa = 0)
        post = [f for f in F if f[0] >= eta0 + dr - 1e-9]
        nx = post[-1][2]
        spread = max(f[2] for f in post) - min(f[2] for f in post) if len(post) > 1 else 0.0
        row = dict(t_f=tfv, a_f=a0, eta_f=eta0, NX_at_tf=last['NX'], ramp_orbits=norb, NX_final=nx, rel=nx / NX_FROZEN - 1, post_spread=spread)
        out['rows'].append(row)
        print('t_f=%5.1f a=%6.2f ramp=%4.1f orbits: N_X final = %.3f  rel = %+.4f%%  (inst N_X(t_f) rel %+.4f%%; post-ramp spread %.1e)' % (
            tfv, a0, norb, nx, 100 * row['rel'], 100 * (last['NX'] / NX_FROZEN - 1), spread), flush=True)
rel = np.array([[r['rel'] for r in out['rows'] if r['ramp_orbits'] == n] for n in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)])
out['range_by_ramp'] = {str(n): [float(rel[i].min()), float(rel[i].max())] for i, n in enumerate((0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0))}
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print(json.dumps(out['range_by_ramp'], indent=1))
