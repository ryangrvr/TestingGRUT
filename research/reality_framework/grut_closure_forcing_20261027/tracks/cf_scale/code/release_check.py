#!/usr/bin/env python3
"""cf_scale P4b consequence: does the heavy-mediator spot point survive with lambda_Y = lambda_Sigma lowered into the
P4b-compatible band (perturbative to M_P and kappa below its boundedness bound)? Own driver around the hw_gates
spot model (research/reality_framework/grut_formation_review_20261025/tracks/hw_gates/code/spotmodel.py and hwgrel.c,
used read-only; the binary is built in this track's scratch). Same pipeline as hw_gates gate 2: Y-relaxed preparation
at delta = pi/2, instant release, window a = 30-60, ledger and case-L proxy R_L(eps) from the sp_critic table (which
hw_critic found 7-8% high: R_L(spot) = 81.5-82, not 87.9). Every other spot coupling is unchanged.
Usage: python3 -B release_check.py <out.json>
"""
import json, os, sys, time
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
os.environ.setdefault('HWG_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                                     'scratchpad/closure/cf_scale/hwg')
HWG = Path(__file__).resolve().parents[4] / 'grut_formation_review_20261025' / 'tracks' / 'hw_gates' / 'code'
sys.path.insert(0, str(HWG))
import spotmodel as S


def one(lam):
    P = S.params(ly=lam, ls=lam)
    z0, info = S.prep('relaxed', np.pi / 2, P)
    t0 = time.time()
    r = S.run(np.array(z0, complex), P, a_end=60.0, rtol=1e-10, atol=1e-13, nsamp=30, win0=30.0)
    w = r['window']; L = S.ledger(P, w); RL = S.R_L(w['eps'])
    x, y, s = S.aligned(P)
    return dict(lam_Y=lam, lam_S=lam, aligned_inflation_minimum_H=[x, y, s], prep=dict(x=info['x'], Y=info['Y'], s=info['s']),
                eps=w['eps'], NX=w['NX'], yield_over_YQ=S.yield_over_YQ(w['NX']), NYH=w['NYH'], NYL=w['NYL'], Nsig=w['Nsig'],
                R_central=L['R_central'], R_L_table=RL, R_L_corrected_minus7pct=None if RL is None else 0.93 * RL,
                pass_factor=None if RL is None else RL / L['R_central'], max_gauss=r['end']['max_gauss'],
                max_ledger=r['end']['max_ledger'], runtime_s=time.time() - t0)


def main(out):
    rows = []
    for lam in [float(v) for v in (sys.argv[2].split(',') if len(sys.argv) > 2 else ['1.0', '0.566', '0.45', '0.37', '0.3', '0.25'])]:
        row = one(lam); rows.append(row)
        print('lam_Y=lam_S=%.3f  eps=%.5f  NX=%.1f  yield=%.3f Y_Q  R=%.3f  R_L=%s  pass x%s  (%.0fs)' % (
            lam, row['eps'], row['NX'], row['yield_over_YQ'], row['R_central'], row['R_L_table'], row['pass_factor'], row['runtime_s']), flush=True)
    Path(out).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))


if __name__ == '__main__':
    main(sys.argv[1])
