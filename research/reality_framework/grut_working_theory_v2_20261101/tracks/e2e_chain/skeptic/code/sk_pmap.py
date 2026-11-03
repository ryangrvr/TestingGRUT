#!/usr/bin/env python3
"""e2e_chain_skeptic (Claude): robustness of e2e_relic's gap fill of p(S/Q) to the relaxed-start interpolation.

e2e_relic's 'relaxed' rows in the gap (S/Q 7-21.8, R = 8-48) take x0 = ln(a_on/a_ref) LINEARLY in R between fx_chain's
heat_onset rows R = 5.703 (x0 = -0.6119) and 41.707 (-0.5565).  The computed x0(R) is concave (slope 3.6e-3 below
R = 5.7, 1.5e-3 across the gap, 1.0e-3 above 41.7), so the rule is a choice.  Here the same network (fx_chain's
m1_kin, READ-ONLY; work dir redirected to this track's scratch) is re-run with
   (i) e2e_relic's rule (regression), (ii) x0 interpolated in ln(1 + R) (concave), (iii) x0 + 0.03 (a bracket),
at S/Q = 7, 10, 14.  Output results/sk_pmap.json.  [C] on the inherited conditional recipe.
"""
import json, math, os, sys, time
import numpy as np
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_chain_skeptic'
os.environ['GRUT_KINETIC_WORKDIR'] = SP + '/kwork'
os.makedirs(SP + '/kwork', exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
TRK = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(TRK, 'fx_chain', 'code', 'kin'))
import m1_kin as MK                                            # noqa: E402  read-only
OUT = os.path.join(HERE, '..', 'results', 'sk_pmap.json')


def W_of_SQ(s):                                                # e2e_relic pmap_runs.py, unchanged
    return ((s + 0.00428) / 0.190255) ** (1 / 0.75013) + 0.97494


ROWS = sorted([(r['R'], math.log(r['onset_untilted']['a_over_aref'])) for r in MK.HO['chain']
               if r['label'].startswith('heavy_release') and r['phi_in_equilibrium']])


def x0_rule(R, rule):
    Rs, xs = [r[0] for r in ROWS], [r[1] for r in ROWS]
    if rule == 'linR':
        return float(np.interp(R, Rs, xs))
    if rule == 'ln1pR':
        return float(np.interp(math.log1p(R), [math.log1p(q) for q in Rs], xs))
    if rule == 'linR+0.03':
        return float(np.interp(R, Rs, xs)) + 0.03
    raise ValueError(rule)


def run(sq, rule):
    R = W_of_SQ(sq) / 13.0412 - 1
    name = 'sk_gas_%s_SQ%g' % (rule, sq)
    MK.add(name, 'm1', R, 'gas', start='relaxed', Rtag='interp')
    orig = MK.x0_for
    MK.x0_for = lambda c: x0_rule(c['R'], rule)
    try:
        t0 = time.time(); o = MK.run_case(name)
    finally:
        MK.x0_for = orig
    return dict(S_per_Q=sq, R=R, rule=rule, x0=o['x0'], p=o['p'], stop=o['stop'], dQ=o['checks']['max_charge_error'],
                seconds=time.time() - t0)


if __name__ == '__main__':
    res = dict(method=__doc__, rows=[])
    ref = {7.0: 0.022011356178076858, 10.0: 0.129946303964089, 14.0: 0.3214958560962796}   # e2e_relic pmap.py rows
    for sq in (7.0, 10.0, 14.0):
        for rule in ('linR', 'ln1pR', 'linR+0.03'):
            r = run(sq, rule); r['p_e2e_relic'] = ref[sq]; r['ratio_to_e2e_relic'] = r['p'] / ref[sq]
            res['rows'].append(r)
            print(json.dumps(r), flush=True)
            json.dump(res, open(OUT, 'w'), indent=1)
