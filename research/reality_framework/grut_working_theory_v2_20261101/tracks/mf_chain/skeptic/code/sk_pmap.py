#!/usr/bin/env python3
"""mf_chain SKEPTIC: fill the gap of the inherited pair-fraction map p(S/Q) between wt_critic's gas map (S/Q <= 6.74) and
inf_chain's case-L runs (S/Q >= 17.1).  The exported/unblocked HF relics (S/Q(1e-6) = 14.5-15.9) and the isocurvature
derivative D sit INSIDE that gap, where chain_mfi.p_of_SQ log-linearly interpolates (slope d ln p/d ln(S/Q) = 3.68).

Method: the reviewed case-L conditional moment network (fx_chain's M1-patched copy of Codex 34f23f8 moment_controls, via
m1_kin.py; imported READ-ONLY, work directory redirected to this skeptic's scratch), gas-only rule, a_ref start, exactly
as inf_chain/code/pmap_kin.py, at R_spot = W(S/Q)/13.0412 - 1.  [C] on a conditional recipe (forced-equilibrium start,
no objects, MB fugacities), inherited, not validated.
Usage: python3 -B sk_pmap.py SQ1 [SQ2 ...]  -> ../results/sk_pmap.json
"""
import json, math, os, sys, time
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_chain_skeptic'
os.environ['GRUT_KINETIC_WORKDIR'] = SP + '/kwork'
HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(HERE, '..', '..', '..', 'fx_chain', 'code', 'kin')
sys.path.insert(0, KIN)
import m1_kin as MK          # noqa: E402  read-only
OUT = os.path.join(HERE, '..', 'results', 'sk_pmap.json')


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def W_of_SQ(s):
    return ((s + 0.00428) / 0.190255) ** (1 / 0.75013) + 0.97494


if __name__ == '__main__':
    try:
        allr = json.load(open(OUT))
    except Exception:
        allr = dict(method=__doc__, rows={})
    for s in sys.argv[1:]:
        sq = float(s)
        R = W_of_SQ(sq) / 13.0412 - 1
        name = 'sk_gas_aref_SQ%g' % sq
        MK.add(name, 'm1', R, 'gas', start='aref', Rtag=None)
        t0 = time.time()
        o = MK.run_case(name)
        row = dict(S_per_Q=sq, R_spot=R, W=13.0412 * (1 + R), p=o['p'], stop=o['stop'], T0_over_m=o['initial']['t'],
                   cond_fraction0=o['initial']['cond_fraction'], max_charge_error=o['checks']['max_charge_error'],
                   max_energy_work_error=o['checks']['max_energy_work_error'], seconds=time.time() - t0)
        allr['rows'][s] = row
        print('S/Q=%.3f R=%.3f p=%.5g stop=%s cond0=%.3f dQ=%.1e %.0fs' % (sq, R, o['p'], o['stop'], row['cond_fraction0'],
              row['max_charge_error'], row['seconds']), flush=True)
        json.dump(allr, open(OUT, 'w'), indent=1, default=float)
