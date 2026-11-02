#!/usr/bin/env python3
"""inf_chain (Claude, own driver): the pair fraction p = X-bar per net charge after freeze-out at LARGE hidden heat,
from the reviewed case-L conditional moment network (Codex 34f23f8 moment_controls over our kin_solve modules), as copied
and M1-patched in tracks/fx_chain/code/kin (m1_kin.py; imported READ-ONLY, nothing there is modified or written: the
work directory is redirected to this track's scratch).  Gas-only recipe (rule a), M1 inputs, start at a_ref
(forced chemical + kinetic equilibrium), energy rho_c (1 + R) in spot units, i.e. W = 13.0412 (1 + R).
This extends the p(S/Q) map used by dq_critic (wt_critic gas-only map to S/Q = 6.74 plus the REVIEW_KINETIC anchor at
S/Q = 21.8) to S/Q ~ 10^2-10^3.  [C] on a conditional recipe (forced equilibrium start, no objects, MB fugacities);
the recipe's validity at these temperatures is not established [E].
Usage: python3 -B pmap_kin.py R1 [R2 ...]   -> merges into ../results/pmap_kin.json
"""
import json, math, os, sys, time
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
os.environ['GRUT_KINETIC_WORKDIR'] = SP + '/wt2/inf_chain/kwork'
HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(HERE, '..', '..', 'fx_chain', 'code', 'kin')
sys.path.insert(0, KIN)
import m1_kin as MK          # noqa: E402  read-only
OUT = os.path.join(HERE, '..', 'results', 'pmap_kin.json')


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def main(Rs):
    try:
        allr = json.load(open(OUT))
    except Exception:
        allr = dict(method=__doc__, rows={})
    for R in Rs:
        name = 'chain_gas_aref_R%g' % R
        MK.add(name, 'm1', R, 'gas', start='aref', Rtag=None)
        t0 = time.time()
        try:
            o = MK.run_case(name)
        except Exception as e:
            print('FAIL R=%g %r' % (R, e), flush=True); allr['rows'][str(R)] = dict(R=R, error=repr(e)); continue
        W = 13.0412 * (1 + R)
        f = o['final']
        row = dict(R_spot=R, W=W, S_per_Q_closed_form=SQ(W), p=o['p'], stop=o['stop'], T_over_m_stop=f['T_over_m'],
                   visible_T_GeV_stop=f['visible_T_GeV'], initial_t=o['initial']['t'], cond_fraction0=o['initial']['cond_fraction'],
                   max_charge_error=o['checks']['max_charge_error'], max_energy_work_error=o['checks']['max_energy_work_error'],
                   seconds=time.time() - t0)
        allr['rows'][str(R)] = row
        print('R=%g W=%.1f S/Q=%.2f p=%.4g stop=%s T0/m=%.3f dQ=%.1e dE=%.1e %.0fs' % (R, W, row['S_per_Q_closed_form'], o['p'], o['stop'],
              o['initial']['t'], row['max_charge_error'], row['max_energy_work_error'], row['seconds']), flush=True)
        json.dump(allr, open(OUT, 'w'), indent=1, default=float)




# ---- self-consistent charge yield (added): with M = 1 + 2p the calibrated charge yield is Y_Q,M1 x 0.9374/M, so the
# visible clock of the network (B.T_of_charge uses K.YQ) is run with K.YQ scaled by 0.9374/M, iterated with p.
_conf0 = MK.configure
YQ_SCALE = dict(f=1.0)


def _conf(variant):
    _conf0(variant)
    MK.K.YQ = MK.K.YQ * YQ_SCALE['f']


MK.configure = _conf


def W_of_SQ(s):
    return ((s + 0.00428) / 0.190255) ** (1 / 0.75013) + 0.97494


def self_consistent(SQt, iters=4):
    """p at S/Q = SQt with K.YQ = Y_Q(M1) * 0.9374/(1 + 2p), fixed-point iteration from the fixed-yield value."""
    R = W_of_SQ(SQt) / 13.0412 - 1
    hist = []; p = None
    for it in range(iters):
        M = 1 + 2 * p if p is not None else 1.0
        YQ_SCALE['f'] = 0.9374 / M if p is not None else 1.0
        name = 'chain_sc_SQ%g_it%d' % (SQt, it)
        MK.add(name, 'm1', R, 'gas', start='aref', Rtag=None)
        o = MK.run_case(name)
        pn = o['p']
        hist.append(dict(it=it, YQ_factor=YQ_SCALE['f'], p=pn, stop=o['stop'], visible_T_GeV_stop=o['final']['visible_T_GeV']))
        if p is not None and abs(pn / p - 1) < 0.01:
            p = pn; break
        p = pn if p is None else math.sqrt(p * pn)          # geometric damping
    YQ_SCALE['f'] = 1.0
    return dict(S_per_Q=SQt, R_spot=R, W=W_of_SQ(SQt), p_fixed_yield=hist[0]['p'], p_self_consistent=hist[-1]['p'], history=hist)


if __name__ == '__main__':
    if sys.argv[1] == 'sc':
        try:
            allr = json.load(open(OUT))
        except Exception:
            allr = dict(method=__doc__, rows={})
        allr.setdefault('self_consistent', {})
        for s in sys.argv[2:]:
            r = self_consistent(float(s))
            allr['self_consistent'][s] = r
            print('S/Q=%s: p fixed-yield %.4g -> self-consistent %.4g  %s' % (s, r['p_fixed_yield'], r['p_self_consistent'],
                  [(h['YQ_factor'], round(h['p'], 4)) for h in r['history']]), flush=True)
            json.dump(allr, open(OUT, 'w'), indent=1, default=float)
    else:
        main([float(x) for x in sys.argv[1:]])
