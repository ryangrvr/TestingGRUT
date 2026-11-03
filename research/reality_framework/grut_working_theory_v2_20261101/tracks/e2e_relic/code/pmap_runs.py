#!/usr/bin/env python3
"""e2e_relic item 1 (Claude, own driver): the reviewed case-L conditional moment network (fx_chain's M1-patched copy of
Codex 34f23f8 moment_controls via m1_kin.py; imported READ-ONLY, work directory redirected to this track's scratch),
gas-only rule, run with BOTH start conventions inside and around the pair-fraction-map gap 6.74 < S/Q < 17.1:
  * 'relaxed' : start at the relaxed onset a_on, x0 = ln(a_on/a_ref) linearly interpolated in R over fx_chain's heat_onset
                rows with phi in equilibrium (R = 0 ... 49.6; the gap R = 8-30 lies between the computed rows R = 5.70
                and R = 41.7, so this is an interpolation, not an extrapolation).  This is the convention of the
                wt_critic/fx_chain gas map (S/Q <= 6.74), i.e. of the M1r map below the gap.
  * 'aref'    : forced equilibrium at a_ref (x0 = 0): inf_chain pmap_kin's and the mf_chain skeptic's convention (>= 6.74).
R from S/Q by the dq_entropy closed form: R = W(S/Q)/13.0412 - 1.  [C] on a conditional recipe (forced equilibrium at the
start, no objects, MB fugacities), inherited, not validated.
Usage: python3 -B pmap_runs.py relaxed|aref SQ1 [SQ2 ...]  -> merges into ../results/pmap_runs.json
"""
import json, math, os, sys, time
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_relic'
os.environ['GRUT_KINETIC_WORKDIR'] = SP + '/kwork'
HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(HERE, '..', '..', 'fx_chain', 'code', 'kin')
sys.path.insert(0, KIN)
import m1_kin as MK          # noqa: E402  read-only
OUT = os.path.join(HERE, '..', 'results', 'pmap_runs.json')
# optional visible-clock yield factor (inf_chain pmap_kin's patch): K.YQ -> K.YQ * f, env E2E_YQ_FACTOR (default 1)
YQF = float(os.environ.get('E2E_YQ_FACTOR', '1'))
_conf0 = MK.configure


def _conf(variant):
    _conf0(variant)
    MK.K.YQ = MK.K.YQ * YQF


MK.configure = _conf


def W_of_SQ(s):
    return ((s + 0.00428) / 0.190255) ** (1 / 0.75013) + 0.97494


if __name__ == '__main__':
    start = sys.argv[1]
    try:
        allr = json.load(open(OUT))
    except Exception:
        allr = dict(method=__doc__, rows={})
    for s in sys.argv[2:]:
        sq = float(s)
        R = W_of_SQ(sq) / 13.0412 - 1
        name = 'e2e_gas_%s_SQ%g_yq%g' % (start, sq, YQF)
        if start == 'relaxed':
            MK.add(name, 'm1', R, 'gas', start='relaxed', Rtag='interp')
        else:
            MK.add(name, 'm1', R, 'gas', start='aref', Rtag=None)
        t0 = time.time()
        o = MK.run_case(name)
        row = dict(start=start, YQ_factor=YQF, S_per_Q=sq, R_spot=R, x0=o['x0'], p=o['p'], stop=o['stop'], T0_over_m=o['initial']['t'],
                   cond_fraction0=o['initial']['cond_fraction'], max_charge_error=o['checks']['max_charge_error'],
                   max_energy_work_error=o['checks']['max_energy_work_error'], seconds=time.time() - t0)
        allr['rows']['%s_%s' % (start, s) + ('' if YQF == 1 else '_yq%g' % YQF)] = row
        print('yq%g %s S/Q=%.3f R=%.3f x0=%+.4f p=%.5g stop=%s cond0=%.3f dQ=%.1e dE=%.1e %.0fs' % (
            YQF, start, sq, R, o['x0'], o['p'], o['stop'], row['cond_fraction0'], row['max_charge_error'],
            row['max_energy_work_error'], row['seconds']), flush=True)
        json.dump(allr, open(OUT, 'w'), indent=1, default=float)
