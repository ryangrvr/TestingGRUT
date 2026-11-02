#!/usr/bin/env python3
"""inf_potential SKEPTIC (Claude): re-run the track's engine copy rpinf.c (built into the skeptic scratch, source read-only)
on the M1r mode-0 reference and on the MODE-4 benchmarks, and record the time series of the homogeneous mediator
diagnostics N_YH, N_YL, N_sigma (columns 10, 11, 12 of the S rows), N_X and E_X.  Purpose: the track's smoke table reports
N_YL only; engine_test.json also stores N_YH, which the NOTES and the rough re-ledger omit.
Arguments are those of the track's engine_test.py (copied here so nothing in the track is executed or written).
Read-only imports: tracks/dq_heat/skeptic/code/dq_model_copy.py (stationary preparation), track bench_*.json.
Usage: python3 -B sk_engine.py M1r HT6t T1 E1 [tend]   Output: ../results/sk_engine.json (merged)
"""
import json, math, os, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]                    # tracks/inf_potential
WT = TRK.parent                          # tracks/
sys.path.insert(0, str(WT / 'dq_heat' / 'skeptic' / 'code'))
import dq_model_copy as fc               # read-only
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_potential_skeptic/build')
OUT = HERE.parent / 'results'
MP = 2.435e18


def build():
    SP.mkdir(parents=True, exist_ok=True)
    subprocess.run(['gcc', '-O2', '-o', str(SP / 'rpinf'), str(TRK / 'code' / 'rpinf.c'), '-lm'], check=True)


def spec_args(z):
    out = ['F']
    for c in z:
        out += ['%.17g' % c.real, '%.17g' % c.imag]
    return out + ['0'] * 6


def common_args(mode, b, gam, eta2, p, G, end, nsamp):
    return [str(mode), '%.17g' % b, '%.17g' % gam, '%.17g' % p['xi'], '%.17g' % eta2, '0', '%.17g' % p['lx'], '%.17g' % p['ly'],
            '%.17g' % p['ls'], '%.17g' % p['kap'], '%.17g' % p['M'], '%.17g' % p['h'], '%.17g' % p['mu'], '0', '0',
            '%.17g' % G.real, '%.17g' % G.imag, '0', '%.17g' % end, '1e-10', '1e-13', str(nsamp), '0', '0']


def parse_S(txt):
    return np.array([list(map(float, l.split()[1:])) for l in txt.splitlines() if l.startswith('S ')])


def run(name, tend, nsamp):
    p = fc.params(0.5)
    z = fc.stationary(fc.U_PREP, p, R=12.0)['z']
    env = dict(os.environ)
    if name == 'M1r':
        eta2 = fc.ETA2; I0 = math.sqrt(6.0) / 3.0
        G = (fc.U_PREP - p['mu']) / (I0 / math.sqrt(eta2))
        args = common_args(0, 3.0, 0.5, eta2, p, G, tend, nsamp) + spec_args(z)
        norm3 = 1.0
    else:
        rec = json.load(open(TRK / 'results' / ('bench_%s.json' % name)))['by_Hi'][0]
        Hi = rec['H_i_GeV']; e2 = (Hi / MP) ** 2; V0 = rec['V0_MP4']
        par = rec['params']
        if name == 'HT6t':
            kind = 3; vev = par['vev']; pw = par['p']; bta = par['beta']
            xt = -(bta / pw) ** (1 / (pw - 1)); ftop = 1 - bta * xt - xt ** pw
            env.update(RPINF_P1='%.17g' % vev, RPINF_P2='%d' % pw, RPINF_P3='%.17g' % bta, RPINF_FTOP='%.17g' % ftop)
        else:
            kind = 1 if name == 'T1' else 2
            M = math.sqrt(6 * par['alpha']) if name == 'T1' else math.sqrt(1.5 * par['alpha'])
            env.update(RPINF_P1='%.17g' % M)
        env.update(RPINF_KIND=str(kind), RPINF_AMP='%.17g' % (V0 / e2), RPINF_PHIMIN='%.17g' % rec['phi_min_MP'],
                   RPINF_I0='%.17g' % rec['phiI_S1_MP'], RPINF_SIGN='%d' % (1 if rec['phidot_S1_over_HiMP'] > 0 else -1))
        G = (fc.U_PREP - p['mu']) / (rec['phiI_S1_MP'] / math.sqrt(e2))
        args = common_args(4, rec['m_I_over_Hi'], 0.5, e2, p, G, tend, nsamp) + spec_args(z)
        norm3 = math.e ** 3                 # S1 -> eps_H = 1 normalisation for comoving numbers
    res = subprocess.run([str(SP / 'rpinf')] + args, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        raise RuntimeError(res.stderr[:400])
    S = parse_S(res.stdout)
    rows = []
    for r in S[::max(1, len(S) // 60)]:
        rows.append(dict(tau=r[0], a=r[1], NX=r[4] / norm3, EX=r[8] / (norm3 ** (4 / 3)), eps=r[9], NYH=r[10] / norm3,
                         NYL=r[11] / norm3, Nsig=r[12] / norm3, absY=r[14], S_over_vev=r[15], ure=r[26], uim=r[27], mL2=r[29]))
    late = S[S[:, 0] >= 0.66 * tend]
    summ = dict(NX=float(np.mean(late[:, 4])) / norm3, EX=float(np.mean(late[:, 8])) / norm3 ** (4 / 3), eps=float(np.mean(late[:, 9])),
                NYH=float(np.mean(late[:, 10])) / norm3, NYL=float(np.mean(late[:, 11])) / norm3, Nsig=float(np.mean(late[:, 12])) / norm3,
                NYH_spread=float(np.std(late[:, 10]) / max(1e-300, abs(np.mean(late[:, 10])))),
                NYL_spread=float(np.std(late[:, 11]) / max(1e-300, abs(np.mean(late[:, 11])))))
    return dict(benchmark=name, tend=tend, normalisation='a = 1 at eps_H = 1 (benchmarks) / at the M1r slice (M1r)', late_mean=summ, series=rows)


def main(names, tend):
    build()
    try:
        allr = json.load(open(OUT / 'sk_engine.json'))
    except Exception:
        allr = dict(method=__doc__, runs={})
    for n in names:
        r = run(n, tend, 600)
        allr['runs']['%s_t%g' % (n, tend)] = r
        s = r['late_mean']
        print('%-5s tend=%g: NX=%.5g EX=%.5g eps=%.4f NYH=%.5g NYL=%.5g Nsig=%.4g (spreads %.2g %.2g)' % (
            n, tend, s['NX'], s['EX'], s['eps'], s['NYH'], s['NYL'], s['Nsig'], s['NYH_spread'], s['NYL_spread']), flush=True)
        for row in r['series'][:: max(1, len(r['series']) // 12)]:
            print('     tau=%7.2f NYH=%.4g NYL=%.4g Nsig=%.3g |Y|=%.3g u=%.3f%+.3fi mL2=%.1f' % (
                row['tau'], row['NYH'], row['NYL'], row['Nsig'], row['absY'], row['ure'], row['uim'], row['mL2']), flush=True)
        json.dump(allr, open(OUT / 'sk_engine.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    args = sys.argv[1:]
    tend = float(args[-1]) if args and args[-1].replace('.', '').isdigit() else 60.0
    names = [a for a in args if not a.replace('.', '').isdigit()] or ['M1r', 'HT6t']
    main(names, tend)
