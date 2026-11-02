#!/usr/bin/env python3
"""inf_chain SKEPTIC (Claude, own code): sensitivity of the Gamma_I = 0.5 releases to WHEN the phenomenological width acts.

The track's engine (rpchain.c MODE 4, started at S1) applies the friction (3H + Gamma) phi' and the radiation source
Gamma phi'^2 from S1, i.e. during the last e-fold of slow roll.  A perturbative decay width acts on oscillating
quanta, not on a slowly rolling condensate (there it is at most a warm-inflation-type dissipation, which needs its
own justification).  The potential track's Python history switched it on at eps_H = 1.  NOTES item 6 records the
10-11% K difference as 'a convention'; here I measure what the same switch does to the RELEASE (N_X, N_YH, N_YL, N_sigma,
E_X), which the NOTES did not test.

Method: a patched copy of the track's engine (generated into my scratch from tracks/inf_chain/code/rpchain.c, read-only):
Gamma_eff = Gamma for ln a >= SKG_LNAON, 0 before (ln a = 1 is eps_H = 1 to < 1e-6 on these backgrounds, a = 1 at S1).
Default SKG_LNAON = -inf reproduces rpchain.c.  The track's own driver (ch_common, imported read-only) prepares the
state; nothing is written into the track's results or scratch (CH_SCRATCH is redirected to mine).
Output ../results/sk_gamma_switch.json.  Label [C] (engine re-run with one stated patch).
"""
import json, os, re, sys, subprocess, pathlib
HERE = pathlib.Path(__file__).resolve().parent
SP = pathlib.Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_chain_skeptic')
os.environ['CH_SCRATCH'] = str(SP)
sys.dont_write_bytecode = True
TRK = HERE.parent.parent.parent
sys.path.insert(0, str(TRK / 'inf_chain' / 'code'))
import numpy as np
import ch_common as ch          # read-only driver
from benchmarks_reg import get_bench


def make_source():
    src = open(TRK / 'inf_chain' / 'code' / 'rpchain.c').read()
    src = src.replace('static double b, gam, xi,', 'static double LNAON = -1e300;\nstatic double b, gam, xi,', 1)
    # rhs: gated width
    a = 'static void rhs(const double *y, double *d) {\n  geo g; geometry(y, &g);\n  double H = g.H;\n'
    assert a in src
    i0 = src.index(a); i1 = src.index('\n}\n', i0)
    body = src[i0 + len(a):i1]
    body = re.sub(r'\bgam\b', 'ge', body)
    src = src[:i0] + a + '  double ge = (y[0] >= LNAON) ? gam : 0.0;\n' + body + src[i1:]
    # diag: gated width in Rdot
    b0 = 'static void diag(double t, const double *y, double *q) {\n  geo g; geometry(y, &g);\n'
    assert b0 in src
    j0 = src.index(b0); j1 = src.index('\n}\n', j0)
    body = src[j0 + len(b0):j1]
    body = re.sub(r'\bgam\b', 'ge', body)
    src = src[:j0] + b0 + '  double ge = (y[0] >= LNAON) ? gam : 0.0;\n' + body + src[j1:]
    # env
    k = src.index('double capT = getenv("RPCH_CAPT")')
    src = src[:k] + 'LNAON = getenv("SKG_LNAON") ? atof(getenv("SKG_LNAON")) : -1e300;\n  ' + src[k:]
    (SP / 'build').mkdir(parents=True, exist_ok=True)
    p = SP / 'build' / 'sk_rpgate.c'
    open(p, 'w').write(src)
    binp = SP / 'build' / 'sk_rpgate'
    subprocess.run(['gcc', '-O2', '-o', str(binp), str(p), '-lm'], check=True)
    return binp


def run(b, lnaon, T=100.0):
    env = b.env()
    if lnaon is not None:
        env['SKG_LNAON'] = '%.17g' % lnaon
    z = b.prep()['z']
    args = ch.engine_args(b, b.gam, T, int(round(T / 0.005)), z=z)
    args[0] = str(BIN)
    r = subprocess.run(args, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[:400])
    S, C, END, I0, W = ch.parse(r.stdout)
    late = S[S[:, 0] >= 0.75 * T]
    out = {nm: float(np.mean(late[:, ch.COLS.index(nm)])) for nm in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig', 'a', 'H')}
    out['K_late'] = float(np.mean(late[:, ch.COLS.index('a')] ** 2 * late[:, ch.COLS.index('H')]))
    out['END'] = END
    # when eps_H = 1 is crossed: w_I = -1/3 crossing of the inflaton
    return out


BIN = None


def main():
    global BIN
    BIN = make_source()
    res = {}
    for key in sys.argv[1:] or ['HT6t', 'T1']:
        b = get_bench(key)
        r0 = run(b, None)
        r1 = run(b, 1.0)
        tr = json.load(open(TRK / 'inf_chain' / 'results' / ('release_%s.json' % key)))['dense']['late_mean']
        rel = {k: r1[k] / r0[k] - 1 for k in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig', 'K_late')}
        reg = {k: r0[k] / tr[k] - 1 for k in ('NX', 'EX', 'eps', 'NYL', 'NYH', 'Nsig')}
        res[key] = dict(gamma_from_S1=r0, gamma_from_end=r1, rel_change=rel, regression_vs_track_dense=reg)
        print(key, 'regression', reg, flush=True)
        print(key, 'switch at eps_H=1 vs S1:', rel, flush=True)
    p = HERE.parent / 'results' / 'sk_gamma_switch.json'
    json.dump(res, open(p, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
