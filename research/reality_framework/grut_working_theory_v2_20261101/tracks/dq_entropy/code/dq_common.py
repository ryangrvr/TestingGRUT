#!/usr/bin/env python3
"""dq_entropy track (Claude, own code; issue #13 decisive pass, 2 November 2026): driver of rpdq.c.

Reuse, stated: the model constants, the exact-Hessian stationary point (fx_common.stationary), the trilinear coupling
(fx_common.trilinear_G), the yield convention (fx_common.yield_over_target) and the decay ledger (fx_common.ledger) are
imported READ-ONLY from tracks/fx_holder/code/fx_common.py.  The engine is our copy rpdq.c of fx_holder's rptri.c (changes
DQ1-DQ3 in its header).  run() below is a copy of fx_common.run with the new arguments (cport, cbin) and the conformal bins.
"""
import os, subprocess, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / 'fx_holder' / 'code'))
import fx_common as fc  # noqa: E402  (read-only reuse)

SP = Path(os.environ.get('DQ_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_entropy'))
BIN = SP / 'build' / 'rpdq'
SRC = HERE / 'rpdq.c'
RES = HERE.parent / 'results'
COLS = fc.COLS
CP_UNIT = 1.0 / (4 * np.pi ** 2)                     # one unit of ln Q of the induced-portal running
LN_MP_MY = np.log(2.435e18 / (30 * fc.H_GEV))        # ln(M_P/M_Y) at the M1 H_i = 13.31


def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-std=c11', '-o', str(BIN), str(SRC), '-lm'], check=True)


def portal_params(p, cport, u):
    """effective (M, h) that absorb the portal at fixed u: M^2 -> M^2 + 2 du, h^2/2 -> h^2/2 - du/(2 ls)."""
    du = cport * (abs(u) ** 2 - p['mu'] ** 2)
    q = dict(p)
    q['M'] = float(np.sqrt(p['M'] ** 2 + 2 * du))
    q['h'] = float(np.sqrt(p['h'] ** 2 - du / p['ls']))
    return q, du


def prepare(lam=0.5, b=3.0, cport=0.0):
    """held preparation at u = 10 i with the portal counterterm included (stationary point of V + xi R |X|^2 + portal)."""
    p = fc.params(lam)
    q, du = portal_params(p, cport, fc.U_PREP)
    st = fc.stationary(fc.U_PREP, q)
    G, I0, ph0 = fc.trilinear_G(st['z'], p, b)
    return p, st, G, du


def run(tag, p, b, gam, z0, G=0j, cport=0.0, end=1000.0, rtol=1e-10, atol=1e-13, nsamp=5000, binw=20.0, cbin=0.0, v0=None, save=True):
    build()
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), '0', g(b), g(gam), g(p['xi']), g(fc.ETA2), g(p['cH']), g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']), g(p['M']),
            g(p['h']), g(p['mu']), g(0.0), g(0.0), g(G.real), g(G.imag), g(cport), g(end), g(rtol), g(atol), str(int(nsamp)), g(binw), g(cbin)]
    v0 = np.zeros(3, complex) if v0 is None else np.asarray(v0, complex)
    args += ['F'] + [g(v) for c in z0 for v in (c.real, c.imag)] + [g(v) for c in v0 for v in (c.real, c.imag)]
    t0 = time.time()
    out = subprocess.run(args, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(tag + ': ' + out.stderr)
    S, Bn, Cn, W, Y, END, I0, ETA = [], [], [], {}, None, None, None, None
    for line in out.stdout.split('\n'):
        if line.startswith('S '):
            S.append([float(v) for v in line.split()[1:]])
        elif line.startswith('B '):
            Bn.append([float(x) for x in line.split()[2:]])
        elif line.startswith('C '):
            Cn.append([float(x) for x in line.split()[2:]])
        elif line.startswith('W') or line.startswith('A3040'):
            v = line.split(); W[v[0]] = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('Y '):
            Y = [float(v) for v in line.split()[1:]]
        elif line.startswith('ETA '):
            ETA = float(line.split()[1])
        elif line.startswith('I0 '):
            v = [float(x) for x in line.split()[1:]]; I0 = dict(I0=v[0], phihat0=v[1], u0=[v[2], v[3]], R0=v[4])
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_H_consistency=float(v[4]),
                                        max_hidden_work=float(v[5]), min_mL2=float(v[6]))
    S = np.array(S); Bn = np.array(Bn); Cn = np.array(Cn)
    rr = dict(tag=tag, S=S, bins=Bn, cbins=Cn, W=W, Y=Y, END=END, I0=I0, eta_end=ETA, seconds=time.time() - t0, lam=p['ly'], b=b,
              gam=gam, end=end, G=[G.real, G.imag], cport=cport, rtol=rtol)
    if save:
        (SP / 'runs').mkdir(parents=True, exist_ok=True)
        np.savez_compressed(SP / 'runs' / (tag + '.npz'), S=S, bins=Bn, cbins=Cn)
    return rr


def bin_col(B, name, off=2):
    return B[:, off + COLS.index(name)]


def window_mean(B, name, lo, hi, off=2, key=0):
    """mean over bins whose [x0, x1] lies inside [lo, hi]; key = column of x0 (0) and x1 (1)."""
    x0, x1 = B[:, key], B[:, key + 1]; vals = bin_col(B, name, off); dx = x1 - x0
    sel = (x0 >= lo - 1e-9) & (x1 <= hi + 1e-9)
    return float(np.sum(vals[sel] * dx[sel]) / np.sum(dx[sel]))
