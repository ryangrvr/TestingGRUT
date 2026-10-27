"""Thin Python wrapper around our SMDR driver (smdr_driver.c). Verifier v_vacua.

SMDR (Martin-Robertson, arXiv:1907.02500, v1.3) is built in the scratchpad by v3_build_smdr.sh; set SMDR_HOME to
override the location.  Conventions (SMDR): V = m2|H|^2 + k|H|^4, k = lambda (same normalization as Buttazzo et al.
and as track `vacua`), beta_X = Q dX/dQ.
"""
import os
import subprocess
import numpy as np

SMDR_HOME = os.environ.get(
    'SMDR_HOME',
    '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/v_vacua_dl/smdr/SMDR-1.3')
EXE = os.path.join(SMDR_HOME, 'smdr_driver')
MPL = 1.2209e19
MPL_RED = 2.435e18


def _fmt(a):
    if isinstance(a, (int, np.integer)):
        return str(int(a))
    if isinstance(a, (float, np.floating)):
        return format(float(a), '.17g')
    return str(a)


def _call(args):
    out = subprocess.run([EXE] + [_fmt(a) for a in args], cwd=SMDR_HOME, capture_output=True, text=True, check=True)
    return out.stdout


def _parse_kv(line):
    tok = line.split()
    return {tok[i]: float(tok[i + 1]) for i in range(1, len(tok) - 1, 2)}


def betas(L, g3, g, gp, yt, yb, ytau, k):
    """SMDR beta functions Q dX/dQ at loop order L (all light Yukawas zero)."""
    txt = _call(['betas', L, (g3), (g), (gp), (yt), (yb), (ytau), (k)])
    d = {}
    for ln in txt.strip().splitlines():
        a, b = ln.split()
        d[a] = float(b)
    return d


def run(L, Q0, c, N=241, Qmin=1e8, Qmax=MPL * 1.0001):
    """Run MSbar couplings c (dict g3,g,gp,yt,yb,ytau,k at Q0) with SMDR at loop order L; return arrays."""
    txt = _call(['run', L, (Q0), (c['g3']), (c['g']), (c['gp']), (c['yt']), (c['yb']),
                 (c['ytau']), (c['k']), N, (Qmin), (Qmax)])
    rows = [_parse_kv(ln) for ln in txt.strip().splitlines() if ln.startswith('pt')]
    return {key: np.array([r[key] for r in rows]) for key in rows[0]}


def fit(Mt, Mh, als):
    """SMDR full-precision fit: MSbar couplings at Q = Mt from pole masses Mt, Mh and alpha_s^(5)(MZ)."""
    txt = _call(['fit', (Mt), (Mh), (als)])
    res = {}
    for ln in txt.strip().splitlines():
        if ln.startswith('fit '):
            res.update(_parse_kv(ln))
        elif ln.startswith('fitMW'):
            res.update(_parse_kv(ln))
    return res


def min_lambda(L, Q0, c, Qcut=MPL, N=241):
    """min over 1e8 <= Q <= Qcut of k(Q): returns (k_min, Q_star, interior?) using a parabola in ln Q."""
    r = run(L, Q0, c, N=N, Qmin=1e8, Qmax=Qcut)
    lnQ, k = np.log(r['Q']), r['k']
    i = int(np.argmin(k))
    if 0 < i < len(k) - 1:
        x = lnQ[i - 1:i + 2]
        y = k[i - 1:i + 2]
        A = np.polyfit(x, y, 2)
        xs = -A[1] / (2 * A[0])
        return float(np.polyval(A, xs)), float(np.exp(xs)), True
    return float(k[i]), float(np.exp(lnQ[i])), False
