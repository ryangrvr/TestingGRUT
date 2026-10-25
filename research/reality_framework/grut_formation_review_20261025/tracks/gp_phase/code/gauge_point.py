#!/usr/bin/env python3
"""gp_phase track (Claude, own code): gauge-point parameters, the inflationary minimum with an
optional extra complex operator, and a wrapper around the compiled release integrator release.c.

Units: H_inf = 4.7e13 GeV = 1. kappa (Codex 'y'), mu real positive by field rephasing.
Extra operator ops (see release.c): 1 S*X^6, 2 S*YX^3, 3 S*Y^2|X|^2, 4 Y*X^3|X|^2.
C = s_op |C| exp(i delta): delta = 0 is aligned with the renormalizable minimum, delta is the
irremovable (rephasing-invariant) CP phase of the extra operator relative to kappa and mu.
"""
import os, subprocess, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, root

HINF = 4.7e13
MRED = 2.435e18
P = dict(lx=1e-3, ly=1.0, ls=1.0, kap=1e-3, M=1e12 / HINF, h=0.1 * 1e12 / HINF, mu=0.5 * 1e12 / HINF, cH=1.0)
SIGN = {0: 1, 1: 1, 2: -1, 3: 1, 4: 1}
DIM = {1: 7, 2: 5, 3: 5, 4: 6}
NAME = {0: 'none', 1: 'Sigma* X^6', 2: 'Sigma* Y X^3', 3: 'Sigma* Y^2 |X|^2', 4: 'Y* X^3 |X|^2'}
SCRATCH = Path(os.environ.get('GP_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/gp_phase/build'))
BIN = SCRATCH / 'release'
SRC = Path(__file__).resolve().parent / 'release.c'


def build():
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        SCRATCH.mkdir(parents=True, exist_ok=True)
        subprocess.run(['gcc', '-O3', '-o', str(BIN), str(SRC), '-lm'], check=True)


def operator_value(op, z):
    X, Y, S = z
    if op == 1: return np.conj(S) * X ** 6
    if op == 2: return np.conj(S) * Y * X ** 3
    if op == 3: return np.conj(S) * Y * Y * abs(X) ** 2
    if op == 4: return np.conj(Y) * X ** 3 * abs(X) ** 2
    return 0j


def potential(z, op=0, C=0j, hubble=0.0, p=P):
    X, Y, S = z
    ax2, ay2, as2 = abs(X) ** 2, abs(Y) ** 2, abs(S) ** 2
    v = (p['lx'] * ax2 ** 2 + p['ly'] * ay2 ** 2 + p['ls'] * (as2 - p['h'] ** 2 / 2) ** 2 + p['M'] ** 2 * ay2
         + 2 * np.real(p['kap'] * np.conj(Y) * X ** 3) + 2 * np.real(p['mu'] * np.conj(S) * Y * Y)
         - hubble * p['cH'] * ax2)
    if op:
        v += 2 * np.real(C * operator_value(op, z))
    return float(np.real(v))


def _unpack(u):
    return np.array([u[0] + 0j, u[1] + 1j * u[2], u[3] + 1j * u[4]])


def force(z, op=0, C=0j, hubble=0.0, p=P):
    """F_i = dV/dz_i^* (same formulas as release.c)."""
    X, Y, S = z
    ax2 = abs(X) ** 2
    Cc = np.conj(C)
    F = np.array([2 * p['lx'] * ax2 * X + 3 * p['kap'] * Y * np.conj(X) ** 2 - hubble * p['cH'] * X,
                  p['M'] ** 2 * Y + 2 * p['ly'] * abs(Y) ** 2 * Y + p['kap'] * X ** 3 + 2 * p['mu'] * S * np.conj(Y),
                  2 * p['ls'] * (abs(S) ** 2 - p['h'] ** 2 / 2) * S + p['mu'] * Y * Y], dtype=complex)
    if op == 1:
        F[0] += 6 * Cc * S * np.conj(X) ** 5; F[2] += C * X ** 6
    elif op == 2:
        F[0] += 3 * Cc * S * np.conj(Y) * np.conj(X) ** 2; F[1] += Cc * S * np.conj(X) ** 3; F[2] += C * Y * X ** 3
    elif op == 3:
        F[0] += 2 * np.real(C * np.conj(S) * Y * Y) * X; F[1] += 2 * Cc * S * np.conj(Y) * ax2; F[2] += C * Y * Y * ax2
    elif op == 4:
        F[0] += C * np.conj(Y) * X ** 4 + 4 * Cc * Y * X * np.conj(X) ** 3; F[1] += C * X ** 3 * ax2
    return F


def _grad(u, op, C, p):
    F = force(_unpack(u), op, C, 1.0, p)
    return np.array([2 * F[0].real, 2 * F[1].real, 2 * F[1].imag, 2 * F[2].real, 2 * F[2].imag])


def inflation_minimum(op=0, C=0j, starts=None, p=P):
    """Minimum of V - c_H H^2 |X|^2 (H=1) over X>0 real (gauge fixed), complex Y, Sigma.
    Returns (z, V_min, gradient norm, number of distinct local minima found)."""
    f = lambda u: potential(_unpack(u), op, C, 1.0, p)
    g = lambda u: _grad(u, op, C, p)
    if starts is None:
        starts = [[23.824, -1.8914, 0.0, -0.26698, 0.0]]
        for th_y in np.linspace(0, 2 * np.pi, 6, endpoint=False):
            for th_s in np.linspace(0, 2 * np.pi, 4, endpoint=False):
                starts.append([23.8, 1.9 * np.cos(th_y), 1.9 * np.sin(th_y), 0.27 * np.cos(th_s), 0.27 * np.sin(th_s)])
    found = []
    for s0 in starts:
        r = minimize(f, s0, jac=g, method='BFGS', options=dict(gtol=1e-11, maxiter=50000))
        u = r.x
        pol = root(g, u, tol=1e-14)
        if pol.success and f(pol.x) <= f(u) + 1e-9 and np.linalg.norm(g(pol.x)) < np.linalg.norm(g(u)):
            u = pol.x
        # Non-renormalizable extra operators make V unbounded at large fields: keep only the local
        # (metastable) minimum in the neighbourhood of the renormalizable one.
        if u[0] < 0 or u[0] > 40 or np.hypot(u[1], u[2]) > 5 or np.hypot(u[3], u[4]) > 2 or not np.isfinite(f(u)):
            continue
        found.append((f(u), u))
    if not found:
        return None, None, None, 0
    found.sort(key=lambda t: t[0])
    distinct = []
    for v, u in found:
        if all(abs(v - w) > 1e-7 for w, _ in distinct):
            distinct.append((v, u))
    v, u = found[0]
    return _unpack(u), v, float(np.linalg.norm(g(u))), len(distinct)


def stationarity(z, op=0, C=0j, p=P):
    """Gradient of V_inf w.r.t. 6 real field components by central differences (scaled)."""
    f = lambda zz: potential(zz, op, C, 1.0, p)
    g = []
    for i in range(3):
        for d in (1.0, 1j):
            e = np.zeros(3, complex); e[i] = d * 1e-6 * max(1.0, abs(z[i]))
            g.append((f(z + e) - f(z - e)) / (2e-6 * max(1.0, abs(z[i]))))
    return np.array(g)


def hessian_eigs(z, op=0, C=0j, p=P, step=1e-4):
    f = lambda u: potential(_unpack(u), op, C, 1.0, p)
    u0 = np.array([z[0].real, z[1].real, z[1].imag, z[2].real, z[2].imag])
    n = 5; Hm = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ej = np.zeros(n); ei[i] = step; ej[j] = step
            Hm[i, j] = Hm[j, i] = (f(u0 + ei + ej) - f(u0 + ei - ej) - f(u0 - ei + ej) + f(u0 - ei - ej)) / (4 * step * step)
    return np.linalg.eigvalsh(Hm)


def run_release(z0, op=0, C=0j, tau_end=1000.0, rtol=1e-11, atol=1e-13, nsamp=60, p=P):
    build()
    args = [str(BIN), str(op), repr(float(np.real(C))), repr(float(np.imag(C))), repr(p['lx']), repr(p['ly']), repr(p['ls']),
            repr(p['kap']), repr(p['mu']), repr(p['M']), repr(p['h']), repr(tau_end), repr(rtol), repr(atol), str(nsamp)]
    args += [repr(float(v)) for zz in z0 for v in (zz.real, zz.imag)]
    out = subprocess.run(args, capture_output=True, text=True)
    rows, end = [], None
    for line in out.stdout.splitlines():
        if line.startswith('S '):
            v = [float(x) for x in line.split()[1:]]
            rows.append(dict(tau=v[0], a3nX=v[1], a3nY=v[2], a3nS=v[3], a4rhoX=v[4], epsX=v[5], gauss_rel=v[6],
                             ledger=v[7], a4rho=v[8], aX=v[9], aY=v[10], aS=v[11]))
        elif line.startswith('END'):
            t = line.split()
            end = dict(steps=int(t[2]), rejects=int(t[4]), max_gauss=float(t[6]), max_ledger=float(t[8]))
        elif line.startswith('FAIL'):
            end = dict(failed_at=float(line.split()[1]))
        elif line.startswith('RUNAWAY'):
            end = dict(runaway_at=float(line.split()[1]))
    return rows, end


def reference_term(op, z, p=P):
    """Magnitude of the renormalizable term an extra operator is compared with at the aligned
    minimum: the kappa term for op 4 (class kappa), the mu term (the only renormalizable Sigma
    phase dependence) for ops 1-3, which are linear in Sigma*."""
    X, Y, S = z
    if op == 4:
        return 2 * p['kap'] * abs(Y) * abs(X) ** 3
    return 2 * p['mu'] * abs(S) * abs(Y) ** 2


def coefficient_for_ratio(op, r, z_aligned, p=P):
    """|C| such that |2 C O| = r * reference_term at the aligned renormalizable minimum."""
    return r * reference_term(op, z_aligned, p) / (2 * abs(operator_value(op, z_aligned)))


def ratio_for_coefficient(op, Cabs, z_aligned, p=P):
    return 2 * Cabs * abs(operator_value(op, z_aligned)) / reference_term(op, z_aligned, p)


def coefficient_from_wilson(op, c, Lam_GeV=MRED):
    d = DIM[op]
    return c / Lam_GeV ** (d - 4) * HINF ** (d - 4)


def wilson(op, Cabs):
    """Convert |C| (H units) to c at Lambda = M_red, and to Lambda (GeV) for c = 1."""
    d = DIM[op]; k = d - 4
    C_gev = Cabs * HINF ** (4 - d)
    return dict(c_at_Mred=C_gev * MRED ** k, Lambda_for_c1_GeV=C_gev ** (-1.0 / k))
