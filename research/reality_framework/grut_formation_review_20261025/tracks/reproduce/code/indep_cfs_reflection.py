#!/usr/bin/env python3
"""INDEPENDENT check (own code) of the discrete causal-CFS layer reflection claimed in ff49e5d.

Own operator, assembled from the equations stated in MATHEMATICS.md / radial_cfs.cpp:
   s^2 U = C_n D( C_f G U ) - m^2 U,   C = (s+alpha)/(s+alpha+sigma),  s=-i omega,
G = staggered 4th-order gradient (odd node ghosts), D = its divergence (even face ghosts),
cubic sigma = smax ((r-R)/W)^3 on nodes and faces, wall at r=L.

Method (different from Codex's exact-root lead matching): drive the time-harmonic problem by a
point source inside a long physical lead, then least-squares fit U_j = a e^{ikjh} + b e^{-ikjh} on a
window far from both source and layer entrance; |R|^2 = |b/a|^2.  Lattice phase k, lattice momentum
p(k) = (2/h) sin(kh/2)(1+sin^2(kh/2)/6), omega^2 = m^2 + p^2.

Bridge test (optional, argv[2] = path of a compiled operator_controls binary built from Codex's
code): the U-response of the executable's full linear time-domain RHS matrix, (M + i omega)Y = -F,
is compared with this operator's (A + omega^2) U = -delta for a small grid.
Usage: python3 indep_cfs_reflection.py OUT.json [operator_controls_binary]
"""
import io, json, subprocess, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


def lattice_p(k, h):
    s = np.sin(k * h / 2)
    return 2 * s / h * (1 + s * s / 6)


def build(n, h, sig_node, sig_face, s, alpha, mass):
    """Return sparse A (size (n-1)x(n-1)) acting on interior nodes 1..n-1."""
    # G: faces 0..n-1 ; face i between nodes i and i+1
    def node_index(j):  # returns (index in 1..n-1 or None, sign)
        if j < 0:
            j, sg = -j, -1.
        else:
            sg = 1.
        if j > n:
            j, sg = 2 * n - j, -sg
        if j == 0 or j == n:
            return None, 0.
        return j, sg
    rows, cols, vals = [], [], []
    for i in range(n):
        for j, c in [(i + 1, 9 / 8), (i, -9 / 8), (i + 2, -1 / 24), (i - 1, 1 / 24)]:
            jj, sg = node_index(j)
            if jj is not None:
                rows.append(i); cols.append(jj - 1); vals.append(sg * c / h)
    G = sp.csr_matrix((vals, (rows, cols)), shape=(n, n - 1), dtype=complex)

    def face_index(f):
        if f < 0:
            f = -1 - f
        if f >= n:
            f = 2 * n - 1 - f
        return f
    rows, cols, vals = [], [], []
    for i in range(1, n):
        for f, c in [(i, 9 / 8), (i - 1, -9 / 8), (i + 1, -1 / 24), (i - 2, 1 / 24)]:
            rows.append(i - 1); cols.append(face_index(f)); vals.append(c / h)
    D = sp.csr_matrix((vals, (rows, cols)), shape=(n - 1, n), dtype=complex)
    Cn = (s + alpha) / (s + alpha + sig_node[1:n])
    Cf = (s + alpha) / (s + alpha + sig_face[:n])
    A = sp.diags(Cn) @ D @ sp.diags(Cf) @ G - mass * mass * sp.identity(n - 1, dtype=complex)
    return A.tocsr()


def profile(n, h, ab, width, smax):
    i = np.arange(n + 1)
    zn = np.clip((i * h - ab) / width, 0, 1)
    zf = np.clip(((i + .5) * h - ab) / width, 0, 1)
    return smax * zn ** 3, smax * zf ** 3


def reflection(k, h, mass, width=90., smax=20., alpha=.1):
    p = lattice_p(k, h); omega = np.hypot(mass, p); s = -1j * omega
    lead = max(200., 3.0 / k)            # physical lead length
    lead = np.ceil(lead / h) * h
    n = int(round((lead + width) / h)); ab = lead
    sn, sf = profile(n, h, ab, width, smax)
    A = build(n, h, sn, sf, s, alpha, mass)
    j0 = int(round(0.1 * lead / h))
    b = np.zeros(n - 1, dtype=complex); b[j0 - 1] = -1.
    U = spsolve((A + omega ** 2 * sp.identity(n - 1)).tocsc(), b)
    j = np.arange(int(round(0.25 * lead / h)), int(round(0.9 * lead / h)))
    x = j * h
    M = np.column_stack([np.exp(1j * k * x), np.exp(-1j * k * x)])
    coef, *_ = np.linalg.lstsq(M, U[j - 1], rcond=None)
    resid = np.linalg.norm(M @ coef - U[j - 1]) / np.linalg.norm(U[j - 1])
    return float(abs(coef[1] / coef[0]) ** 2), float(resid)


def bridge(binary, count=400, h=.25, mass=1.):
    raw = subprocess.check_output([binary, 'matrix', str(count), str(h), str(mass)], text=True)
    Mx = np.loadtxt(io.StringIO(raw), delimiter=',')
    n = count
    nU = n - 1
    out = []
    for k in [.05, .3, 1.]:
        p = lattice_p(k, h); omega = np.hypot(mass, p); s = -1j * omega
        sn, sf = profile(n, h, n * h / 2, n * h / 2, 20.)
        A = build(n, h, sn, sf, s, .1, mass)
        j0 = n // 8
        b = np.zeros(nU, dtype=complex); b[j0 - 1] = -1.
        U_mine = spsolve((A + omega ** 2 * sp.identity(nU)).tocsc(), b)
        F = np.zeros(len(Mx), dtype=complex); F[nU + j0 - 1] = 1.    # source in the V equation
        Y = np.linalg.solve(Mx + 1j * omega * np.eye(len(Mx)), -F)
        U_exe = Y[:nU]
        out.append(dict(k=k, h=h, nodes=n - 1, relative_difference=float(np.linalg.norm(U_exe - U_mine) / np.linalg.norm(U_mine))))
    return out


def main(outpath, binary=None):
    rows = []
    for h in [.5, .25, .125, .0625]:
        for mass in [1., np.sqrt(6.)]:
            for k in [.01, .02, .03, .05, .1, .2, .5, 1., 2., 3.]:
                R, res = reflection(k, h, mass)
                rows.append(dict(h=h, mass=float(mass), k=k, reflected_energy=R, lead_fit_residual=res))
                print(h, round(mass, 4), k, '%.4e' % R, '%.1e' % res, flush=True)
    soft = []
    for k in [.005, .01, .02, .03, .05, .1, .2]:
        R, res = reflection(k, 2., 1.)
        soft.append(dict(h=2., mass=1., k=k, reflected_energy=R, lead_fit_residual=res))
    out = dict(label='INDEPENDENT (own operator assembly, source-driven lead fit)', rows=rows, cold_grid_h2=soft,
               max_charged=max(r['reflected_energy'] for r in rows if r['mass'] == 1.),
               max_neutral=max(r['reflected_energy'] for r in rows if r['mass'] != 1.))
    if binary:
        out['bridge_to_executable_rhs'] = bridge(binary)
        print(out['bridge_to_executable_rhs'])
    json.dump(out, open(outpath, 'w'), indent=1)
    print('max charged', out['max_charged'], 'max neutral', out['max_neutral'])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
