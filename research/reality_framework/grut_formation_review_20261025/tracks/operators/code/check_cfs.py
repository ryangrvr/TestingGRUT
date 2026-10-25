#!/usr/bin/env python3
"""Item 1 (and the lead/reflection part of item 2). Own code.
(a) resolvent identity: the zero-memory direct system reproduces exactly
    s^2 U - s U0 - V0 = C_n D C_f G U - m^2 U   (homogeneous linear vacuum response);
(b) cross-check of my assembled matrix against the RHS matrix exported by Codex's C++
    (REPRODUCTION, labelled; needs --codex-binary);
(c) finite-grid spectra in the production geometry (R=360, L=450 and R=720, L=810);
(d) reflection from an independent frequency-domain method (point source near r=0 in
    the actual radial grid, two-wave least-squares fit; no lead-root matching).
"""
import argparse, json, sys, subprocess, io, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
sys.path.insert(0, str(Path(__file__).parent))
from lattice import *

OUT = Path(__file__).resolve().parents[1] / 'results'


def resolvent_check():
    rng = np.random.default_rng(7)
    rows = []
    for (n, h, m, R, smax) in [(60, .5, 1.0, 18., 20.), (80, .25, np.sqrt(6), 12., 20.), (50, 1.0, 1.0, 30., 3.)]:
        A, ops = cfs_system(n, h, m, R, sigma_max=smax)
        u, v, p, q = ops['slices']
        U0, V0 = rng.normal(size=n - 1), rng.normal(size=n - 1)
        Y0 = np.zeros(A.shape[0]); Y0[u] = U0; Y0[v] = V0
        for s in [0.3 + 1.1j, -0.05 + 0.7j, 1e-3 + 2.5j, 0.8 - 0.2j]:
            Y = np.linalg.solve(s * np.eye(A.shape[0]) - A, Y0)
            M = stretched_operator(s, ops['G'], ops['D'], ops['sn'], ops['sf'], 0.1)
            Us = np.linalg.solve((s * s + m * m) * np.eye(n - 1) - M, s * U0 + V0)
            rows.append(dict(n=n, h=h, m=float(m), s=str(s), rel_diff=float(np.max(abs(Y[u] - Us)) / np.max(abs(Us)))))
    return rows


def codex_matrix_compare(binary):
    """REPRODUCTION: compare with Codex's exported C++ RHS matrix (setup: ab=n h/2)."""
    rows = []
    for count, h, mass in [(64, 1.0, 1.0), (128, .5, np.sqrt(6.))]:
        raw = subprocess.check_output([binary, 'matrix', str(count), str(h), repr(float(mass))], text=True)
        C = np.loadtxt(io.StringIO(raw), delimiter=',')
        n = count
        R = n * h / 2
        A, ops = cfs_system(n, h, mass, R, sigma_max=20.0)
        u, v, p, q = ops['slices']
        N = n - 1
        idx = list(range(0, N)) + list(range(N, 2 * N)) + [2 * N + f for f in range(n // 2, n)] + [2 * N + n + (j - 1) for j in range(n // 2 + 1, n)]
        Am = A[np.ix_(idx, idx)]
        rows.append(dict(count=count, h=h, mass=float(mass), shape=list(C.shape),
                         max_abs_diff=float(np.max(abs(Am - C))), max_abs_entry=float(np.max(abs(C)))))
    return rows


def spectra():
    rows = []
    for (L, R, h) in [(450., 360., .5), (450., 360., .25), (810., 720., .5)]:
        for m in (1.0, np.sqrt(6.)):
            n = int(round(L / h))
            A, ops = cfs_system(n, h, m, R)
            u, v, p, q = ops['slices']
            # drop memory states with sigma == 0 (they obey y' = -alpha y with zero data; exact restriction)
            keep_p = np.where(ops['sf'] > 0)[0] + p.start
            keep_q = np.where(ops['sn'] > 0)[0] + q.start
            idx = np.r_[np.arange(u.start, u.stop), np.arange(v.start, v.stop), keep_p, keep_q]
            Ar = A[np.ix_(idx, idx)]
            t0 = time.time()
            lam = np.linalg.eigvals(Ar)
            amp = {}
            for dt in (.03125, .015625):
                x = dt * lam
                amp[str(dt)] = float(np.max(abs(1 + x + x * x / 2 + x ** 3 / 6 + x ** 4 / 24)))
            phys = lam[np.argsort(-lam.real)[:5]]
            rows.append(dict(L=L, R=R, h=h, m=float(m), dim=int(len(idx)), max_real=float(lam.real.max()),
                             top5=[str(np.round(z, 8)) for z in phys], min_real=float(lam.real.min()),
                             RK4_max_amplification=amp, seconds=time.time() - t0))
            print('spectrum', rows[-1], flush=True)
    return rows


def sparse_ops(n, h):
    G = sp.csr_matrix(G_matrix_sparse(n, h))
    D = -G.T.tocsr()
    return G, D


def G_matrix_sparse(n, h):
    rows, cols, vals = [], [], []
    for f in range(n):
        for j, c in ((f + 1, C1), (f, -C1), (f + 2, C3), (f - 1, -C3)):
            idx, s = node_fold(j, n)
            if idx is not None:
                rows.append(f); cols.append(idx - 1); vals.append(s * c / h)
    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n - 1))


def reflection_fd(k, h, m, L=450., R=360., smax=20., alpha=.1, src_r=8., fit=(40., 300.)):
    n = int(round(L / h))
    G, D = sparse_ops(n, h)
    sn, sf = profiles(n, h, R, smax)
    p = lattice_momentum(k, h)
    omega = np.sqrt(m * m + p * p)
    s = -1j * omega
    Cn = (s + alpha) / (s + alpha + sn)
    Cf = (s + alpha) / (s + alpha + sf)
    M = sp.diags(Cn) @ D @ sp.diags(Cf) @ G
    A = ((s * s + m * m) * sp.eye(n - 1) - M).tocsc()
    r = np.arange(1, n) * h
    f = np.exp(-(r - src_r) ** 2 / 2.0)
    U = spl.spsolve(A, f)
    resid = np.max(abs(A @ U - f)) / np.max(abs(f))
    sel = (r >= fit[0]) & (r <= fit[1])
    j = np.arange(1, n)[sel]
    B = np.column_stack([np.exp(1j * k * h * j), np.exp(-1j * k * h * j)])
    coef, *_ = np.linalg.lstsq(B, U[sel], rcond=None)
    fitres = np.max(abs(B @ coef - U[sel])) / np.max(abs(U[sel]))
    Rr = coef[1] / coef[0]
    ideal = np.exp(-smax * (L - R) * p * omega / (alpha ** 2 + omega ** 2))
    return dict(k=float(k), h=h, m=float(m), L=L, R=R, reflected_energy=float(abs(Rr) ** 2),
                ideal_continuum_wall=float(ideal), solve_residual=float(resid), fit_residual=float(fitres))


def reflection_table():
    rows = []
    for h in (.5, .25, .125):
        for m in (1.0, np.sqrt(6.)):
            for k in (.002, .005, .01, .02, .03, .05, .1, 1.0, 3.0):
                if k * h > 2.5:
                    continue
                rows.append(reflection_fd(k, h, m))
                print('refl', rows[-1], flush=True)
    # zero-layer control: total reflection (|R|=1)
    ctrl = [reflection_fd(k, .25, 1.0, smax=0.0) for k in (.01, .1, 1.0)]
    return rows, ctrl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--codex-binary', default=None)
    ap.add_argument('--skip-spectra', action='store_true')
    a = ap.parse_args()
    out = {}
    out['resolvent_identity'] = resolvent_check()
    print('resolvent max rel diff', max(r['rel_diff'] for r in out['resolvent_identity']), flush=True)
    if a.codex_binary:
        out['REPRODUCTION_codex_cpp_matrix_vs_mine'] = codex_matrix_compare(a.codex_binary)
        print(out['REPRODUCTION_codex_cpp_matrix_vs_mine'], flush=True)
    out['reflection_fd_source_fit'], out['zero_layer_control'] = reflection_table()
    print('zero layer', out['zero_layer_control'], flush=True)
    if not a.skip_spectra:
        out['spectra'] = spectra()
    (OUT / 'cfs_checks.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
