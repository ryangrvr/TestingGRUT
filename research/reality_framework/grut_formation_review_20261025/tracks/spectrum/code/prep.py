#!/usr/bin/env python3
"""Initial data for kg_radial.c (spectrum track, Claude, 25 Oct 2026).

prepared : psi = 0.18 exp(-r^2/7200), psi_t = -i psi, S_t = 0, S the static neutral solution
           computed by Newton iteration on the solver's own 6th-order grid operator (so the neutral
           field starts exactly in discrete equilibrium with |psi|^2).
packet   : linear-mode validation data, an outgoing particle packet plus an outgoing antiparticle
           packet (exact discrete positive/negative-frequency projections in the sine basis), with
           exactly known N_X, N_Xbar, E and a neutral packet.
Binary layout: int64 N, float64 h, then x,y,z,vx,vy,vz each of length N+1 (U=r psi, Z=r(S-1)).
"""
import argparse, json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.fft import dst


def lap6_matrix(N, h):
    """6th-order Laplacian on interior nodes 1..N-1 with odd ghosts at 0 and N."""
    n = N - 1
    c = {0: -49 / 18, 1: 1.5, 2: -3 / 20, 3: 1 / 90}
    rows, cols, vals = [], [], []
    for i in range(1, N):
        for off in range(-3, 4):
            j = i + off
            w = c[abs(off)] / h**2
            if j < 0:
                j, w = -j, -w
            elif j > N:
                j, w = 2 * N - j, -w
            if j == 0 or j == N:
                continue
            rows.append(i - 1); cols.append(j - 1); vals.append(w)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def static_Z(r, f, L6):
    """Solve L6 Z = r S (3(S^2-1) + f^2), S = 1 + Z/r, by Newton."""
    Z = -r * f**2 / 6.0
    for it in range(50):
        S = 1 + Z / r
        F = L6 @ Z - r * S * (3 * (S * S - 1) + f * f)
        dF = 3 * (S * S - 1) + f * f + 6 * S * S  # d/dZ of r S(...) = (3(S^2-1)+f^2) + S*6S
        J = L6 - sp.diags(dF)
        dZ = spl.spsolve(J.tocsc(), -F)
        Z += dZ
        if np.max(np.abs(dZ)) < 1e-15:
            break
    S = 1 + Z / r
    res = np.max(np.abs(L6 @ Z - r * S * (3 * (S * S - 1) + f * f)))
    return Z, res, it


def write(path, N, h, arrs):
    with open(path, 'wb') as fo:
        np.array([N], dtype=np.int64).tofile(fo)
        np.array([h], dtype=np.float64).tofile(fo)
        for a in arrs:
            full = np.zeros(N + 1); full[1:N] = a
            full.tofile(fo)


def prepared(L, h, path):
    N = int(round(L / h)); r = np.arange(1, N) * h
    f = 0.18 * np.exp(-r * r / 7200.0)
    L6 = lap6_matrix(N, h)
    Z, res, it = static_Z(r, f, L6)
    U = r * f; V = -1j * U
    write(path, N, h, [U.real, U.imag, Z, V.real, V.imag, np.zeros_like(r)])
    W = 4 * np.pi * h
    Q = W * np.sum(-np.imag(np.conj(U) * V))
    return dict(N=N, h=h, L=L, static_newton_residual=float(res), newton_iterations=it, Q0=float(Q))


def packet(L, h, path, spec):
    """Outgoing particle + antiparticle + neutral packets in linear mode with exact mode numbers."""
    N = int(round(L / h)); r = np.arange(1, N) * h
    n = np.arange(1, N); k = n * np.pi / L
    K2 = (49 / 18 - 3 * np.cos(k * h) + 0.3 * np.cos(2 * k * h) - (1 / 45) * np.cos(3 * k * h)) / h**2
    w1 = np.sqrt(1 + K2); w6 = np.sqrt(6 + K2)
    T = lambda a: dst(a, type=1, norm='ortho')  # orthonormal, self-inverse
    U = np.zeros(N - 1, complex); V = np.zeros(N - 1, complex)
    for (amp, c, wdt, k0, sign) in spec['charged']:
        g = amp * np.exp(-(r - c)**2 / (2 * wdt**2))
        u = g * np.exp(1j * sign * k0 * r)          # sign=+1 particle (outgoing e^{-iwt+ikr}); -1 antiparticle
        um = T(u.real) + 1j * T(u.imag)
        vm = -1j * sign * w1 * um                    # exact +/- frequency projection
        U += T(um.real) + 1j * T(um.imag); V += T(vm.real) + 1j * T(vm.imag)
    Zr = np.zeros(N - 1); Zv = np.zeros(N - 1)
    for (amp, c, wdt, k0) in spec['neutral']:
        g = amp * np.exp(-(r - c)**2 / (2 * wdt**2))
        u = g * np.exp(1j * k0 * r); um = T(u.real) + 1j * T(u.imag); vm = -1j * w6 * um
        Zr += np.real(T(um.real) + 1j * T(um.imag)); Zv += np.real(T(vm.real) + 1j * T(vm.imag))
    write(path, N, h, [U.real, U.imag, Zr, V.real, V.imag, Zv])
    um = T(U.real) + 1j * T(U.imag); vm = T(V.real) + 1j * T(V.imag)
    a = 0.5 * (um + 1j * vm / w1); b = 0.5 * (um - 1j * vm / w1)
    W = 4 * np.pi * h
    NX = W * np.sum(w1 * abs(a)**2); NXb = W * np.sum(w1 * abs(b)**2)
    zm = T(Zr); zv = T(Zv); Es = W * np.sum(0.5 * (zv**2 + w6**2 * zm**2)); Ns = W * np.sum(0.5 * (zv**2 + w6**2 * zm**2) / w6)
    return dict(N=N, h=h, L=L, NX=float(NX), NXbar=float(NXb), E_charged=float(W * np.sum(w1**2 * (abs(a)**2 + abs(b)**2))),
                Q=float(NX - NXb), E_neutral=float(Es), N_neutral=float(Ns), spec=spec)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('kind', choices=['prepared', 'packet'])
    ap.add_argument('--L', type=float, required=True)
    ap.add_argument('--h', type=float, required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    if a.kind == 'prepared':
        info = prepared(a.L, a.h, a.out)
    else:
        spec = dict(charged=[(0.02, 150., 25., 0.40, +1), (0.01, 150., 40., 0.15, +1), (0.004, 170., 25., 0.55, -1),
                             (0.01, 160., 60., 0.06, +1)],
                    neutral=[(0.003, 150., 20., 0.8)])
        info = packet(a.L, a.h, a.out, spec)
    print(json.dumps(info))
    with open(a.out + '.json', 'w') as fo:
        json.dump(info, fo, indent=1)
