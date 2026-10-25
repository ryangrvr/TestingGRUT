#!/usr/bin/env python3
"""Independent three-channel conversion rate  anti-X + object(Q) -> X + object(Q-2).

Claude review, track `conversion`.  Own code and method:
  * profile: my own BVP solve (profile.py), not the packet CSV;
  * radial system: the three-channel equation of bb211d1 MATHEMATICS.md sec. 2
    (quoted model: U = rho (u, v, sqrt2 w), K^2 = diag((nu+Om)^2-1, (nu-Om)^2-1,
    nu^2-6), V real symmetric) -- I re-derived the linearisation (NOTES.md);
  * integration: adaptive DOP853 of the three regular solutions together, with
    QR re-orthonormalisation every 1.5 length units (Godunov-type stabilisation
    against the growing closed s channel);
  * matching: Riccati-Bessel decomposition of the open channels and removal of
    the growing modified-Bessel component of the closed channel at rho_m;
    S matrix from flux-normalised in/out amplitudes.
Cross section sigma_{v->u} = pi (2l+1)/k_v^2 |S_uv|^2 and rate sigma*k_v/e_v.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn, spherical_yn, spherical_in, spherical_kn

A, B = 0.1, 3.0
HBARC_GEV_CM = 1.973269804e-14
C_CM_S = 2.99792458e10
M_X = 0.0924089637375362
UNIT_SV = (HBARC_GEV_CM / M_X) ** 2 * C_CM_S       # cm^3/s per m_X^-2 (times v in c)
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/conversion')

d = np.load(SCR / 'my_profile_Q200000.npz')
OM = float(d['Omega'])
_r = d['r']
_F = CubicSpline(_r, d['F'])
_S = CubicSpline(_r, d['S'])


def potential(rho):
    F = _F(rho); S = _S(rho)
    v11 = S * S - 1 + 2 * A * F * F
    v12 = A * F * F
    v13 = np.sqrt(2) * S * F
    v33 = B * (3 * S * S - 1) + F * F - 6
    return np.array([[v11, v12, v13], [v12, v11, v13], [v13, v13, v33]])


def regular_start(ell, K2, rho0):
    M0 = potential(0.0) - np.diag(K2)
    U = np.eye(3) + M0 * rho0 ** 2 / (2 * (2 * ell + 3))
    Up = ((ell + 1) * np.eye(3) + (ell + 3) * M0 * rho0 ** 2 / (2 * (2 * ell + 3))) / rho0
    return U, Up   # common factor rho0^ell(+1) dropped


def propagate(ell, K2, rho_m, rtol=1e-11, seg=1.5):
    rho0 = 1e-3
    U, Up = regular_start(ell, K2, rho0)
    Y = np.vstack([U, Up])               # 6 x 3
    cent = ell * (ell + 1)
    Kd = np.diag(K2)

    def f(rho, y):
        Ym = y.reshape(6, 3)
        W = cent / rho ** 2 * np.eye(3) + potential(rho) - Kd
        return np.vstack([Ym[3:], W @ Ym[:3]]).ravel()

    edges = np.r_[rho0, np.arange(seg, rho_m, seg), rho_m]
    edges = np.unique(edges)
    for a, b in zip(edges[:-1], edges[1:]):
        s = solve_ivp(f, (a, b), Y.ravel(), method='DOP853', rtol=rtol, atol=1e-300)
        Y = s.y[:, -1].reshape(6, 3)
        q, _ = np.linalg.qr(Y)
        Y = q
    return Y


def smatrix(ell, e_anti, rho_m=40.0, rtol=1e-11):
    nu = e_anti + OM
    K2 = np.array([(nu + OM) ** 2 - 1, (nu - OM) ** 2 - 1, nu ** 2 - 6])
    assert K2[0] > 0 and K2[1] > 0 and K2[2] < 0
    Y = propagate(ell, K2, rho_m, rtol)
    U, Up = Y[:3], Y[3:]
    k = np.sqrt(K2[:2]); kap = np.sqrt(-K2[2])
    # Open channels: U = A jhat + Bc nhat (Riccati-Bessel, x = k rho)
    Amat = np.zeros((2, 3)); Bmat = np.zeros((2, 3))
    for i in range(2):
        x = k[i] * rho_m
        jh = x * spherical_jn(ell, x); jhp = k[i] * (spherical_jn(ell, x) + x * spherical_jn(ell, x, True))
        nh = x * spherical_yn(ell, x); nhp = k[i] * (spherical_yn(ell, x) + x * spherical_yn(ell, x, True))
        W = jh * nhp - jhp * nh
        Amat[i] = (U[i] * nhp - Up[i] * nh) / W
        Bmat[i] = (jh * Up[i] - jhp * U[i]) / W
    # Closed channel: U = g khat + dlt ihat; remove growing ihat component.
    x = kap * rho_m
    ih = x * spherical_in(ell, x); ihp = kap * (spherical_in(ell, x) + x * spherical_in(ell, x, True))
    kh = x * spherical_kn(ell, x); khp = kap * (spherical_kn(ell, x) + x * spherical_kn(ell, x, True))
    Wc = kh * ihp - khp * ih
    dlt = (kh * Up[2] - khp * U[2]) / Wc
    # null space of the 1x3 row dlt -> 3x2 combination matrix
    _, _, vt = np.linalg.svd(dlt.reshape(1, 3))
    N = vt[1:].T
    Aa = Amat @ N; Bb = Bmat @ N             # 2x2: channel x solution
    out_amp = (-1j * Aa - Bb) / 2 * np.sqrt(k)[:, None]
    in_amp = (1j * Aa - Bb) / 2 * np.sqrt(k)[:, None]
    S = out_amp @ np.linalg.inv(in_amp)
    unit = np.max(np.abs(S.conj().T @ S - np.eye(2)))
    recip = abs(S[0, 1] - S[1, 0])
    return S, dict(unitarity=float(unit), reciprocity=float(recip))


def rate(p_anti, ellmax=None, rho_m=40.0, rtol=1e-11, tol_tail=1e-9, verbose=False):
    e = np.sqrt(1 + p_anti ** 2)
    total = 0.0; parts = []
    ell = 0
    while True:
        S, chk = smatrix(ell, e, rho_m, rtol)
        P = abs(S[0, 1]) ** 2
        sig = np.pi * (2 * ell + 1) / p_anti ** 2 * P
        parts.append(dict(ell=ell, P_conversion=float(P), sigma_v_cm3_s=float(sig * p_anti / e * UNIT_SV), **chk))
        total += sig
        if verbose:
            print(ell, P, chk, flush=True)
        ell += 1
        if ellmax is not None and ell > ellmax:
            break
        if ellmax is None and ell >= 4 and (parts[-1]['sigma_v_cm3_s'] + parts[-2]['sigma_v_cm3_s']) < tol_tail * total * p_anti / e * UNIT_SV:
            break
        if ell > 25:
            break
    return float(total * p_anti / e * UNIT_SV), parts


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'spot'
    res = Path(__file__).resolve().parents[1] / 'results'
    if mode == 'spot':
        out = []
        for p in [1e-5, 0.1, 0.5047378541243651, 0.7071067811865476]:
            t = time.time()
            sv, parts = rate(p)
            ctrl = {}
            # numerical controls at the two physically important momenta
            if p in (1e-5, 0.7071067811865476):
                ctrl['rho_m_30'] = rate(p, ellmax=len(parts) - 1, rho_m=30.0)[0]
                ctrl['rtol_1e-9'] = rate(p, ellmax=len(parts) - 1, rtol=1e-9)[0]
            out.append(dict(p_over_m=p, sigma_v_cm3_s=sv, controls=ctrl, partials=parts, seconds=time.time() - t))
            print(p, sv, ctrl, time.time() - t, flush=True)
        (res / 'scatter_spot.json').write_text(json.dumps(dict(Omega=OM, unit_cm3_s=UNIT_SV, rows=out), indent=1) + '\n')
    elif mode == 'scan':
        lo, hi, n = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
        tag = sys.argv[5]
        ps = np.linspace(lo, hi, n)
        rows = []
        for p in ps:
            sv, parts = rate(p, rtol=1e-10)
            rows.append(dict(p_over_m=float(p), sigma_v_cm3_s=sv,
                             max_unitarity=max(q['unitarity'] for q in parts),
                             ellmax=len(parts) - 1, P_by_ell=[q['P_conversion'] for q in parts]))
            print(p, sv, flush=True)
        (SCR / f'scan_{tag}.json').write_text(json.dumps(rows) + '\n')


if __name__ == '__main__':
    main()
