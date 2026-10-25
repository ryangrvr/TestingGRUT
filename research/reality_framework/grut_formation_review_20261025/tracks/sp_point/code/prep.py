"""sp_point: which preparation is physical? (own code)

a) Conditional radial re-minimisation at the displaced phase (own global radial solver; s -> h/sqrt2 exactly at
   delta = pi/2) and its release; sensitivity to Codex's slightly unconverged s = 0.0063102554.
b) Conditional preparation across delta in (0, pi).
c) Interpolation r(f) = r_aligned + f (r_cond - r_aligned) at delta = pi/2.
d) Inflationary relaxation: start from Codex's supplied state inside de Sitter (H = 1, -c_H|X|^2 on), evolve N e-folds
   with the full potential, then release (velocities carried over). Shows how late the preparation must act.
e) Energetics / de Sitter statistics of the supplied displacement: Delta V, Mahalanobis distance in the Gaussian
   Starobinsky-Yokoyama equilibrium around the aligned minimum, light-mode kick per e-fold.

Usage: python3 prep.py ../results/prep.json
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import solve_ivp
import spmodel as m

H = m.H_STAR
C = 0.025


def rel(z, v=None, tang=False):
    r = m.run(H, np.asarray(z, complex), a_end=1000, cpar=C, windows=(1, 300, 600, 1000), ns=8, v0=v)
    w = r['windows']['600-1000']
    return dict(eps=w['epsX'], NX=w['NX'], EX=w['EX'], JY=w['JY'], YX_over_YQ=w['YX_over_YQ'],
                eps300=r['windows']['300-600']['epsX'], pairs=w['pairs_per_charge'])


def F_infl(z, cH=1.0):
    M, h, mu = m.pars(H)
    x, y, s = z
    return np.array([2 * m.LX * abs(x)**2 * x + 3 * m.KAP * y * np.conj(x)**2 - cH * x,
                     M * M * y + 2 * m.LY * abs(y)**2 * y + m.KAP * x**3 + 2 * mu * s * np.conj(y),
                     2 * m.LS * (abs(s)**2 - h * h / 2) * s + mu * y * y])


def desitter(z0, N):
    def rhs(t, u):
        z = u[0:6:2] + 1j * u[1:6:2]
        v = u[6::2] + 1j * u[7::2]
        a = -3 * v - F_infl(z)
        out = np.empty(12)
        out[0:6:2], out[1:6:2] = v.real, v.imag
        out[6::2], out[7::2] = a.real, a.imag
        return out
    u0 = np.zeros(12)
    u0[0:6:2], u0[1:6:2] = z0.real, z0.imag
    if N == 0:
        return z0, np.zeros(3, complex)
    sol = solve_ivp(rhs, [0, N], u0, method='DOP853', rtol=1e-12, atol=1e-14)
    u = sol.y[:, -1]
    return u[0:6:2] + 1j * u[1:6:2], u[6::2] + 1j * u[7::2]


def invariants(z):
    x, y, s = z
    return dict(theta_kappa=float(np.angle(np.conj(y) * x**3)), theta_mu=float(np.angle(np.conj(s) * y * y)),
                absX=float(abs(x)), absY=float(abs(y)), absS=float(abs(s)))


def main(out):
    res = {}
    t0 = time.time()
    xa, ya, sa = m.minimum(H, 1.0)
    z_sup = np.array([xa, -ya, -sa * np.exp(1j * np.pi / 2)], complex)
    xc, yc, sc = m.minimum(H, 1.0, cos_mu=-np.cos(np.pi / 2))
    z_cond = np.array([xc, -yc, -sc * np.exp(1j * np.pi / 2)], complex)
    z_cond_codex = np.array([23.823180454537532, -1.8901458070297728, -0.0063102553808654935 * 1j], complex)
    with ProcessPoolExecutor(max_workers=4) as pool:
        A = list(pool.map(rel, [z_sup, z_cond, z_cond_codex]))
        res['a'] = dict(aligned_radii=[xa, ya, sa], conditional_radii=[xc, yc, sc],
                        supplied=A[0], conditional=A[1], conditional_with_codex_radii=A[2],
                        codex_conditional_abs_eps=0.0516)
        print('a', res['a'], time.time() - t0, flush=True)
        ds = np.linspace(0.05, np.pi - 0.05, 61)
        zs = []
        for d in ds:
            x, y, s = m.minimum(H, 1.0, cos_mu=-np.cos(d))
            zs.append(np.array([x, -y, -s * np.exp(1j * d)], complex))
        B = list(pool.map(rel, zs))
        res['b'] = [dict(delta=float(d), radii=[float(abs(z[0])), float(abs(z[1])), float(abs(z[2]))], **b)
                    for d, z, b in zip(ds, zs, B)]
        print('b max|eps|', max(abs(b['eps']) for b in B), time.time() - t0, flush=True)
        fs = np.linspace(0, 1, 21)
        zs = [z_sup + f * (z_cond - z_sup) for f in fs]
        Cc = list(pool.map(rel, zs))
        res['c'] = [dict(f=float(f), **c) for f, c in zip(fs, Cc)]
        print('c', [(round(f, 2), round(c['eps'], 4)) for f, c in zip(fs, Cc)], time.time() - t0, flush=True)
        Ns = [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0]
        states = [desitter(z_sup, N) for N in Ns]
        D = list(pool.map(rel, [s[0] for s in states], [s[1] for s in states]))
        res['d'] = [dict(N_efolds=N, invariants=invariants(s[0]), absv=[float(abs(v)) for v in s[1]], **d)
                    for N, s, d in zip(Ns, states, D)]
        print('d', [(N, round(d['eps'], 4), round(s[0][2].real, 4), round(s[0][2].imag, 4)) for N, s, d in zip(Ns, states, D)],
              time.time() - t0, flush=True)
    # e) energetics and de Sitter statistics
    z_al = np.array([xa, -ya, -sa], complex)
    dV_sup = m.V(z_sup, H, 1.0) - m.V(z_al, H, 1.0)
    dV_cond = m.V(z_cond, H, 1.0) - m.V(z_al, H, 1.0)
    def Vr(u):
        return m.V(z_al + (u[0::2] + 1j * u[1::2]), H, 1.0)
    e = 1e-4
    Hs = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            di = np.zeros(6); dj = np.zeros(6); di[i] = e; dj[j] = e
            Hs[i, j] = (Vr(di + dj) - Vr(di - dj) - Vr(-di + dj) + Vr(-di - dj)) / (4 * e * e)
    Hp = Hs / 2  # canonical phi = sqrt2 (Re, Im)
    w, U = np.linalg.eigh(Hp)
    keep = w > 1e-8 * w.max()
    def mahal(z):
        du = z - z_al
        u = np.ravel([[q.real, q.imag] for q in du])
        phi = np.sqrt(2) * u
        c = U.T @ phi
        return float(np.sqrt(np.sum(c[keep]**2 * w[keep] / (3 / (8 * np.pi**2)))))
    M_, h_, mu_ = m.pars(H)
    res['e'] = dict(dV_supplied_over_H4=dV_sup, dV_conditional_over_H4=dV_cond,
                    twice_mu_S_Y2=2 * mu_ * sa * ya * ya, dS_scale_3_over_8pi2=3 / (8 * np.pi**2),
                    hessian_mass2_over_H2=w.tolist(), mahalanobis_supplied=mahal(z_sup),
                    mahalanobis_conditional=mahal(z_cond),
                    gaussian_density_ratio_supplied=float(np.exp(-mahal(z_sup)**2 / 2)),
                    sigma_phase_kick_per_efold=1 / (2 * np.pi * np.sqrt(2) * sa),
                    sigma_transverse_rms=float(np.sqrt(3 / (8 * np.pi**2) / (2 * 2 * m.LS * (sa * sa - h_ * h_ / 2) / 2)) / np.sqrt(2)),
                    masses=m.masses(H, 1.0),
                    note='Gaussian (quadratic) de Sitter equilibrium around the aligned minimum; estimate only.')
    print('e', res['e'], flush=True)
    json.dump(dict(H_GeV=H, results=res, integrator='Yoshida-8 c=%g; de Sitter stage scipy DOP853 rtol 1e-12' % C),
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
