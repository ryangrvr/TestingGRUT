#!/usr/bin/env python3
"""cf_scale_skeptic: P4-S fixed points by total-degree homotopy continuation (own code).

cf_scale counted real fixed points of the six hidden quartics (lam_X, lam_Y, lam_S, l_XY, l_XS, l_YS) at e = e* with a
700-start multistart (which cannot certify completeness). Here every isolated complex solution is tracked from the
2^6 = 64 Bezout start points of G_i = z_i^2 - 1 with the gamma trick, H(z,t) = (1-t) gamma G(z) + t F(z), Euler-
predictor / Newton-corrector with adaptive steps. Finite endpoints are polished with Newton at t = 1 and deduplicated.
Polynomials: the one-loop betas (checked in mv_beta_check.py) at kappa = 0, phi sector = 0, plus 16 pi^2 A_lam lam_i.

Analytic statement checked alongside [D]: the lam_S equation reads 20 L^2 + (16pi^2 A - 432 e^2) L + C = 0 with
C = 7776 e^4 + l_XS^2 + l_YS^2 >= 7776 e^4 > 0 for real portals. For A >= 0 a real root needs
16pi^2 A - 432 e^2 >= sqrt(80*7776) e^2, and then both roots are negative. So no real fixed point with lam_S > 0 exists
for any A_lam >= 0 and any e* > 0; V4 is then unbounded below along pure Sigma.
Usage: python3 -B p4s_homotopy.py <out.json>
"""
import json, sys, itertools
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
K = 16 * np.pi**2
E_STAR = 0.31818              # cf_scale e*_hidden at f_g = 9.83e-3 (recomputed below)


def F(z, A, e):
    lX, lY, lS, xy, xs, ys = z
    e2, e4 = e * e, e**4
    return np.array([
        20 * lX**2 - 12 * e2 * lX + 6 * e4 + xs**2 + xy**2 + K * A * lX,
        20 * lY**2 - 108 * e2 * lY + 486 * e4 + xy**2 + ys**2 + K * A * lY,
        20 * lS**2 - 432 * e2 * lS + 7776 * e4 + xs**2 + ys**2 + K * A * lS,
        4 * xy**2 + 8 * xy * (lX + lY) - 60 * e2 * xy + 108 * e4 + 2 * xs * ys + K * A * xy,
        4 * xs**2 + 8 * xs * (lS + lX) - 222 * e2 * xs + 432 * e4 + 2 * xy * ys + K * A * xs,
        4 * ys**2 + 8 * ys * (lS + lY) - 270 * e2 * ys + 3888 * e4 + 2 * xs * xy + K * A * ys])


def JF(z, A, e):
    lX, lY, lS, xy, xs, ys = z
    e2 = e * e
    return np.array([
        [40 * lX - 12 * e2 + K * A, 0, 0, 2 * xy, 2 * xs, 0],
        [0, 40 * lY - 108 * e2 + K * A, 0, 2 * xy, 0, 2 * ys],
        [0, 0, 40 * lS - 432 * e2 + K * A, 0, 2 * xs, 2 * ys],
        [8 * xy, 8 * xy, 0, 8 * xy + 8 * (lX + lY) - 60 * e2 + K * A, 2 * ys, 2 * xs],
        [8 * xs, 0, 8 * xs, 2 * ys, 8 * xs + 8 * (lS + lX) - 222 * e2 + K * A, 2 * xy],
        [0, 8 * ys, 8 * ys, 2 * xs, 2 * xy, 8 * ys + 8 * (lS + lY) - 270 * e2 + K * A]], dtype=complex)


def track(z0, A, e, gamma):
    G = lambda z: z**2 - 1
    JG = lambda z: np.diag(2 * z)
    H = lambda z, t: (1 - t) * gamma * G(z) + t * F(z, A, e)
    Hz = lambda z, t: (1 - t) * gamma * JG(z) + t * JF(z, A, e)
    Ht = lambda z, t: -gamma * G(z) + F(z, A, e)
    z = z0.astype(complex); t = 0.0; dt = 1e-3
    while t < 1.0:
        dt = min(dt, 1.0 - t)
        try:
            dz = -np.linalg.solve(Hz(z, t), Ht(z, t))
        except np.linalg.LinAlgError:
            return None, 'singular'
        zp = z + dt * dz; tn = t + dt
        ok = False
        for _ in range(6):
            try:
                step = np.linalg.solve(Hz(zp, tn), H(zp, tn))
            except np.linalg.LinAlgError:
                break
            zp = zp - step
            if np.linalg.norm(step) < 1e-11 * (1 + np.linalg.norm(zp)):
                ok = True
                break
        if ok and np.linalg.norm(zp - z) < 0.2 * (1 + np.linalg.norm(z)):
            z, t = zp, tn
            dt = min(dt * 1.6, 0.05)
        else:
            dt *= 0.5
            if dt < 1e-14:
                return z, 'stalled'
        if np.linalg.norm(z) > 1e8:
            return z, 'infinity'
    for _ in range(20):
        try:
            z = z - np.linalg.solve(JF(z, A, e), F(z, A, e))
        except np.linalg.LinAlgError:
            break
    return z, 'finite'


def V4min_hidden(z, rng, n=200):
    """min of V4 on the unit sphere of the six charged real components, by sampling magnitudes (V4 depends only on
    |X|^2, |Y|^2, |S|^2 when kappa = 0): min over the simplex of a quadratic form in (a_X, a_Y, a_S)."""
    lX, lY, lS, xy, xs, ys = z
    Mq = np.array([[lX, xy / 2, xs / 2], [xy / 2, lY, ys / 2], [xs / 2, ys / 2, lS]])
    W = rng.dirichlet(np.ones(3), size=20000)
    W = np.vstack([W, np.eye(3), [[.5, .5, 0], [.5, 0, .5], [0, .5, .5]]])
    vals = np.einsum('ni,ij,nj->n', W, Mq, W)
    return float(vals.min())


def main(out):
    rng = np.random.default_rng(7)
    fg = (0.6154 * np.sqrt(3 / 5))**2 * (41 / 6) / K
    e = np.sqrt(K * fg / (46 / 3))
    res = dict(method=__doc__, f_g=fg, e_star=e, A_crit_analytic=(432 + np.sqrt(80 * 7776)) * e**2 / K, scan=[])
    gamma = np.exp(1j * 2.1374)
    starts = [np.array(s, float) for s in itertools.product((1.0, -1.0), repeat=6)]
    for A in (0.0237, 0.1, 0.5, 0.78, 0.79, 1.0, 2.83, 10.0):
        sols, status = [], {}
        for s0 in starts:
            z, st = track(s0, A, e, gamma)
            status[st] = status.get(st, 0) + 1
            if st == 'finite' and np.max(np.abs(F(z, A, e))) < 1e-9 * (1 + np.max(np.abs(z))**2):
                if not any(np.allclose(z, w, rtol=1e-7, atol=1e-10) for w in sols):
                    sols.append(z)
        real = [z.real for z in sols if np.max(np.abs(z.imag)) < 1e-8 * (1 + np.max(np.abs(z)))]
        rows = []
        for z in real:
            rows.append(dict(z=dict(zip(['lam_X', 'lam_Y', 'lam_S', 'l_XY', 'l_XS', 'l_YS'], z.tolist())),
                             V4min_simplex=V4min_hidden(z, rng)))
        res['scan'].append(dict(A_lam=A, path_status=status, n_finite_distinct=len(sols), n_real=len(real),
                                n_real_with_lamS_pos=sum(r['z']['lam_S'] > 0 for r in rows),
                                n_bounded_below=sum(r['V4min_simplex'] > 0 for r in rows),
                                max_lamS_real=max([r['z']['lam_S'] for r in rows]) if rows else None, real_fps=rows))
        print('A=%.4g paths %s | finite distinct %d | real %d | lam_S>0: %d | bounded below: %d | max lam_S %s' % (
            A, status, len(sols), len(real), res['scan'][-1]['n_real_with_lamS_pos'], res['scan'][-1]['n_bounded_below'],
            res['scan'][-1]['max_lamS_real']), flush=True)
    print('e* = %.5f  A_crit = %.4f' % (e, res['A_crit_analytic']))
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main(sys.argv[1])
