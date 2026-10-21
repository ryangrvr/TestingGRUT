#!/usr/bin/env python3
"""Verifier v_wbound: planar kink-antikink collisions (phi + complex X, v = 1, lambda = 0.01, g^2 = 0.01,
lambda_X = 0.001), written independently of wbound/lattice.py.  Question: does the debris left between the receding
walls have a *locally time-averaged* w above 1/3 (not only an instantaneous excess), i.e. do the walls produce
w-bar > 1/3 states locally?

Setup: vacuum -1 outside, +1 between walls at z = -/+100, speeds +/-u; X seeded with white noise 1e-4.
KDK leapfrog, 2nd-order Laplacian, periodic box of 1200.  Window diagnostics: space-time integrals of rho and of
3p - rho = sum fdot^2 - sum f'^2 - 4U, and of Phi0 = lam (phi^2 - 1), over |z| < Zw for t in [t1, t2] chosen after
the receding walls have left the window.  Writes results/collide_1d.json.
"""
import json
import math
import os
import sys

import numpy as np

LAM, G2, LX = 0.01, 0.01, 0.001
DW = math.sqrt(2 / LAM)  # kink width 1/m with phi = tanh(z / DW)


def run(u, ppw=12, Lbox=1200.0, T=None, seed=7, windows=(30.0, 100.0)):
    g = 1 / math.sqrt(1 - u * u)
    dx = DW / g / ppw
    N = int(Lbox / dx)
    dt = 0.4 * dx
    z = (np.arange(N) - N / 2) * dx
    tcol = 100.0 / u
    if T is None:
        T = tcol + 440.0
    def kink(s):
        a = np.clip(g * s / DW, -300, 300)
        return np.tanh(a), (g / DW) / np.cosh(a) ** 2
    k1, d1 = kink(z + 100.0)
    k2, d2 = kink(z - 100.0)
    phi = k1 - k2 - 1.0
    pphi = -u * d1 - u * d2
    rng = np.random.default_rng(seed)
    x1 = 1e-4 * rng.normal(size=N)
    x2 = 1e-4 * rng.normal(size=N)
    p1 = np.zeros(N)
    p2 = np.zeros(N)

    def lap(f):
        return (np.roll(f, -1) - 2 * f + np.roll(f, 1)) / dx**2

    def forces(ph, a1, a2):
        X2 = 0.5 * (a1 * a1 + a2 * a2)
        fp = lap(ph) - LAM * ph * (ph * ph - 1) - 2 * G2 * ph * X2
        m = G2 * ph * ph + 2 * LX * X2
        return fp, lap(a1) - m * a1, lap(a2) - m * a2

    def dens(ph, a1, a2, q, q1, q2):
        X2 = 0.5 * (a1 * a1 + a2 * a2)
        U = LAM / 4 * (ph * ph - 1) ** 2 + G2 * ph * ph * X2 + LX * X2 * X2
        gr = sum(((np.roll(f, -1) - f) / dx) ** 2 for f in (ph, a1, a2))
        kin = q * q + q1 * q1 + q2 * q2
        return 0.5 * kin + 0.5 * gr + U, kin - gr - 4 * U

    nsteps = int(round(T / dt))
    every = max(1, int(round(0.5 / dt)))
    t1 = tcol + 150.0
    acc = {Zw: dict(rho=0.0, tr=0.0, phi0=0.0, n=0) for Zw in windows}
    masks = {Zw: np.abs(z) < Zw for Zw in windows}
    box = dict(rho=0.0, tr=0.0, n=0)
    maxphi_post, maxphi_win = 0.0, 0.0
    series = []
    fp, f1, f2 = forces(phi, x1, x2)
    for n in range(nsteps + 1):
        t = n * dt
        if n % every == 0:
            rho, tr = dens(phi, x1, x2, pphi, p1, p2)
            if t > tcol + 20:
                maxphi_post = max(maxphi_post, float(np.max(np.abs(phi))))
                box["rho"] += rho.sum() * dx
                box["tr"] += tr.sum() * dx
                box["n"] += 1
            if t >= t1:
                for Zw, mk in masks.items():
                    a = acc[Zw]
                    a["rho"] += rho[mk].sum() * dx
                    a["tr"] += tr[mk].sum() * dx
                    a["phi0"] += (LAM * (phi[mk] ** 2 - 1)).sum() * dx
                    a["n"] += 1
                maxphi_win = max(maxphi_win, float(np.max(np.abs(phi[masks[max(windows)]]))))
            if n % (every * 20) == 0:
                mk = masks[max(windows)]
                series.append(dict(t=t, E_box=float(rho.sum() * dx), E_win=float(rho[mk].sum() * dx),
                                   w_win_inst=float((tr[mk].sum() / rho[mk].sum() + 1) / 3) if rho[mk].sum() > 0 else None,
                                   maxphi_win=float(np.max(np.abs(phi[mk]))),
                                   X2max=float(np.max(0.5 * (x1 * x1 + x2 * x2)))))
        if n == nsteps:
            break
        pphi += 0.5 * dt * fp
        p1 += 0.5 * dt * f1
        p2 += 0.5 * dt * f2
        phi += dt * pphi
        x1 += dt * p1
        x2 += dt * p2
        fp, f1, f2 = forces(phi, x1, x2)
        pphi += 0.5 * dt * fp
        p1 += 0.5 * dt * f1
        p2 += 0.5 * dt * f2
    out = dict(u=u, gamma=g, dx=dx, N=N, dt=dt, T=T, t_collision=tcol, window_t_from=t1,
               max_absphi_after_collision=maxphi_post, max_absphi_in_window=maxphi_win,
               box_post_wbar=(box["tr"] / box["rho"] + 1) / 3)
    for Zw, a in acc.items():
        out[f"win{int(Zw)}"] = dict(E_mean=a["rho"] / a["n"], wbar_direct=(a["tr"] / a["rho"] + 1) / 3,
                                    wbar_from_Phi0=(a["phi0"] / a["rho"] + 1) / 3, n=a["n"])
    out["series"] = series
    return out


if __name__ == "__main__":
    us = [float(s) for s in (sys.argv[1:] or ["0.9", "0.99", "0.999"])]
    res = []
    for u in us:
        r = run(u)
        res.append(r)
        print(json.dumps({k: v for k, v in r.items() if k != "series"}), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "results", "collide_1d.json"), "w") as fh:
        json.dump(res, fh, indent=1)
