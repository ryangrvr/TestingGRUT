#!/usr/bin/env python3
"""Verifier v_wbound: spherically symmetric 3D kinetic lumps (phi only, v = 1, lambda = 0.01), flat space.
Independent test of the wbound dilatation theorem  int_0^T (3P - E) dt = D(T) - D(0) - E T <= E R0 + |D(0)|
and of how long a coherent over-barrier region keeps wbar_T > 1/3, at radii (R0 up to 400) that a 3D Cartesian
lattice cannot reach.  Spherical symmetry suppresses non-radial instabilities (and there is no X here), so these runs
are the *most favourable* case for a long-lived excess.

Discretization: cell centres r_i = (i+1/2) dr, shell weights V_i = r_i^2 dr, face areas A_{i+1/2} = r_{i+1/2}^2;
E = sum V_i (pi^2/2 + U) + sum A dr (dphi/dr)^2/2 is conserved by the semi-discrete flow; leapfrog (KDK).
3P - E = sum V_i (pi^2 - 4U) - sum A dr (dphi/dr)^2 ; D = -sum r pi dphi/dr (face-averaged).  Outer edge beyond
R0 + T so nothing reflects.  Writes results/radial_lump.json.
"""
import json
import math
import os
import sys

import numpy as np

LAM = 0.01


def run(R0, A, profile="cos2", dr=0.2, dt=0.05, T=None, sample=1.0):
    if T is None:
        T = 3 * R0 + 150
    Rmax = R0 + T + 60
    N = int(Rmax / dr)
    r = (np.arange(N) + 0.5) * dr
    rf = (np.arange(N - 1) + 1.0) * dr  # interior faces
    Vi = r**2 * dr
    Af = rf**2
    phi = np.ones(N)
    if profile == "cos2":
        pi = np.where(r < R0, A * np.cos(0.5 * math.pi * r / R0) ** 2, 0.0)
    else:  # flat top with cos^2 edge of width 8
        e = 8.0
        pi = np.where(r < R0 - e, A, np.where(r < R0, A * np.cos(0.5 * math.pi * (r - (R0 - e)) / e) ** 2, 0.0))

    def force(ph):
        flux = Af * (ph[1:] - ph[:-1]) / dr
        lap = np.zeros(N)
        lap[:-1] += flux
        lap[1:] -= flux
        return lap / Vi - LAM * ph * (ph * ph - 1)

    def diag(ph, p):
        g = (ph[1:] - ph[:-1]) / dr
        U = LAM / 4 * (ph * ph - 1) ** 2
        E = np.sum(Vi * (0.5 * p * p + U)) + np.sum(Af * dr * 0.5 * g * g)
        tr = np.sum(Vi * (p * p - 4 * U)) - np.sum(Af * dr * g * g)
        pf = 0.5 * (p[1:] + p[:-1])
        Dm = -np.sum(Af * dr * rf * pf * g)
        return E, tr, Dm

    nsteps = int(round(T / dt))
    every = max(1, int(round(sample / dt)))
    ts, Es, trs, Ds, maxphi = [], [], [], [], []
    F = force(phi)
    for n in range(nsteps + 1):
        if n % every == 0:
            E, tr, Dm = diag(phi, pi)
            ts.append(n * dt)
            Es.append(E)
            trs.append(tr)
            Ds.append(Dm)
            maxphi.append(float(np.max(np.abs(phi))))
        if n == nsteps:
            break
        pi += 0.5 * dt * F
        phi += dt * pi
        F = force(phi)
        pi += 0.5 * dt * F
    ts, Es, trs, Ds = map(np.array, (ts, Es, trs, Ds))
    E0 = Es[0]
    # running integral W(T) = int_0^T (3P - E) dt  (trapezoid)
    W = np.concatenate([[0.0], np.cumsum(0.5 * (trs[1:] + trs[:-1]) * np.diff(ts))])
    wbarT = np.where(ts > 0, (W / np.maximum(ts, 1e-30) + E0) / (3 * E0), np.nan)  # wbar_T = <P>/E
    above = np.where(wbarT[1:] > 1 / 3)[0]
    last_above = float(ts[1:][above[-1]]) if len(above) else 0.0
    # dilatation identity: D(T) - D(0) = int 3P dt = W + E T
    dil_res = np.max(np.abs((Ds - Ds[0]) - (W + E0 * ts))) / (E0 * ts[-1])
    bound = E0 * R0 + abs(Ds[0])
    res = dict(R0=R0, A=A, profile=profile, dr=dr, dt=dt, T=T, E0=E0, E_drift=float(np.ptp(Es) / E0),
               local_eps_center_over_lam=0.5 * A * A / LAM, D0=float(Ds[0]),
               max_W_over_E_R0=float(np.max(W) / (E0 * R0)), W_end_over_E_R0=float(W[-1] / (E0 * R0)),
               bound_violated=bool(np.max(W) > bound * (1 + 1e-3)), dilatation_identity_rel=float(dil_res),
               wbarT_max=float(np.nanmax(wbarT[1:])), t_of_wbarT_max=float(ts[1:][np.nanargmax(wbarT[1:])]),
               wbarT_end=float(wbarT[-1]), last_t_wbarT_above_third=last_above,
               last_above_over_R0=last_above / R0, max_absphi=float(max(maxphi)),
               wbarT_samples={str(float(t)): float(w) for t, w in zip(ts[1::max(1, len(ts) // 40)], wbarT[1::max(1, len(ts) // 40)])})
    return res


if __name__ == "__main__":
    jobs = [(30, 0.2, "cos2", 0.2), (100, 0.2, "cos2", 0.2), (400, 0.2, "cos2", 0.2), (100, 0.2, "cos2", 0.1),
            (100, 0.2063, "flat", 0.2), (400, 0.2063, "flat", 0.2), (100, 0.02, "cos2", 0.2)]
    out = []
    for R0, A, prof, dr in jobs:
        r = run(R0, A, prof, dr=dr, dt=min(0.05, dr / 4))
        out.append(r)
        print(json.dumps({k: v for k, v in r.items() if k != "wbarT_samples"}), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "results", "radial_lump.json"), "w") as fh:
        json.dump(out, fh, indent=1)
