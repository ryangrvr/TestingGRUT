#!/usr/bin/env python3
"""Wall-core interaction energy at fixed X field (Claude review, track gp_walls).

A planar phi wall (kink v tanh(z/w), w = sqrt(2/lam)/v) and the Q=200000 core of core_profile.py,
with the core centre a distance d from the wall plane. Two estimates of
  U(d) = E[wall + core] - E[core] - sigma * Area   (X field held fixed => charge fixed)
1. product trial  S = S_Q(r) tanh(z/w): a valid configuration at the same charge, so the true
   minimum lies below it; U_trial(d) < 0 is therefore a lower bound on the attraction.
2. relaxed phi at d = 0: minimise the phi functional on an axisymmetric grid with X frozen
   (L-BFGS), boundary values -1/+1 at z -> -/+ infinity. Still an upper bound on U (X frozen).
Energies in GeV via E = (v/g) * int d^3rho e.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize

G, VPHI, LAMP = 0.1, 0.924089637375362, 0.03
B = LAMP / G**2
W = np.sqrt(2.0 / B)            # wall half-width in rho units (= m_X * sqrt(2/lam)/v)
EUNIT = VPHI / G


def load(npz):
    d = np.load(npz)
    x, F, S = d['x'], d['F'], d['S']
    return CubicSpline(x, F), CubicSpline(x, S), float(d['Om'])


def kink_density(z):
    t = np.tanh(z / W)
    return 0.5 * (1 - t * t)**2 / W**2 + 0.25 * B * (t * t - 1)**2


def sigma_dimless():
    z = np.linspace(-30, 30, 200001)
    return np.trapezoid(kink_density(z), z)


def U_trial(Fs, Ss, d, P=40.0, Z=45.0, h=0.02):
    rc = np.arange(h / 2, P, h)
    z = np.arange(-Z, Z + h / 2, h)
    RC, ZZ = np.meshgrid(rc, z, indexing='ij')
    r = np.sqrt(RC**2 + (ZZ - d)**2)
    rr = np.minimum(r, Fs.x[-1])
    SQ, dSQ, F = Ss(rr), Ss(rr, 1), Fs(rr)
    SQ = np.where(r > Fs.x[-1], 1.0, SQ); dSQ = np.where(r > Fs.x[-1], 0.0, dSQ); F = np.where(r > Fs.x[-1], 0.0, F)
    t = np.tanh(ZZ / W); sech2 = 1 - t * t
    St = SQ * t
    dr_c = t * dSQ * RC / r
    dz = t * dSQ * (ZZ - d) / r + SQ * sech2 / W
    e_trial = 0.5 * (dr_c**2 + dz**2) + 0.5 * St**2 * F**2 + 0.25 * B * (St**2 - 1)**2
    e_core = 0.5 * dSQ**2 + 0.5 * SQ**2 * F**2 + 0.25 * B * (SQ**2 - 1)**2
    diff = e_trial - e_core - kink_density(ZZ)
    return EUNIT * 2 * np.pi * np.sum(diff * RC) * h * h


def relaxed_d0(Fs, Ss, P=30.0, Z=30.0, h=0.25):
    """Minimise the phi energy with X frozen on an axisymmetric cell-centred grid.
    Outer boundary (radius P, z=+/-Z): wall-only profile tanh(z/W) for the wall runs, 1 for the core run.
    U = E[wall+core] - E[core] - E[wall], all with the same grid, so discretisation errors largely cancel."""
    rc = np.arange(h / 2, P, h)
    z = np.arange(-Z + h / 2, Z, h)
    RC, ZZ = np.meshgrid(rc, z, indexing='ij')
    r = np.sqrt(RC**2 + ZZ**2)
    F = np.where(r < Fs.x[-1], Fs(np.minimum(r, Fs.x[-1])), 0.0)
    SQ = np.where(r < Ss.x[-1], Ss(np.minimum(r, Ss.x[-1])), 1.0)
    nr, nz = RC.shape
    rf = rc + h / 2                                  # radial faces (last one is the outer boundary)
    wz = np.concatenate([RC, RC[:, -1:]], axis=1)    # weights for z-faces

    def make(bout, blo, bhi):
        def E_and_grad(flat, Fe):
            S = flat.reshape(nr, nz)
            gr = np.diff(np.concatenate([S, bout[None, :]], axis=0), axis=0) / h
            gz = np.diff(np.concatenate([np.full((nr, 1), blo), S, np.full((nr, 1), bhi)], axis=1), axis=1) / h
            pot = 0.5 * S**2 * Fe**2 + 0.25 * B * (S**2 - 1)**2
            E = h * h * (0.5 * np.sum(gr**2 * rf[:, None]) + 0.5 * np.sum(gz**2 * wz) + np.sum(pot * RC))
            g = np.zeros_like(S)
            t = gr * rf[:, None]
            g -= t; g[1:] += t[:-1]
            tz = gz * wz
            g -= tz[:, 1:]; g += tz[:, :-1]
            g *= h
            g += h * h * (S * Fe**2 + B * (S**2 - 1) * S) * RC
            return E, g.ravel()
        return E_and_grad

    wall_f = make(np.tanh(z / W), np.tanh((z[0] - h) / W), np.tanh((z[-1] + h) / W))
    core_f = make(np.ones(nz), 1.0, 1.0)
    opts = dict(maxiter=50000, ftol=1e-16, gtol=1e-11)
    runs = {}
    for label, fun, Fe, S0 in [('wall_only', wall_f, 0 * F, np.tanh(ZZ / W)),
                               ('wall_plus_core', wall_f, F, SQ * np.tanh(ZZ / W)),
                               ('core_only', core_f, F, SQ)]:
        m = minimize(fun, S0.ravel(), args=(Fe,), jac=True, method='L-BFGS-B', options=opts)
        runs[label] = m
    U = EUNIT * 2 * np.pi * (runs['wall_plus_core'].fun - runs['core_only'].fun - runs['wall_only'].fun)
    sigma_grid = runs['wall_only'].fun / (0.5 * P**2)
    Ufrozen_trial = EUNIT * 2 * np.pi * (wall_f(( SQ * np.tanh(ZZ / W)).ravel(), F)[0] - core_f(SQ.ravel(), F)[0]
                                         - wall_f(np.tanh(ZZ / W).ravel(), 0 * F)[0])
    return dict(grid_h=h, box_P=P, box_Z=Z, U_relaxed_d0_GeV=float(U), U_trial_same_grid_GeV=float(Ufrozen_trial),
                wall_only_sigma_dimless_on_grid=float(sigma_grid),
                converged={k: bool(v.success) for k, v in runs.items()},
                iterations={k: int(v.nit) for k, v in runs.items()})


def main(npz, outpath):
    Fs, Ss, Om = load(npz)
    sig = sigma_dimless()
    sig_GeV3 = sig * VPHI**2 * (G * VPHI)
    rows = []
    for d in [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12]:
        rows.append(dict(d_rho=d, d_GeVinv=d / (G * VPHI), U_trial_GeV=float(U_trial(Fs, Ss, d))))
    conv = dict(h002=rows[0]['U_trial_GeV'], h001=float(U_trial(Fs, Ss, 0.0, h=0.01)),
                bigbox=float(U_trial(Fs, Ss, 0.0, P=55.0, Z=60.0, h=0.02)))
    relaxed = [relaxed_d0(Fs, Ss, h=hh) for hh in (0.25, 0.15)]
    rhalf = float(Ss.x[np.argmin(np.abs(Ss(Ss.x) - 0.5))])
    geo = sig_GeV3 * np.pi * (rhalf / (G * VPHI))**2
    out = dict(sigma_dimless=sig, sigma_GeV3=sig_GeV3, sigma_analytic_GeV3=2 * np.sqrt(2 * LAMP) * VPHI**3 / 3,
               wall_width_rho=W, wall_width_GeVinv=W / (G * VPHI),
               U_trial_rows=rows, U_trial_convergence_d0=conv, relaxed_d0=relaxed,
               geometric_estimate_sigma_pi_rhalf2_GeV=float(geo),
               object_mass_GeV=17064.537454, binding_margin_GeV=1417.255,
               note='U<0 is attraction. X frozen: the true attraction at fixed charge is at least |U_relaxed| >= |U_trial|.')
    Path(outpath).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
