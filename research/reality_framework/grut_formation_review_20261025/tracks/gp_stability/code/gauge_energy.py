#!/usr/bin/env python3
"""Gauge-sector stability checks for the 7737487/b5174e5 point (track gp_stability, own code).

Uses the own fixed-charge profile (profile_fixed_charge.py).  Computes:
  * the Higgsed vector mass and screening length vs the object's size;
  * the screened (Yukawa) electrostatic energy with a REAL-SPACE radial Green's function,
    exact for a piecewise-linear source (Codex used a Fourier/Parseval route);
  * the local (large-m_A) limit  int n^2 /(72 h^2)  and its leading gradient correction;
  * an UNSCREENED Coulomb upper bound (two routes: Green's function and enclosed charge);
  * central potentials vs the binding per quantum, and thin-wall Q_max scalings;
  * heavy vacuum masses vs the object's total energy;
  * the residual Z6 algebra.
Conventions: D = d + i e q A, Gauss  div E = e (n_X + 3 n_Y + 6 n_Sigma), Heaviside-Lorentz.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid

from profile_fixed_charge import solve, fields, observables, V, MX, G

E_G = 1e-7          # new gauge coupling e
H = 1e11            # GeV, canonical Higgs amplitude, |Sigma| = h/sqrt2
M_Y = 1e12
MU = 5e11
LAM_S = 1.0
HBARC_M = 1.973269804e-16   # GeV m


def cell_weights(alpha, m, d):
    """Integrals over one cell of length d of exp(-m(d-t)) * [1, t/d] and exp(-m t) * [1, t/d]."""
    if alpha < 1e-4:
        a0 = d * (1 - alpha / 2 + alpha ** 2 / 6)
        a1 = d * (0.5 - alpha / 6 + alpha ** 2 / 24)    # int e^{-m(d-t)} t/d dt
        b1 = d * (0.5 - alpha / 3 + alpha ** 2 / 8)     # int e^{-m t} t/d dt
    else:
        em = -np.expm1(-alpha)            # 1 - e^{-alpha}
        a0 = em / m
        a1 = (d / m - em / m ** 2) / d
        b1 = (em / m ** 2 - d * np.exp(-alpha) / m) / d
    return a0, a1, b1


def yukawa_potential(r, s, m):
    """Solve (-lap + m^2) phi = s for a spherical source sampled on a uniform grid r (r[0]=0)."""
    d = r[1] - r[0]
    g = r * s
    if m == 0:
        enc = cumulative_trapezoid(r * r * s, r, initial=0.0)
        out = cumulative_trapezoid(g[::-1], -r[::-1], initial=0.0)[::-1]   # int_r^inf r' s
        rs = np.where(r > 0, r, 1.0)
        return np.where(r > 0, enc / rs, 0.0) + out
    alpha = m * d
    a0, a1, b1 = cell_weights(alpha, m, d)
    decay = np.exp(-alpha)
    n = len(r)
    I1 = np.zeros(n); I2 = np.zeros(n)
    for i in range(n - 1):        # forward: int_0^r e^{-m(r-r')} g
        I1[i + 1] = decay * I1[i] + g[i] * (a0 - a1) + g[i + 1] * a1
    b0 = a0                        # int_0^d e^{-m t} dt equals a0
    for i in range(n - 2, -1, -1):  # backward: int_r^inf e^{-m(r'-r)} g
        I2[i] = decay * I2[i + 1] + g[i] * (b0 - b1) + g[i + 1] * b1
    em = np.exp(-m * r)
    J1 = em * cumulative_trapezoid(em * g, r, initial=0.0)
    tot = np.trapezoid(em * g, r)
    J2 = em * (tot - cumulative_trapezoid(em * g, r, initial=0.0))
    rs = np.where(r > 0, r, 1.0)
    phi = (I1 - J1 + I2 - J2) / (2 * m * rs)
    # r = 0 limit: phi(0) = int_0^inf e^{-m r'} r' s dr'
    phi[0] = np.trapezoid(np.exp(-m * r) * g, r)
    return phi


def main():
    L = 70.0
    sol = solve(L, 3001)
    obs = observables(sol, L)
    Om = obs['Omega']
    r = np.linspace(0, L, 140001)
    F, Fp, S, Sp, q = fields(sol, r)
    s = F * F
    I4 = 4 * np.pi * simpson(r * r * s * s, x=r)
    pref = E_G ** 2 * V ** 4 * Om ** 2 / (2 * MX ** 3)       # GeV per unit of int F^2 phi d^3rho

    mA = 6 * E_G * H
    mbar = mA / MX
    # radii (in GeV^-1) of the X charge distribution
    n_shell = 4 * np.pi * r * r * s
    Qcum = cumulative_trapezoid(n_shell, r, initial=0.0); Qcum /= Qcum[-1]
    r50 = float(np.interp(0.5, Qcum, r)) / MX
    r90 = float(np.interp(0.9, Qcum, r)) / MX
    r99 = float(np.interp(0.99, Qcum, r)) / MX
    rrms = float(np.sqrt(simpson(r ** 4 * s, x=r) / simpson(r * r * s, x=r))) / MX
    half = float(r[np.argmax(F < 0.5 * F[0])]) / MX

    # local (large m_A) bound and its gradient correction
    N2 = Om ** 2 * V ** 4 * I4 / MX                            # int n^2 d^3x  [GeV^3]
    dE_local = N2 / (72 * H * H)
    grad_ratio = obs['int_gradF2_sq_d3rho'] / obs['int_F4_d3rho']   # in rho units
    rel_corr_leading = grad_ratio / mbar ** 2

    # real-space Yukawa ratios vs Codex's Fourier ratios [Q]
    codex_ratios = {0.1: 0.08165947240771287, 1.0: 0.8248639415303187, 10.0: 0.9976011667001629,
                    1000.0: 0.9999997590066865, mbar: 0.9999999999994883}
    ratios = []
    for m in [0.1, 1.0, 10.0, 1000.0, mbar]:
        phi = yukawa_potential(r, s, m)
        R = m * m * 4 * np.pi * simpson(r * r * s * phi, x=r) / I4
        ratios.append(dict(mA_over_mX=m, screened_over_local_realspace=float(R),
                           one_minus_ratio_series=float(grad_ratio / m ** 2),
                           codex_fourier_ratio_Q=codex_ratios[m],
                           screened_energy_GeV_at_this_mass=float(R * pref * I4 / m ** 2)))
    dE_screened = dE_local * (1 - rel_corr_leading)

    # unscreened Coulomb (m_A -> 0) upper bound, two routes
    phiC = yukawa_potential(r, s, 0.0)
    EC_green = pref * 4 * np.pi * simpson(r * r * s * phiC, x=r)
    Qphys = q / G ** 2                                          # enclosed X number (q from own BVP state)
    Qr = 200000.0 * Qcum
    rs = np.where(r > 0, r, 1.0)
    integrand = np.where(r > 0, (Qr / rs) ** 2, 0.0)
    EC_shell = E_G ** 2 * MX / (8 * np.pi) * (simpson(integrand, x=r) + 200000.0 ** 2 / L)
    R_uniform = np.sqrt(5.0 / 3.0) * rrms                       # uniform sphere with same rms radius
    EC_uniform = 0.6 * (E_G * 200000) ** 2 / (4 * np.pi * R_uniform)
    ePhi0_C = E_G ** 2 * MX * V ** 2 * Om * phiC[0]             # e*Phi(0) unscreened [GeV]
    n0 = Om * MX * V ** 2 * s[0]                                # GeV^3
    ePhi0_scr = n0 / (36 * H * H)                               # e*Phi(0) screened (local) [GeV]
    gap = (1 - Om) * MX
    margin = obs['binding_margin_GeV']
    b_per = margin / 200000
    # thin-wall unscreened scaling: E_C/Q ~ Q^(2/3) at fixed density
    Qmax_evap = 200000 * (b_per / (EC_green / 200000)) ** 1.5
    Qmax_phi = 200000 * (gap / ePhi0_C) ** 1.5

    # heavy masses
    mS = np.sqrt(2 * LAM_S) * H
    split = np.sqrt(2) * MU * H                     # mu * Sigma * (Y^2 + h.c.), Sigma = h/sqrt2
    mY = [np.sqrt(M_Y ** 2 - split), np.sqrt(M_Y ** 2 + split)]
    Eobj = obs['E_GeV']
    heavy = dict(vector_GeV=mA, Higgs_radial_GeV=mS, mediator_GeV=mY,
                 min_heavy_GeV=min(mA, mS, *mY), object_energy_GeV=Eobj,
                 min_heavy_over_object_energy=min(mA, mS, *mY) / Eobj,
                 objects_needed_to_exceed_mA=int(np.ceil(mA / Eobj)),
                 vector_mass_shift_from_core_X_GeV2=E_G ** 2 * (V * F[0]) ** 2,
                 induced_Y_amplitude_in_core_GeV=1e-3 * (V * F[0] / np.sqrt(2)) ** 3 / M_Y ** 2)

    # Z6 algebra: residual charge = (Q_X + 3 Q_Y) mod 6
    Q = 200000
    z6 = dict(Q_mod6=Q % 6, min_cores_neutral=next(k for k in range(1, 7) if (k * Q) % 6 == 0),
              extra_X_quanta_to_neutralize=(-Q) % 6, extra_antiX_quanta_alternative=Q % 6,
              AB_phase_over_2pi=(Q % 6) / 6, AB_sin2=float(np.sin(np.pi * (Q % 6) / 6) ** 2),
              note='Compact Gauss law: Q_X + 3 Q_Y + 6 Q_Sigma = 0, so (Q_X + 3 Q_Y) = 0 mod 6 for the whole slice.')

    out = dict(
        profile=obs,
        radii=dict(r_half_charge_GeVinv=r50, r90_GeVinv=r90, r99_GeVinv=r99, r_rms_GeVinv=rrms,
                   r_F_half_GeVinv=half, r_half_m=r50 * HBARC_M),
        screening=dict(mA_GeV=mA, mA_over_mX=mbar, screening_length_GeVinv=1 / mA,
                       screening_length_m=HBARC_M / mA, mA_times_r_half=mA * r50,
                       int_n2_d3x_GeV3=N2, dE_local_GeV=dE_local,
                       relative_gradient_correction=rel_corr_leading, dE_screened_GeV=dE_screened,
                       dE_over_margin=dE_screened / margin, realspace_vs_fourier=ratios,
                       local_quartic_equivalent=dict(
                           delta_lambda_eff=(Om * MX) ** 2 / (18 * H * H),
                           ratio_to_lambda_X=(Om * MX) ** 2 / (18 * H * H) / 1e-3,
                           note='int n^2/(72h^2) = int omega^2 f^4/(72 h^2) acts like a quartic omega^2/(18h^2)*f^4/4.')),
        unscreened=dict(EC_green_GeV=EC_green, EC_enclosed_charge_GeV=EC_shell,
                        EC_uniform_sphere_same_rms_GeV=EC_uniform, EC_over_margin=EC_green / margin,
                        EC_over_screened=EC_green / dE_screened,
                        ePhi0_unscreened_GeV=ePhi0_C, ePhi0_screened_GeV=ePhi0_scr,
                        gap_mX_minus_omega_GeV=gap, ePhi0_over_gap=ePhi0_C / gap,
                        binding_per_quantum_GeV=b_per,
                        Qmax_thin_wall_evaporation_estimate=Qmax_evap,
                        Qmax_thin_wall_potential_estimate=Qmax_phi),
        heavy=heavy, z6=z6)
    path = Path(__file__).resolve().parents[1] / 'results' / 'gauge_energy.json'
    path.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
