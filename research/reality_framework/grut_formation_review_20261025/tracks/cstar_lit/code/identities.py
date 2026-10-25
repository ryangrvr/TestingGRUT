#!/usr/bin/env python3
"""Independent checks of the elementary identities that 8986d34's C* documents rely on, plus a
Monte Carlo reproduction of the Harlow-Usatyuk-Zhao (2501.02359) code formulas that define one
published 'observer dictionary'. Own code (track cstar_lit); no Codex code is imported or run.

Writes results/identities.json.   Usage:  python3 -I code/identities.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(20261025)

# ---------------------------------------------------------------- constants (CODATA 2018 / IAU)
c = 299792458.0
hbar = 1.054571817e-34
G = 6.67430e-11
Mpc = 3.0856775814913673e22
GeV = 1.602176634e-10          # J
hbar_GeV_s = 6.582119569e-25


def de_sitter(H0_km_s_Mpc=67.4, omega_L=0.685):
    """S_dS two ways: (a) area law pi c^5/(hbar G H^2); (b) reduced-Planck 8 pi^2 M^2/H^2.
    Also Bousso's N = 3 pi/(Lambda l_P^2) with Lambda = 3 H^2/c^2, and rho_L = 3 M^2 H^2."""
    H = H0_km_s_Mpc * 1e3 / Mpc * np.sqrt(omega_L)                 # 1/s
    S_area = np.pi * c**5 / (hbar * G * H**2)
    Mred_GeV = np.sqrt(hbar * c**5 / (8 * np.pi * G)) / GeV         # reduced Planck mass in GeV
    H_GeV = hbar_GeV_s * H
    S_red = 8 * np.pi**2 * Mred_GeV**2 / H_GeV**2
    lP2 = hbar * G / c**3
    Lam = 3 * H**2 / c**2
    N_bousso = 3 * np.pi / (Lam * lP2)
    rho = 3 * Mred_GeV**2 * H_GeV**2
    return dict(H_Lambda_per_s=H, M_reduced_GeV=Mred_GeV, S_area_nats=S_area, S_reducedPlanck_nats=S_red,
                N_bousso=N_bousso, rho_Lambda_GeV4=rho,
                identity_rho_eq_24pi2M4_over_S_relerr=abs(rho * S_red / (24 * np.pi**2 * Mred_GeV**4) - 1),
                note='All three entropy expressions are the same quantity; Bousso N-bound uses Lambda>0 and is trivial as Lambda->0.')


def closed_de_sitter_constraint():
    """Friedmann (Hamiltonian) constraint for k=+1: (adot/a)^2 + 1/a^2 = rho/(3M^2) = H^2.
    a(t)=cosh(Ht)/H, check by finite differences on a grid (H=1 units). Also acceleration equation
    addot/a = -(rho+3p)/(6M^2) with p=-rho gives +H^2."""
    t = np.linspace(-3, 3, 6001)
    a = np.cosh(t)
    ad = np.gradient(a, t, edge_order=2)
    add = np.gradient(ad, t, edge_order=2)
    con = (ad / a) ** 2 + 1 / a**2 - 1.0
    acc = add / a - 1.0
    sl = slice(5, -5)
    return dict(max_constraint_residual=float(np.max(np.abs(con[sl]))),
                max_acceleration_residual=float(np.max(np.abs(acc[sl]))),
                pressure_over_density=-1.0,
                meaning='Compact S^3 slices, positive vacuum density, p=-rho: all GR constraints hold for any H>0. '
                        'Closure fixes neither rho_Lambda nor p=0.')


def power_law_fading():
    """a ~ t^2: H=2/t, w from a ~ t^(2/(3(1+w))), rho ~ H^2 ~ t^-2, addot>0, event horizon
    d_E = a(t) * int_t^inf dt'/a(t') = t (numerical quadrature at t=1..10)."""
    w = 2 / (3 * 2.0) - 1
    out = []
    for t0 in [1.0, 3.0, 10.0]:
        # int_t0^T t'^-2 dt' on a log grid up to T = 1e8 t0, plus analytic tail 1/T
        s = np.linspace(np.log(t0), np.log(1e8 * t0), 400001)
        tp = np.exp(s)
        f = tp ** -2 * tp                      # dt = t ds
        integral = float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(s))) + 1 / tp[-1]
        out.append(dict(t=t0, event_horizon_over_t=float(t0**2 * integral / t0)))
    return dict(w=w, accelerating=True, rho_power=-2, horizons=out,
                meaning='Fading density with w>-1 can keep an event horizon; no Minkowski future follows.')


def gauss_law_lattice(n=12, trials=5):
    """Compact T^3 lattice: for arbitrary link fields E_i(x), sum_x div E = 0 identically
    (periodic differences telescope). Hence Q_total = sum rho = 0 for any solution of div E = rho.
    Second check: Poisson eq lap(phi) = -rho on T^3 has no solution when sum rho != 0 (zero mode)."""
    res = []
    for _ in range(trials):
        E = rng.normal(size=(3, n, n, n))
        div = sum(E[i] - np.roll(E[i], 1, axis=i) for i in range(3))
        res.append(float(abs(div.sum()) / np.abs(div).sum()))
    rho = np.zeros((n, n, n)); rho[0, 0, 0] = 1.0           # single net charge
    rho_k = np.fft.fftn(rho)
    return dict(max_relative_total_divergence=max(res), single_charge_zero_mode=float(abs(rho_k[0, 0, 0])),
                meaning='Total gauge charge on a compact boundaryless space vanishes; a lone net charge has no '
                        'Gauss-law solution (the k=0 Fourier mode of rho must vanish). Global (ungauged) charges '
                        'are not constrained by this.')


def finite_clock_toy():
    """Page-Wootters-type toy: H = Z(x)1 - 1(x)Z. Kernel dimension and relational conditioning."""
    Z = np.diag([1.0, -1.0]); I = np.eye(2)
    Hc = np.kron(Z, I) - np.kron(I, Z)
    ker = int(np.sum(np.abs(np.linalg.eigvalsh(Hc)) < 1e-12))
    psi = np.array([1, 0, 0, 1], complex) / np.sqrt(2)
    sx = np.array([[0, 1], [1, 0]])
    vals = []
    for th in np.linspace(0, np.pi, 5):
        clock_bra = np.array([1, np.exp(-1j * th)]) / np.sqrt(2)
        s = clock_bra @ psi.reshape(2, 2); s = s / np.linalg.norm(s)
        vals.append([float(th), float(np.real(np.conj(s) @ sx @ s)), float(np.cos(th))])
    return dict(constraint_kernel_dimension=ker, conditional_sigma_x_vs_cos=vals,
                meaning='A constraint alone does not force a one-dimensional physical space; the one-state claim '
                        'comes from the nonperturbative path integral / factorization, not from H=0.')


def gram_rank_test(k=6, trials=200):
    """UWZ/HUZ criterion: PSD Gram matrix M has rank one iff Tr(M)^n = Tr(M^n) for all n.
    Check on random rank-1 and rank-2 PSD matrices."""
    out = {}
    for r in (1, 2):
        dev = []
        for _ in range(trials):
            A = rng.normal(size=(k, r)) + 1j * rng.normal(size=(k, r))
            M = A @ A.conj().T
            dev.append(max(abs(np.trace(M) ** n - np.trace(np.linalg.matrix_power(M, n))) / abs(np.trace(M)) ** n
                           for n in (2, 3, 4)))
        out[f'rank{r}_max_relative_violation'] = float(np.max(dev)) if r == 1 else float(np.min(dev))
    return out


def haar_unit_vector(d):
    v = rng.normal(size=d)
    return v / np.linalg.norm(v)


def huz_code(d_M=4, d_f=1, d_obs_list=(1, 2, 4, 8, 16), samples=40000):
    """Monte Carlo of the HUZ code. Without observer: V|psi> = sqrt(d) <0|O(|psi>|psi0>) (a scalar),
    d = d_M d_f. Only u = O^T|0> (uniform on S^{d-1}) enters. Fluctuation of <phi|V^dag V|psi> - <phi|psi>
    should be ~ 1 + |<phi*|psi>|^2 (HUZ eq. 1.3). With cloned observer of d_obs equally weighted pointer
    states, |omega> = sum_a |a>_Ob' |a>_Ob / sqrt(d_obs), V-hat: H_M -> H_Ob' and the fluctuation should
    be ~ e^{-S_2} + Tr(w w^T)|<phi*|psi>|^2 = (1 + |<phi*|psi>|^2)/d_obs (HUZ eq. 4.2), up to O(1/d).
    Also: rank of V-hat on any set of k states is <= d_obs exactly (image lies in H_Ob')."""
    phi = np.zeros(d_M, complex); phi[0] = 1
    psi = np.zeros(d_M, complex); psi[0] = 1 / np.sqrt(2); psi[1] = 1j / np.sqrt(2)   # <phi*|psi> = 1/sqrt2
    ip = np.vdot(phi, psi); ipstar = np.vdot(phi.conj(), psi)
    rows = []
    for d_obs in d_obs_list:
        d = d_obs * d_M * d_f
        acc = 0.0
        for _ in range(samples):
            u = haar_unit_vector(d).reshape(d_obs, d_M, d_f)[:, :, 0] if d_f == 1 else None
            # amplitude for observer pointer a: sqrt(d) * u[a] . psi (psi0 = first f basis state)
            amp_psi = np.sqrt(d) * (u @ psi) / np.sqrt(d_obs)
            amp_phi = np.sqrt(d) * (u @ phi) / np.sqrt(d_obs)
            val = np.vdot(amp_phi, amp_psi)           # sum_a conj(phi_a) psi_a = <phi|Vhat^dag Vhat|psi>
            acc += abs(val - ip) ** 2
        mc = acc / samples
        pred_leading = (1 + abs(ipstar) ** 2) / d_obs
        # exact finite-d Haar value: E[d^2 u_i u_j u_k u_l] = d/(d+2) (sym. deltas); evaluated directly
        exact = (d / (d + 2)) * (abs(ip) ** 2 + (1 / d_obs) * (1 + abs(ipstar) ** 2)) - abs(ip) ** 2
        rows.append(dict(d_obs=d_obs, S2_obs=float(np.log(d_obs)), d=d, monte_carlo=float(mc),
                         huz_leading=float(pred_leading), exact_finite_d=float(exact)))
    # rank bound: map k = 3*d_obs generic states through one fixed V-hat
    rank_rows = []
    for d_obs in (2, 4):
        dM = 3 * d_obs
        d = d_obs * dM
        u = haar_unit_vector(d).reshape(d_obs, dM)
        X = rng.normal(size=(dM, dM)) + 1j * rng.normal(size=(dM, dM))   # dM generic states
        images = np.sqrt(d) * (u @ X) / np.sqrt(d_obs)                    # columns: V-hat|x_j> in C^{d_obs}
        rank_rows.append(dict(d_obs=d_obs, n_input_states=dM,
                              rank_of_images=int(np.linalg.matrix_rank(images, tol=1e-9))))
    return dict(fluctuations=rows, rank_bound=rank_rows,
                meaning='Effective inner products are reliable only to ~e^{-S_Ob}, and the effective space has '
                        'dimension <= e^{S_Ob} (dimension of the clone space). An e^{10^122}-state effective '
                        'description in this dictionary needs an observer with S_Ob >~ 10^122 nats, i.e. of order '
                        'the de Sitter entropy itself.')


def main():
    out = dict(
        scope='Elementary identities and a reproduction of published code formulas; no cosmological inference.',
        de_sitter=de_sitter(),
        closed_de_sitter=closed_de_sitter_constraint(),
        power_law_fading=power_law_fading(),
        gauss_law=gauss_law_lattice(),
        finite_clock=finite_clock_toy(),
        gram_rank=gram_rank_test(),
        huz_code=huz_code(),
        one_dimensional_algebra=dict(dim_End_C=1, max_orthogonal_nonzero_projectors=1),
    )
    (ROOT / 'results' / 'identities.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
