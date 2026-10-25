#!/usr/bin/env python3
"""Independent stationary-profile solver for the C1-HP-03/v1 two-field object.

Own formulation (Claude review, track `conversion`): fixed-frequency boundary
value problem in the regular variables u = rho*F and z = rho*(S-1), so there
is no 1/rho term.  Charge and energy are computed by quadrature of the
solution, the frequency for a given charge is found by root finding, and
dE/dQ is obtained by finite differences of (E(Omega), Q(Omega)).

Dimensionless equations (model definition quoted from bb211d1 MATHEMATICS.md
section 2, a = lambda_X/g^2 = 0.1, b = lambda_phi/g^2 = 3):
    F'' + 2F'/rho = (S^2 + a F^2 - Omega^2) F
    S'' + 2S'/rho = (b (S^2 - 1) + F^2) S
Charge  g^2 Q = 4 pi Omega int rho^2 F^2,   energy E = (v/g) * eb,
eb = 4 pi int rho^2 [ (F'^2+S'^2)/2 + Omega^2 F^2/2 + b(S^2-1)^2/4
                      + S^2 F^2/2 + a F^4/4 ].
These are the quoted model equations; the solver and quadrature are mine.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq

A, B = 0.1, 3.0
G = 0.1
V_PHI = 0.924089637375362      # GeV, quoted PARAMETER_POINT.json
M_X = G * V_PHI                # GeV
L = 60.0


def rhs(r, y, om):
    u, up, z, zp = y
    rs = np.where(r > 0, r, 1.0)
    F = np.where(r > 0, u / rs, up)
    S = np.where(r > 0, 1 + z / rs, 1 + zp)
    upp = (S * S + A * F * F - om * om) * u
    zpp = (B * (S * S - 1) + F * F) * S * r
    return np.vstack([up, upp, zp, zpp])


def bc(ya, yb, om):
    kF = np.sqrt(1 - om * om)
    # u ~ exp(-kF r) at large r; z decays faster than u, use z(L)=0.
    return np.array([ya[0], ya[2], yb[1] + kF * yb[0], yb[2]])


def guess(r, F0=2.62, S0=0.006, R=4.3, w=0.9):
    t = 1 / (1 + np.exp((r - R) / w))
    F = F0 * t
    S = 1 - (1 - S0) * t
    u = r * F
    z = r * (S - 1)
    return np.vstack([u, np.gradient(u, r), z, np.gradient(z, r)])


def solve(om, init=None, tol=1e-10, n=2001):
    r = np.linspace(0, L, n)
    y0 = guess(r) if init is None else init(r)
    sol = solve_bvp(lambda x, y: rhs(x, y, om), lambda a, b: bc(a, b, om), r, y0,
                    tol=tol, max_nodes=200000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def observables(sol, om, n=60001):
    r = np.linspace(0, L, n)
    u, up, z, zp = sol.sol(r)
    rs = np.where(r > 0, r, 1.0)
    F = np.where(r > 0, u / rs, up)
    S = np.where(r > 0, 1 + z / rs, 1 + zp)
    # derivatives of F and S from u, z: F' = (u' - F)/r, S' = (z' - (S-1))/r
    Fp = np.where(r > 0, (up - F) / rs, 0.0)
    Sp = np.where(r > 0, (zp - (S - 1)) / rs, 0.0)
    dens = 0.5 * (Fp ** 2 + Sp ** 2) + 0.5 * om ** 2 * F ** 2 + B / 4 * (S * S - 1) ** 2 \
        + 0.5 * S * S * F * F + A / 4 * F ** 4
    w = np.full(n, r[1] - r[0]); w[0] = w[-1] = w[0] / 2   # trapezoid; refined by Simpson below
    from scipy.integrate import simpson
    qn = 4 * np.pi * om * simpson(r * r * F * F, x=r)
    eb = 4 * np.pi * simpson(r * r * dens, x=r)
    grad = 4 * np.pi * simpson(r * r * 0.5 * (Fp ** 2 + Sp ** 2), x=r)
    kin = 4 * np.pi * simpson(r * r * 0.5 * om ** 2 * F ** 2, x=r)
    pot = 4 * np.pi * simpson(r * r * (B / 4 * (S * S - 1) ** 2 + 0.5 * S * S * F * F + A / 4 * F ** 4), x=r)
    virial = (grad + 3 * (pot - kin)) / eb     # Derrick: E_grad + 3(E_pot - E_kin) = 0
    return dict(Omega=om, qn=qn, Q=qn / G ** 2, eb=eb, E_GeV=V_PHI / G * eb,
                E_over_Qm=eb / qn, F0=float(F[0]), S0=float(S[0]), virial=virial,
                rms_residual=float(np.max(sol.rms_residuals)), nodes=int(sol.x.size))


class Continuation:
    """Keep the last solution as the initial guess for nearby frequencies."""

    def __init__(self):
        self.last = None

    def __call__(self, om):
        init = None if self.last is None else (lambda r, s=self.last: s.sol(r))
        sol = solve(om, init)
        self.last = sol
        return sol, observables(sol, om)


def omega_for_charge(cont, Q, lo, hi):
    target = G ** 2 * Q
    return brentq(lambda om: cont(om)[1]['qn'] - target, lo, hi, xtol=1e-14)


def packet_profile_check(csv_path, om):
    """Quadrature of the packet's published profile CSV (data check, not my solve)."""
    from scipy.integrate import simpson
    d = np.loadtxt(csv_path, delimiter=',', skiprows=1)
    r, F, Fp, S, Sp = d.T
    dens = 0.5 * (Fp ** 2 + Sp ** 2) + 0.5 * om ** 2 * F ** 2 + B / 4 * (S * S - 1) ** 2 \
        + 0.5 * S * S * F * F + A / 4 * F ** 4
    qn = 4 * np.pi * om * simpson(r * r * F * F, x=r)
    eb = 4 * np.pi * simpson(r * r * dens, x=r)
    return dict(Q=qn / G ** 2, E_over_Qm=eb / qn, E_GeV=V_PHI / G * eb, rho_max=float(r[-1]), points=len(r))


def main():
    out = {}
    cont = Continuation()
    # Locate the frequency for Q = 200000 starting from a coarse bracket.
    s, o = cont(0.86)
    print('seed', o, flush=True)
    om = omega_for_charge(cont, 200000., 0.855, 0.867)
    sol, base = cont(om)
    print('Q=200000', base, flush=True)
    out['Q200000'] = base
    # dE/dQ by symmetric finite differences in Omega (two step sizes).
    fd = []
    for h in [2e-4, 5e-5]:
        _, op = cont(om + h)
        _, om_ = cont(om - h)
        _, o0 = cont(om)
        fd.append(dict(dOmega=h, dEdQ_over_m=(op['eb'] - om_['eb']) / (op['qn'] - om_['qn']),
                       dQ_dOmega=(op['Q'] - om_['Q']) / (2 * h)))
    out['dE_dQ_finite_difference'] = fd
    # Complete-conversion charge loss 2 * 6260.9014 pairs (quoted pair count).
    pairs = 12521.802854958796 / 2
    Qlow = 200000. - 2 * pairs
    om_low = omega_for_charge(cont, Qlow, om, om + 0.03)
    _, low = cont(om_low)
    out['Q_after_complete_conversion'] = low
    dE = (base['eb'] - low['eb']) * V_PHI / G
    out['exact_core_energy_release_GeV'] = dE
    out['frozen_profile_release_GeV_2omega'] = 2 * om * M_X * pairs
    out['release_ratio_exact_over_frozen'] = dE / (2 * om * M_X * pairs)
    # Charge-weighted mean outgoing chemical energy 2*omega averaged over the depletion
    out['omega_over_m_after_complete_conversion'] = om_low
    print(json.dumps(out, indent=1), flush=True)
    if len(sys.argv) > 1:
        csv = Path(sys.argv[1])
        out['packet_profile_csv_quadrature_DATA_CHECK'] = packet_profile_check(csv, 0.8611951308217225)
        print(out['packet_profile_csv_quadrature_DATA_CHECK'])
    res = Path(__file__).resolve().parents[1] / 'results' / 'profile_omega.json'
    res.write_text(json.dumps(out, indent=2) + '\n')
    # Save my own Q=200000 profile for the scattering solver (scratch, not results).
    scratch = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/conversion')
    r = np.linspace(0, L, 120001)
    u, up, z, zp = sol.sol(r)
    rs = np.where(r > 0, r, 1.0)
    F = np.where(r > 0, u / rs, up)
    S = np.where(r > 0, 1 + z / rs, 1 + zp)
    np.savez(scratch / 'my_profile_Q200000.npz', r=r, F=F, S=S, Omega=om)


if __name__ == '__main__':
    main()
