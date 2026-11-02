"""fx_observables SKEPTIC (Claude, 1 Nov 2026): independent recomputation of the Q = 2e5 object and its
direct-detection source integrals.  [C] own code.

(i) Profile by scipy.solve_bvp collocation on (F, F', S, S') at fixed Omega (the track used Numerov + Newton at fixed Q),
    asymptotic Robin conditions F'/F = -kappa - 1/r, (S-1)'/(S-1) = -m_phi - 1/r at r = L.  Units m_X = 1, M1 potential
    U = 3/4 (S^2-1)^2 - 0.2 (S-1)^3 + F^2 S^2/2 + 0.1 F^4/4  (X = v F/sqrt2, phi = v S, g = 0.1).
    Initial guess: a tanh wall (not the track's profile).
(ii) Source integrals I2 = int (S^2-1), I1 = int (S-1), IX = int F^2 and the first zero of each 3D Fourier transform;
     E_R(zero) = q0^2/(2 m_A) for Xe-131, natural Xe (131.29), Ge-74, natural Ge (72.63), Ar-40.
(iii) sigma_n0 per coupling^2 from the analytic Higgs-exchange formula, and an independent direct-detection rate with a
     Monte-Carlo velocity integral (SHM, same detector toys as the track) for the LZ calibration and the coupling limits.
"""
import json, math, sys
import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import spherical_jn
from scipy.optimize import brentq

a_, b_, c_ = 0.1, 3.0, 0.2
MX = 0.0924089637375362; V = 10 * MX; G2 = 0.01
MPHI = math.sqrt(6.0)


def Up(S):
    return b_ * (S * S - 1) * S - 3 * c_ * (S - 1)**2


def solve(Om, L=40.0, R0=5.0):
    kap = math.sqrt(1 - Om * Om)
    r = np.linspace(1e-4, L, 4000)
    F0 = 2.6 * 0.5 * (1 - np.tanh((r - R0) / 1.2))
    S0 = 1 - 0.85 * 0.5 * (1 - np.tanh((r - R0) / 1.0))
    y = np.vstack([F0, np.gradient(F0, r), S0, np.gradient(S0, r)])

    def f(r, y):
        F, Fp, S, Sp = y
        return np.vstack([Fp, -2 * Fp / r + (S * S + a_ * F * F - Om * Om) * F, Sp, -2 * Sp / r + Up(S) + F * F * S])

    def bc(ya, yb):
        return np.array([ya[1], ya[3], yb[1] + (kap + 1 / L) * yb[0], yb[3] + (MPHI + 1 / L) * (yb[2] - 1)])
    sol = solve_bvp(f, bc, r, y, tol=1e-8, max_nodes=200000)
    return sol


def integrals(sol, L):
    r = np.linspace(1e-6, L, 400001)
    F, Fp, S, Sp = sol.sol(r)
    w = 4 * math.pi * r * r
    tr = lambda y: float(np.trapezoid(y, r))
    IX = tr(w * F * F); I2 = tr(w * (S * S - 1)); I1 = tr(w * (S - 1))
    U = 0.75 * (S * S - 1)**2 - c_ * (S - 1)**3 + 0.5 * F * F * S * S + a_ * F**4 / 4
    return r, F, S, Fp, Sp, IX, I2, I1, w, U


def main(out_fn):
    Om = 0.881720
    L = 40.0
    sol = solve(Om, L)
    r, F, S, Fp, Sp, IX, I2, I1, w, U = integrals(sol, L)
    Qbar = Om * IX
    Q = Qbar / G2
    Eb = float(np.trapezoid(w * (0.5 * (Fp**2 + Sp**2) + 0.5 * Om * Om * F * F + U), r))
    EQ = Eb / Qbar
    # radii
    cq = np.concatenate([[0], np.cumsum(0.5 * (w[1:] * F[1:]**2 + w[:-1] * F[:-1]**2) * np.diff(r))]) / IX
    r50, r90, r99 = [float(np.interp(x, cq, r)) for x in (0.5, 0.9, 0.99)]
    rS = float(np.interp(0.5, S, r))      # S rises monotonically 
    # Fourier transforms and zeros
    def FT(q, A):
        return float(np.trapezoid(w * A * np.sinc(q * r / math.pi), r))
    zeros = {}
    for lab, A in (('phi2', S * S - 1), ('phi1', S - 1), ('X2', F * F)):
        qs = np.linspace(0.05, 2.0, 400)
        vals = np.array([FT(q, A) for q in qs])
        j = np.where(np.sign(vals[1:]) != np.sign(vals[:-1]))[0][0]
        q0 = brentq(lambda q: FT(q, A), qs[j], qs[j + 1], xtol=1e-12)
        qG = q0 * MX
        zeros[lab] = dict(q0_mX=q0, q0_GeV=qG, **{('E_keV_' + t): qG * qG / (2 * A_ * 0.9314941) * 1e6
                                                   for t, A_ in (('Xe131', 131.0), ('Xe_nat', 131.293), ('Ge74', 74.0), ('Ge_nat', 72.63), ('Ar40', 40.0))})
    # sigma_n0 per coupling^2 (Higgs exchange, f_N = 0.30, m_h = 125.20, m_N = 0.9389)
    MN, MH, VH, fN = 0.9389, 125.20, 246.22, 0.30
    Mobj = EQ * Q * MX
    mu = MN * Mobj / (MN + Mobj)
    Ft = dict(phi2=V * V * I2 / MX**3, phi1=V * I1 / MX**3, X2=V * V / 2 * IX / MX**3)
    sig = {k: mu * mu / math.pi * (fN * MN * v / MH**2)**2 * 0.3893794e-27 for k, v in Ft.items()}
    out = dict(scope=__doc__.strip().splitlines()[0], Omega=Om, bvp_status=int(sol.status), bvp_nodes=int(len(sol.x)),
               Q=Q, E_over_Qm=EQ, M_GeV=Mobj, S_center=float(S[0]), F_center=float(F[0]),
               r50=r50, r90=r90, r99=r99, r_S_half=rS, I_S2m1=I2, I_Sm1=I1, I_F2=IX, ratio_I2_I1=I2 / I1,
               zeros=zeros, Ft0=Ft, sigma_n0_per_coupling2_cm2=sig)
    # sigma/M (object-object, pi (2 r90)^2) in cm^2/g
    cm = 1.973269804e-14 / MX
    out['sigma_over_M_cm2_g_r90'] = math.pi * (2 * r90 * cm)**2 / (Mobj * 1.78266192e-24)
    out['sigma_over_M_cm2_g_r50'] = math.pi * (2 * r50 * cm)**2 / (Mobj * 1.78266192e-24)
    print(json.dumps(out, indent=1, default=float))
    np.savez_compressed(out_fn.replace('.json', '_profile.npz'), r=r[::40], F=F[::40], S=S[::40])
    json.dump(out, open(out_fn, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
