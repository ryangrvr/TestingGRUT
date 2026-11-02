"""qm_micro common module (Claude, 2 Nov 2026).  Own code.

M1 light sector in the rescaled units of tracks/m1_vacuum (m_X = g v = 1):
    L = (1/g^2) [ 1/2 |d psi|^2 + 1/2 (d S)^2 - U(S, |psi|) ],   X = psi m_X/(g sqrt2),  phi = v S,
    U = Uphi(S) + 1/2 S^2 |psi|^2 + a/4 |psi|^4,  Uphi = b/4 (S^2-1)^2 - c (S-1)^3,
    a = lambda_X/g^2 = 0.1, b = lambda_phi/g^2 = 3, c = eps/g^2 = 0.2, 1/g^2 = 100.
Charge density n = 100 Im(psi* d_t psi) with the sign chosen so that psi = F e^{-i Omega t} has n = 100 Omega F^2.
Contents:
  * bulk(): the zero-pressure bulk Q-matter state (S_B, F_B, Omega_b) solved to 1e-15;
  * planar_profile(): the 1D interface (bulk at x -> -inf, vacuum at x -> +inf) at Omega = Omega_b by sparse Newton
    on a uniform grid with a pinning condition F(0) = F_B/2 and Omega as the extra unknown (which must come out as
    Omega_b to exponential accuracy);
  * surface tension sigma = 100 int (F'^2 + S'^2) dx (planar, at coexistence).
"""
import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.optimize import fsolve

a_, b_, c_ = 0.1, 3.0, 0.2
G2 = 0.01
MPHI = math.sqrt(6.0)


def Uphi(S):
    return b_ / 4 * (S * S - 1) ** 2 - c_ * (S - 1) ** 3


def Uphi_p(S):
    return b_ * (S * S - 1) * S - 3 * c_ * (S - 1) ** 2


def Uphi_pp(S):
    return b_ * (3 * S * S - 1) - 6 * c_ * (S - 1)


def U(S, F):
    return Uphi(S) + 0.5 * S * S * F * F + a_ / 4 * F ** 4


def bulk():
    def eqs(x):
        S, F, Om = x
        return [S * S + a_ * F * F - Om * Om, Uphi_p(S) + F * F * S, 0.5 * Om * Om * F * F - U(S, F)]
    S, F, Om = fsolve(eqs, [0.1512, 2.4065, 0.7759], xtol=1e-13)
    # polish with Newton in extended precision-free form
    for _ in range(5):
        r = np.array(eqs([S, F, Om]))
        J = np.array([[2 * S, 2 * a_ * F, -2 * Om],
                      [Uphi_pp(S) + F * F, 2 * F * S, 0.0],
                      [-(Uphi_p(S) + S * F * F), Om * Om * F - S * S * F - a_ * F ** 3, Om * F * F]])
        d = np.linalg.solve(J, -r)
        S, F, Om = S + d[0], F + d[1], Om + d[2]
    return dict(S=S, F=F, Omega=Om, n=100 * Om * F * F, resid=float(np.max(np.abs(eqs([S, F, Om])))))


def planar_profile(Lb=60.0, Lv=60.0, h=0.01, x_pin=0.0, tol=1e-12, verbose=False):
    """Planar interface: bulk for x < 0 side (x = -Lb), vacuum at x = +Lv.  Returns dict with x, F, S, Omega."""
    B = bulk()
    SB, FB, OmB = B['S'], B['F'], B['Omega']
    N = int(round((Lb + Lv) / h))
    x = -Lb + h * np.arange(N + 1)
    xi = x[1:-1]
    n = xi.size
    ipin = int(np.argmin(np.abs(xi - x_pin)))
    prof = 0.5 * (1 - np.tanh(xi / 1.5))
    F = FB * prof
    S = 1 - (1 - SB) * prof
    Om = OmB
    k = h * h / 12.0
    for it in range(80):
        # Numerov residuals for F'' = G_F, S'' = G_S with Dirichlet values at both ends
        Fp = np.r_[FB, F, 0.0]; Sp = np.r_[SB, S, 1.0]
        GF = (Sp ** 2 + a_ * Fp ** 2 - Om ** 2) * Fp
        GS = Uphi_p(Sp) + Fp ** 2 * Sp
        RF = Fp[2:] - 2 * Fp[1:-1] + Fp[:-2] - k * (GF[2:] + 10 * GF[1:-1] + GF[:-2])
        RS = Sp[2:] - 2 * Sp[1:-1] + Sp[:-2] - k * (GS[2:] + 10 * GS[1:-1] + GS[:-2])
        Rpin = F[ipin] - 0.5 * FB
        R = np.r_[RF, RS, Rpin]
        nr = np.max(np.abs(R))
        if verbose:
            print(it, nr, Om)
        if nr < tol:
            break
        # Jacobian
        dGF_dF = Sp ** 2 + 3 * a_ * Fp ** 2 - Om ** 2
        dGF_dS = 2 * Sp * Fp
        dGS_dS = Uphi_pp(Sp) + Fp ** 2
        dGS_dF = 2 * Fp * Sp
        dGF_dOm = -2 * Om * Fp
        idx = np.arange(n)
        rows, cols, vals = [], [], []
        def add(r, c, v):
            rows.append(r); cols.append(c); vals.append(v)
        for (off_r, off_c, D, base) in ((0, 0, dGF_dF, True), (0, n, dGF_dS, False), (n, n, dGS_dS, True), (n, 0, dGS_dF, False)):
            d0 = D[1:-1]; dm = D[:-2]; dp = D[2:]
            add(idx + off_r, idx + off_c, (-2 if base else 0) - 10 * k * d0)
            add(idx[:-1] + off_r, idx[1:] + off_c, (1 if base else 0) - k * dp[:-1])
            add(idx[1:] + off_r, idx[:-1] + off_c, (1 if base else 0) - k * dm[1:])
        dRF_dOm = -k * (dGF_dOm[2:] + 10 * dGF_dOm[1:-1] + dGF_dOm[:-2])
        add(idx, np.full(n, 2 * n), dRF_dOm)
        add(np.array([2 * n]), np.array([ipin]), np.array([1.0]))
        J = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * n + 1, 2 * n + 1))
        d = spsolve(J, -R)
        lam = 1.0
        F, S, Om = F + lam * d[:n], S + lam * d[n:2 * n], Om + lam * d[2 * n]
    else:
        raise RuntimeError('profile Newton did not converge: %g' % nr)
    Fa = np.r_[FB, F, 0.0]; Sa = np.r_[SB, S, 1.0]
    # surface tension (excess grand potential per area) two ways: equipartition and direct
    dF = np.gradient(Fa, h, edge_order=2); dS = np.gradient(Sa, h, edge_order=2)
    sig_equip = 100 * np.trapezoid(dF ** 2 + dS ** 2, x)
    Veff = U(Sa, Fa) - 0.5 * Om * Om * Fa * Fa
    sig_direct = 100 * np.trapezoid(0.5 * dF ** 2 + 0.5 * dS ** 2 + Veff, x)
    # Gibbs dividing surface for charge: position where the charge excess vanishes
    nloc = 100 * Om * Fa * Fa
    xg = -Lb + np.trapezoid(nloc, x) / B['n']
    return dict(x=x, F=Fa, S=Sa, Omega=Om, Omega_b=OmB, S_B=SB, F_B=FB, n_B=B['n'], sigma_equip=sig_equip,
                sigma_direct=sig_direct, x_gibbs=xg, newton_its=it, resid=nr, h=h)


if __name__ == '__main__':
    import time
    t0 = time.time()
    B = bulk()
    print('bulk', B)
    P = planar_profile(verbose=True)
    print('Omega - Omega_b = %.3e  sigma = %.6f / %.6f  x_gibbs = %.4f  t=%.1fs' % (P['Omega'] - P['Omega_b'], P['sigma_equip'], P['sigma_direct'], P['x_gibbs'], time.time() - t0))
    P2 = planar_profile(h=0.02)
    print('h=0.02: sigma = %.6f  Om-Omb=%.2e' % (P2['sigma_equip'], P2['Omega'] - P2['Omega_b']))
