"""VER0-B blind reproduction: DETERMINISTIC NUMERICAL SANITY CHECK ONLY (not part of any proof).

Single pre-chosen time t* = 0.5 (no scan over t), protocol P1.
Gibbs expectations by deterministic quadrature: Gauss-Hermite (probabilists') in p,
trapezoid on [-L, L] in x with weight exp(-x^2/2 - x^4/4) (spectrally accurate for this decay).
Checks:
 (a) m2 = E_mu[x^2] and Var_mu(x^2) = 1 - m2 - m2^2 > 0;
 (b) K(t*) = 3 Cov(X0(t*)^2, Y(t*)) from the variational ODE vs the t^7..t^11 series;
 (c) kappa3(X^eps(t*)) / eps for a few eps from the full nonlinear forced flow -> K(t*),
     and [kappa3/eps - K]/eps^2 roughly constant (odd parity => O(eps^3) remainder);
 (d) P0: kappa3(X(t*)) = 0 to quadrature precision.
"""
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.integrate import solve_ivp

TSTAR = 0.5
L, NX = 7.0, 281
xs = np.linspace(-L, L, NX)
wx = np.exp(-xs**2 / 2 - xs**4 / 4); wx[0] *= 0.5; wx[-1] *= 0.5
ps, wp = hermegauss(48)
X, P = np.meshgrid(xs, ps, indexing='ij')
W = np.outer(wx, wp); W /= W.sum()
x0, p0 = X.ravel(), P.ravel(); w = W.ravel(); n = x0.size

m2 = np.sum(w * x0**2); m4 = np.sum(w * x0**4)
print(f"(a) m2 = {m2:.12f}, m4 = {m4:.12f}, 1-m2 = {1-m2:.12f}, Var(x^2) = {m4-m2**2:.12f}")


def q(t):
    if t <= np.pi:
        u = t / np.pi
        return 10*u**3 - 15*u**4 + 6*u**5
    return 1.0


def rhs_var(t, s):
    x, v, y, yv = s[:n], s[n:2*n], s[2*n:3*n], s[3*n:]
    return np.concatenate([v, -x - x**3, yv, q(t) - (1 + 3*x**2)*y])


sol = solve_ivp(rhs_var, (0, TSTAR), np.concatenate([x0, p0, np.zeros(n), np.zeros(n)]),
                method='DOP853', rtol=1e-12, atol=1e-14)
XT, YT = sol.y[:n, -1], sol.y[2*n:3*n, -1]
K = 3 * (np.sum(w * XT**2 * YT) - np.sum(w * XT**2) * np.sum(w * YT))
pi = np.pi
c = {7: 3*(m2**2+m2-1)/(28*pi**3), 8: -9*(m2**2+m2-1)/(112*pi**4),
     9: (pi**2*m2**2+12*m2**2+12*m2+19*pi**2*m2-12-pi**2)/(672*pi**5),
     10: -(m2**2+19*m2-1)/(1120*pi**4),
     11: (2*m2**2+11*pi**2*m2**2+38*m2+38*pi**2*m2-60*pi**2-2)/(12320*pi**5)}
ser = sum(ck * TSTAR**k for k, ck in c.items())
print(f"(b) K(t*) quadrature = {K:.6e};  series through t^11 = {ser:.6e};  leading c7 t^7 = {c[7]*TSTAR**7:.6e}")
print(f"    E[X0(t*)^2] = {np.sum(w*XT**2):.12f} (should equal m2 by Gibbs invariance)")


def kappa3_forced(eps, prot):
    def rhs(t, s):
        x, v = s[:n], s[n:]
        f = eps * (q(t) if prot == 'P1' else 0.0)
        return np.concatenate([v, -x - x**3 + f])
    s = solve_ivp(rhs, (0, TSTAR), np.concatenate([x0, p0]), method='DOP853', rtol=1e-12, atol=1e-14)
    Xe = s.y[:n, -1]
    mu = np.sum(w * Xe); c2 = np.sum(w * (Xe - mu)**2); c3 = np.sum(w * (Xe - mu)**3)
    return c3, c2


for eps in [0.4, 0.2, 0.1]:
    k3, k2 = kappa3_forced(eps, 'P1')
    print(f"(c) eps={eps:4}: kappa3/eps = {k3/eps:.6e}, (kappa3/eps - K)/eps^2 = {(k3/eps - K)/eps**2:.4e}, kappa2 = {k2:.10f}")
k3, k2 = kappa3_forced(0.3, 'P0')
print(f"(d) P0: kappa3 = {k3:.3e}, kappa2 = {k2:.12f}")
