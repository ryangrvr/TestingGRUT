"""Independent verification core for the thermal_potential module (key v_thermal).

Written from scratch; does not import or copy module code.  Units v = 1 (v = 1000 GeV).
Real fields (phi, q1, q2, D), X = (q1 + i q2)/sqrt(2).
V0 = C + mp phi^2/2 + mx R^2/2 + md D^2/2 + lp phi^4/4 + lx R^4/4 + ld D^4/4
     + c phi^2 R^2/2 + d phi^2 D^2/2 + al R^2 D^2/2,  R^2 = q1^2 + q2^2
with c = g_X^2, d = h_D^2.  (g_X^2 phi^2 |X|^2 = c phi^2 R^2/2, lambda_X |X|^4 = lx R^4/4.)
"""
import numpy as np
from scipy import integrate, special
from scipy.integrate import solve_ivp

PI = np.pi
L16 = 16 * PI**2

# ---------------------------------------------------------------- couplings
BENCH = dict(lp=0.01, lx=0.001, c=0.1**2, d=0.01**2, ld=0.001, al=0.0,
             mp=-0.01, mx=0.0, md=0.1**2)
MU0 = 0.1


def quartic_tensor(p):
    """lambda_ijkl = d^4 V / d f_i d f_j d f_k d f_l, basis (phi,q1,q2,D), built by
    finite symbolic rule from the monomials (exact for a quartic polynomial)."""
    lp, lx, c, d, ld, al = p['lp'], p['lx'], p['c'], p['d'], p['ld'], p['al']
    # V4(f) as a function; extract 4th derivatives by polarization of a homogeneous quartic
    def V4(f):
        ph, q1, q2, D = f
        R2 = q1 * q1 + q2 * q2
        return (lp * ph**4 / 4 + lx * R2**2 / 4 + ld * D**4 / 4 + c * ph**2 * R2 / 2
                + d * ph**2 * D**2 / 2 + al * R2 * D**2 / 2)
    n = 4
    L = np.zeros((n, n, n, n))
    E = np.eye(n)
    # polarization identity for symmetric 4-linear form: d^4V = sum over signs
    import itertools
    for i, j, k, l in itertools.product(range(n), repeat=4):
        s = 0.0
        for e1, e2, e3, e4 in itertools.product([1, -1], repeat=4):
            s += e1 * e2 * e3 * e4 * V4(e1 * E[i] + e2 * E[j] + e3 * E[k] + e4 * E[l])
        # sum_{signs} prod(e) V4(sum e_a x_a) = 2^4 * 4! * (1/4!) d^4V(x_i,x_j,x_k,x_l) * ... normalization:
        L[i, j, k, l] = s / 16.0
    return L


def mass_matrix(p):
    return np.diag([p['mp'], p['mx'], p['mx'], p['md']])


def extract(L, M, C=0.0):
    return dict(lp=L[0, 0, 0, 0] / 6, lx=L[1, 1, 1, 1] / 6, c=L[0, 0, 1, 1] / 2, d=L[0, 0, 3, 3] / 2,
                ld=L[3, 3, 3, 3] / 6, al=L[1, 1, 3, 3] / 2, mp=M[0, 0], mx=M[1, 1], md=M[3, 3])


def betas_generic(L, M):
    """Standard one-loop scalar betas: 16pi^2 b(l_ijkl) = l_ijmn l_mnkl + 2 perms,
    16pi^2 b(m2_ij) = l_ijkl m2_kl."""
    t1 = np.einsum('ijmn,mnkl->ijkl', L, L)
    bL = t1 + np.einsum('ikmn,mnjl->ijkl', L, L) + np.einsum('ilmn,mnjk->ijkl', L, L)
    bM = np.einsum('ijkl,kl->ij', L, M)
    return bL / L16, bM / L16


def run_couplings(p, mu, mu0=MU0):
    if abs(np.log(mu / mu0)) < 1e-15:
        return dict(p)
    L0 = quartic_tensor(p)
    M0 = mass_matrix(p)
    y0 = np.concatenate([L0.ravel(), M0.ravel()])

    def rhs(t, y):
        L = y[:256].reshape(4, 4, 4, 4)
        M = y[256:].reshape(4, 4)
        bL, bM = betas_generic(L, M)
        return np.concatenate([bL.ravel(), bM.ravel()])
    sol = solve_ivp(rhs, [0.0, np.log(mu / mu0)], y0, rtol=1e-12, atol=1e-16, method='DOP853')
    y = sol.y[:, -1]
    return extract(y[:256].reshape(4, 4, 4, 4), y[256:].reshape(4, 4))


# ---------------------------------------------------------------- Debye masses
def debye(p):
    """Pi_a/T^2 = (1/24) d_a d_a Tr(Hessian), from Tr H written out by hand."""
    lp, lx, c, d, ld, al = p['lp'], p['lx'], p['c'], p['d'], p['ld'], p['al']
    # Tr H = phi^2 (3lp + 2c + d) + R^2 (c + 4 lx + al) + D^2 (d + 3 ld + 2 al) + const
    return dict(phi=2 * (3 * lp + 2 * c + d) / 24, X=2 * (c + 4 * lx + al) / 24, D=2 * (d + 3 * ld + 2 * al) / 24)


def debye_numeric(p, h=1e-3):
    """Independent numeric check: Tr Hessian by finite differences, then its second derivative / 24."""
    def V(f):
        ph, q1, q2, D = f
        R2 = q1 * q1 + q2 * q2
        return (p['mp'] * ph**2 / 2 + p['mx'] * R2 / 2 + p['md'] * D**2 / 2 + p['lp'] * ph**4 / 4
                + p['lx'] * R2**2 / 4 + p['ld'] * D**4 / 4 + p['c'] * ph**2 * R2 / 2
                + p['d'] * ph**2 * D**2 / 2 + p['al'] * R2 * D**2 / 2)
    E = np.eye(4)

    def trH(f):
        f = np.asarray(f, float)
        return sum((V(f + h * E[i]) - 2 * V(f) + V(f - h * E[i])) / h**2 for i in range(4))
    f0 = np.array([0.3, 0.2, 0.1, 0.25])
    out = []
    H2 = 0.05
    for a in range(4):
        out.append((trH(f0 + H2 * E[a]) - 2 * trH(f0) + trH(f0 - H2 * E[a])) / H2**2 / 24)
    return out


# ---------------------------------------------------------------- thermal function J_B
ZETA = special.zeta
GE = np.euler_gamma
LN_AB = np.log(16 * PI**2) + 1.5 - 2 * GE
KMAX = 40
# my derivation: n != 0 Matsubara modes, order m^{2k}, k >= 3:
#   V_k = (T/2) (-1)^{k+1}/k m^{2k} Gamma(k-3/2)/((4pi)^{3/2} Gamma(k)) * 2 zeta(2k-3) / (2 pi T)^{2k-3}
#   J coefficient of (y^2)^k = V_k * 2 pi^2 / T^4 / (m^2/T^2)^k
_SER = []
for k in range(3, KMAX + 3):
    ck = (0.5 * (-1)**(k + 1) / k * special.gamma(k - 1.5) / ((4 * PI)**1.5 * special.gamma(k))
          * 2 * ZETA(2 * k - 3) / (2 * PI)**(2 * k - 3) * 2 * PI**2)
    _SER.append((k, ck))


def _jb_series(y2):
    """real part, valid |y| < 2 pi; y2 may be negative (analytic continuation, Re)."""
    y2 = np.asarray(y2, float)
    a = np.abs(y2)
    with np.errstate(divide='ignore', invalid='ignore'):
        lg = np.where(a > 0, np.log(np.where(a > 0, a, 1.0)) - LN_AB, 0.0)
    cube = np.where(y2 > 0, np.sqrt(np.clip(y2, 0, None))**3, 0.0)
    out = -PI**4 / 45 + PI**2 * y2 / 12 - PI / 6 * cube - y2**2 / 32 * lg
    for k, ck in _SER:
        out = out + ck * y2**k
    return out


def _jbp_series(y2):
    y2 = np.asarray(y2, float)
    a = np.abs(y2)
    with np.errstate(divide='ignore', invalid='ignore'):
        lg = np.where(a > 0, np.log(np.where(a > 0, a, 1.0)) - LN_AB, 0.0)
    sq = np.where(y2 > 0, np.sqrt(np.clip(y2, 0, None)), 0.0)
    out = PI**2 / 12 - PI / 4 * sq - y2 / 16 * lg - y2 / 32
    for k, ck in _SER:
        out = out + k * ck * y2**(k - 1)
    return out


def _jb_bessel(y2, kmax=80):
    y = np.sqrt(y2)
    s = 0.0
    for k in range(1, kmax + 1):
        s = s - y2 / k**2 * special.kn(2, k * y)
    return s


def _jbp_bessel(y2, kmax=80):
    y = np.sqrt(y2)
    s = 0.0
    for k in range(1, kmax + 1):
        s = s + y * special.k1(k * y) / (2 * k)
    return s


Y2SPLIT = 4.0


def Jb(y2):
    y2 = np.asarray(y2, float)
    out = np.empty_like(y2)
    m = y2 <= Y2SPLIT
    out[m] = _jb_series(y2[m])
    if np.any(~m):
        out[~m] = _jb_bessel(y2[~m])
    if np.any(y2 < -Y2SPLIT):
        raise ValueError('y2 < -4 not supported')
    return out


def Jbp(y2):
    y2 = np.asarray(y2, float)
    out = np.empty_like(y2)
    m = y2 <= Y2SPLIT
    out[m] = _jbp_series(y2[m])
    if np.any(~m):
        out[~m] = _jbp_bessel(y2[~m])
    if np.any(y2 < -Y2SPLIT):
        raise ValueError('y2 < -4 not supported')
    return out


def Jb_quad(y2):
    """direct quadrature of the defining integral; Re part for y2 < 0."""
    if y2 >= 0:
        f = lambda x: x * x * np.log1p(-np.exp(-np.sqrt(x * x + y2)))
        return integrate.quad(f, 0, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    a = np.sqrt(-y2)
    g = lambda x: x * x * np.log(2 * abs(np.sin(np.sqrt(max(a * a - x * x, 0.0)) / 2)) + 1e-300)
    f = lambda x: x * x * np.log1p(-np.exp(-np.sqrt(max(x * x - a * a, 0.0)) - 1e-300))
    # integrable log singularity at x = a on both sides
    I1 = integrate.quad(g, 0, a, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    I2 = integrate.quad(f, a, a + 1, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    I3 = integrate.quad(f, a + 1, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)[0]
    return I1 + I2 + I3


def Jbp_quad(y2):
    assert y2 >= 0
    f = lambda x: x * x / (2 * np.sqrt(x * x + y2)) / np.expm1(np.sqrt(x * x + y2)) if (x * x + y2) > 0 else 0.5 * x
    return integrate.quad(f, 0, np.inf, limit=400, epsabs=1e-15, epsrel=1e-13)[0]


# ---------------------------------------------------------------- single-dof pieces
def re_pow32(m2):
    m2 = np.asarray(m2, float)
    return np.where(m2 > 0, np.sqrt(np.clip(m2, 0, None))**3, 0.0)


def re_pow12(m2):
    m2 = np.asarray(m2, float)
    return np.where(m2 > 0, np.sqrt(np.clip(m2, 0, None)), 0.0)


def cw(m2, mu):
    m2 = np.asarray(m2, float)
    a = np.abs(m2)
    with np.errstate(divide='ignore', invalid='ignore'):
        lg = np.where(a > 0, np.log(np.where(a > 0, a, 1.0) / mu**2), 0.0)
    return m2**2 * (lg - 1.5) / (64 * PI**2)


def cwp(m2, mu):
    m2 = np.asarray(m2, float)
    a = np.abs(m2)
    with np.errstate(divide='ignore', invalid='ignore'):
        lg = np.where(a > 0, np.log(np.where(a > 0, a, 1.0) / mu**2), 0.0)
    return m2 * (lg - 1.0) / (32 * PI**2)


def vth(m2, T):
    return T**4 / (2 * PI**2) * Jb(np.asarray(m2, float) / T**2)


def vthp(m2, T):
    return T**2 / (2 * PI**2) * Jbp(np.asarray(m2, float) / T**2)


# ---------------------------------------------------------------- V_eff along phi (X = D = 0)
def masses_phi(phi, p):
    s = np.asarray(phi, float)**2
    return [(1, p['mp'] + 3 * p['lp'] * s, 3 * p['lp']),
            (2, p['mx'] + p['c'] * s, p['c']),
            (1, p['md'] + p['d'] * s, p['d'])]   # (n, m2, dm2/ds)


def V0_phi(phi, p):
    s = np.asarray(phi, float)**2
    return p['mp'] * s / 2 + p['lp'] * s**2 / 4


def Veff(phi, T, p, mu, scheme):
    Pi = debye(p)
    PIs = [Pi['phi'] * T**2, Pi['X'] * T**2, Pi['D'] * T**2]
    out = V0_phi(phi, p)
    for idx, (n, m2, _) in enumerate(masses_phi(phi, p)):
        P = PIs[idx]
        if scheme == 'oneloop':
            out = out + n * (cw(m2, mu) + vth(m2, T))
        elif scheme == 'AE':
            out = out + n * (cw(m2, mu) + vth(m2, T) - T / (12 * PI) * (re_pow32(m2 + P) - re_pow32(m2)))
        elif scheme == 'Parwani':
            out = out + n * (cw(m2 + P, mu) + vth(m2 + P, T))
        elif scheme == 'DR':
            if idx == 0:
                out = out + n * (cw(m2, mu) + vth(m2, T) + T / (12 * PI) * re_pow32(m2))
            else:
                out = out + n * (cw(m2, mu) + vth(m2, T) - T / (12 * PI) * (re_pow32(m2 + P) - re_pow32(m2)))
    return out


def dVds(phi, T, p, mu, scheme):
    """dV/d(phi^2); curvature at 0 is 2 dV/ds."""
    Pi = debye(p)
    PIs = [Pi['phi'] * T**2, Pi['X'] * T**2, Pi['D'] * T**2]
    s = np.asarray(phi, float)**2
    out = p['mp'] / 2 + p['lp'] * s / 2
    for idx, (n, m2, dm) in enumerate(masses_phi(phi, p)):
        P = PIs[idx]
        if scheme == 'oneloop':
            f = cwp(m2, mu) + vthp(m2, T)
        elif scheme == 'AE':
            f = cwp(m2, mu) + vthp(m2, T) - T / (8 * PI) * (re_pow12(m2 + P) - re_pow12(m2))
        elif scheme == 'Parwani':
            f = cwp(m2 + P, mu) + vthp(m2 + P, T)
        elif scheme == 'DR':
            if idx == 0:
                f = cwp(m2, mu) + vthp(m2, T) + T / (8 * PI) * re_pow12(m2)
            else:
                f = cwp(m2, mu) + vthp(m2, T) - T / (8 * PI) * (re_pow12(m2 + P) - re_pow12(m2))
        out = out + n * dm * f
    return out
