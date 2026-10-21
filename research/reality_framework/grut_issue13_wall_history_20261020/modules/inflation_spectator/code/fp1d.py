"""1D Starobinsky-Yokoyama Fokker-Planck for the coarse-grained phi (X = D = 0), units H = 1, time in e-folds.

    dP/dN = d/dy [ w'(y)/3 P ] + D d^2P/dy^2,   D = 1/(8 pi^2),   P_eq ~ exp(-s), s = (8 pi^2/3) w.

Discretisation: "square-root approximation" (SqRA) reversible Markov chain on a cell-centred grid,
k_{i->j} = (D/h^2) exp(-(s_j - s_i)/2), which satisfies detailed balance with pi_i ~ exp(-s_i) exactly.
For an even potential the phi-odd sector is the half chain y>0 with an extra killing rate 2D/h^2 on the
cell next to y=0. Then, for a start z on the + side,

    S(z,N) = E_z[sgn phi_N] = survival probability of the killed half chain,  p_-(z,N) = (1 - S)/2,

and the smallest odd eigenvalue Lambda_1 is the smallest eigenvalue of the killed half chain.
Large barriers: Lambda_1 = 1/E_qsd[tau] with the exact birth-death mean absorption time
T_m = sum_{j<=m} (sum_{l>=j} pi_l)/(pi_j r_j), evaluated in log space (no cancellation).
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import logsumexp
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from common import DIFF, BETA_FP


def make_grid(sp, h=None, smax=90.0):
    """Half grid y_m = (m+1/2) h, m = 0..M-1, extending to where s - s_min > smax."""
    ym = sp.ymin()
    s0 = sp.s(np.array([ym]))[0]
    # resolution: well width, barrier curvature, and an absolute cap
    scales = [0.02]
    if ym > 0:
        s2 = BETA_FP * 2 * sp.lp * ym ** 2
        scales.append(1.0 / np.sqrt(s2) / 25.0)
    if sp.a2 < 0:
        scales.append(1.0 / np.sqrt(BETA_FP * abs(sp.a2)) / 25.0)
    if h is None:
        h = min(scales)
    y = ym + 1.0
    while sp.s(np.array([y]))[0] - s0 < smax:
        y *= 1.15
    M = int(np.ceil(y / h))
    yy = (np.arange(M) + 0.5) * h
    return yy, h


def odd_chain(sp, yy, h):
    s = sp.s(yy)
    k = DIFF / h ** 2
    up = np.zeros_like(s)
    dn = np.zeros_like(s)
    up[:-1] = k * np.exp(-(s[1:] - s[:-1]) / 2)
    dn[1:] = k * np.exp(-(s[:-1] - s[1:]) / 2)
    kill = np.zeros_like(s)
    kill[0] = 2 * k
    return s, up, dn, kill, k


def odd_spectrum(sp, yy, h):
    s, up, dn, kill, k = odd_chain(sp, yy, h)
    diag = up + dn + kill
    off = -k * np.ones(len(s) - 1)
    lam, vec = eigh_tridiagonal(diag, off)
    return s, lam, vec


def lambda1_mfpt(sp, yy, h):
    """Lambda_1 = 1/E_pi[tau] for the killed half chain, all in log space."""
    s, up, dn, kill, k = odd_chain(sp, yy, h)
    logpi = -s
    # log(pi_j r_j): r_0 = kill[0], r_j = dn[j]
    logr = np.log(np.where(np.arange(len(s)) == 0, kill[0], np.maximum(dn, 1e-300)))
    # tail sums log sum_{l>=j} pi_l
    rev = np.logaddexp.accumulate(logpi[::-1])[::-1]
    term = rev - (logpi + logr)              # log of each summand in T_m
    logT = np.logaddexp.accumulate(term)     # log T_m
    logZ = logsumexp(logpi)
    logE = logsumexp(logpi + logT) - logZ
    return np.exp(-logE), logE


def p_minus_from_spectrum(s, lam, vec, m0, N):
    """p_-(N) for a start in half-grid cell m0 (exact spectral sum; accurate to ~1e-13 absolute)."""
    w = np.exp(-(s - s[m0]) / 2)              # sqrt(pi_m / pi_m0)
    A = vec.T @ w
    c = vec[m0, :] * A
    return 0.5 * np.sum(c * (-np.expm1(-lam * N)))


def equilibrium_stats(sp, yy, h):
    s = sp.s(yy)
    p = np.exp(-(s - s.min()))
    Z = p.sum() * h
    m2 = (p * yy ** 2).sum() * h / Z
    m4 = (p * yy ** 4).sum() * h / Z
    ipk = int(np.argmax(p))
    return {"phi2_over_H2": m2, "phi4_over_H4": m4, "peak_y": yy[ipk],
            "P_peak_over_P_origin": float(np.exp(-(s[ipk] - sp.s(np.array([0.0]))[0])))}


def p_minus_delta(sp, yy, h, m0, Ns):
    """p_-(N) for a delta start in half-grid cell m0, accurate in RELATIVE terms for any barrier.

    Reversibility maps the forward density to the backward semigroup: p_j(t)/pi_j = [e^{At} e_m0]_j / pi_m0,
    with A the killed half-chain generator. The absorbed probability is
        F(N) = kappa_0 (pi_0/pi_m0) int_0^N [e^{At} e_m0]_0 dt,
    where [e^{At} e_m0]_0 is the O(1e-2) probability of going from the barrier cell to cell m0, so no
    exponentially small number is ever formed by subtraction. p_- = F/2."""
    from scipy.sparse import bmat, csr_matrix
    s, up, dn, kill, k = odd_chain(sp, yy, h)
    M = len(s)
    diag = -(up + dn + kill)
    A = diags([dn[1:], diag, up[:-1]], [-1, 0, 1], format="csr")
    e0 = csr_matrix(([1.0], ([0], [0])), shape=(1, M))
    Aug = bmat([[A, None], [e0, csr_matrix((1, 1))]], format="csr")
    v0 = np.zeros(M + 1)
    v0[m0] = 1.0
    Ns = list(Ns)
    out = {}
    for N in Ns:
        vN = expm_multiply(Aug * N, v0)
        I = vN[-1]
        logF = np.log(kill[0]) - (s[0] - s[m0]) + np.log(max(I, 1e-300))
        out[N] = {"p_minus": 0.5 * np.exp(logF), "log10_p_minus": (logF + np.log(0.5)) / np.log(10.0),
                  "int_u0": I}
    return out


def p_minus_delta_curve(sp, yy, h, m0, Nmax=60.0, dn=0.25, atol=1e-24):
    """log p_-(n) on n = 0, dn, ..., Nmax for a delta start in cell m0 (same identity as p_minus_delta,
    integrated with a stiff BDF solver; validated against expm_multiply and the spectral sum)."""
    from scipy.sparse import bmat, csr_matrix
    from scipy.integrate import solve_ivp
    s, up, dn_, kill, k = odd_chain(sp, yy, h)
    M = len(s)
    A = diags([dn_[1:], -(up + dn_ + kill), up[:-1]], [-1, 0, 1], format="csr")
    e0 = csr_matrix(([1.0], ([0], [0])), shape=(1, M))
    Aug = bmat([[A, None], [e0, csr_matrix((1, 1))]], format="csr")
    v0 = np.zeros(M + 1)
    v0[m0] = 1.0
    ns = np.arange(0.0, Nmax + 1e-9, dn)
    sol = solve_ivp(lambda t, u: Aug @ u, (0.0, Nmax), v0, method="BDF", jac=Aug, rtol=1e-10, atol=atol,
                    t_eval=ns)
    I = sol.y[-1]
    floor = I < 100 * atol
    with np.errstate(divide="ignore"):
        logp = np.log(0.5 * kill[0]) - (s[0] - s[m0]) + np.log(np.where(I > 0, I, np.nan))
    return ns, logp, floor
