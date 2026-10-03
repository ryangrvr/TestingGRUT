"""BRI1-PF4Q core (BRI1_PF4Q.md sec. 3): deterministic quadrature of Gibbs averages over the clamped single Duffing
oscillator.  NUMERICAL EVIDENCE ONLY -- not certified.  No Monte Carlo, no sampling of any kind.
Independent code path, not independent reviewer.

Method A: variational equations + tensor probabilists' Gauss-Hermite in (x0,p0) with factor exp(-x0^4/4) in the weight.
Method B: nonlinear clamped oscillator at small eps + uniform trapezoid in (x0,p0); K = kappa3(X^eps)/eps (odd in eps).
"""
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.integrate import solve_ivp

PI = np.pi
TAU = (PI, 1.5 * PI, 2 * PI)
T_V = PI / 8                                   # pre-declared validation time (code validation only)


def smooth(u):
    return 10 * u**3 - 15 * u**4 + 6 * u**5


def q_of(proto, t):
    if proto == 0:
        return 0.0
    if t <= PI:
        return smooth(t / PI)
    return 1.0 if proto == 1 else 1.0 + smooth((t - PI) / PI)


# ---------------------------------------------------------------- quadrature rules
def rule_gh(n, prune=1e-40):
    """Tensor probabilists' Gauss-Hermite with weight exp(-x^2/2 - x^4/4 - p^2/2); ratio-normalised; pruned."""
    z, w = hermegauss(n)
    X, P = np.meshgrid(z, z, indexing="ij"); W = np.outer(w * np.exp(-z**4 / 4), w)
    keep = W >= prune * W.max()
    return X[keep], P[keep], W[keep] / W[keep].sum(), float(W[~keep].sum() / W.sum())


def rule_trap(h, Lx=4.5, Lp=8.5):
    # exactly symmetric under (x,p) -> (-x,-p) (needed for exact oddness of kappa3 in eps on the quadrature measure)
    xs = h * np.arange(-int(Lx / h), int(Lx / h) + 1); ps = h * np.arange(-int(Lp / h), int(Lp / h) + 1)
    X, P = np.meshgrid(xs, ps, indexing="ij"); W = np.exp(-(X**2 / 2 + X**4 / 4 + P**2 / 2))
    return X.ravel(), P.ravel(), (W / W.sum()).ravel()


# ---------------------------------------------------------------- ODE runs (restart at t = pi)
def _run(rhs, y0, times, tol):
    out = {}; seg = [(0.0, PI), (PI, 2 * PI)]; y = y0
    for (a, b) in seg:
        te = sorted(t for t in times if a < t <= b + 1e-15)
        if b not in te: te.append(b)
        sol = solve_ivp(rhs, (a, b), y, method="DOP853", rtol=tol, atol=tol, t_eval=te)
        assert sol.success, sol.message
        for k, t in enumerate(sol.t):
            out[float(t)] = sol.y[:, k]
        y = sol.y[:, -1]
    return out


def variational(X0, P0, times, tol=1e-13):
    """State rows: x, p, y(P1), v(P1), y(P2), v(P2).  Returns {t: (x, y1_P1, y1_P2)}."""
    M = len(X0)
    def rhs(t, s):
        s = s.reshape(6, M); x, p, y1, v1, y2, v2 = s; k = 1 + 3 * x * x
        return np.concatenate([p, -x - x**3, v1, -k * y1 + q_of(1, t), v2, -k * y2 + q_of(2, t)])
    y0 = np.concatenate([X0, P0, np.zeros(4 * M)])
    out = _run(rhs, y0, times, tol)
    return {t: (v.reshape(6, M)[0], v.reshape(6, M)[2], v.reshape(6, M)[4]) for t, v in out.items()}


def nonlinear(X0, P0, eps, proto, times, tol=1e-12):
    M = len(X0)
    def rhs(t, s):
        x, p = s[:M], s[M:]
        return np.concatenate([p, -x - x**3 + eps * q_of(proto, t)])
    out = _run(rhs, np.concatenate([X0, P0]), times, tol)
    return {t: v[:M] for t, v in out.items()}


# ---------------------------------------------------------------- statistics
def k3(W, A, B, C):
    """Joint third cumulant of weighted samples (exact algebra on the quadrature measure)."""
    ma, mb, mc = W @ A, W @ B, W @ C
    return W @ ((A - ma) * (B - mb) * (C - mc))


def c_term(W, xa, xb, yc):
    return W @ (xa * xb * yc) - (W @ (xa * xb)) * (W @ yc)


IDX = [(a, b, c) for a in range(3) for b in range(a, 3) for c in range(b, 3)]   # 10 unique components


def K_from_variational(W, sol, proto):
    x = [sol[t][0] for t in TAU]; y = [sol[t][proto] for t in TAU]
    return {(a, b, c): c_term(W, x[a], x[b], y[c]) + c_term(W, x[a], x[c], y[b]) + c_term(W, x[b], x[c], y[a])
            for (a, b, c) in IDX}


def K_from_fd(W, X0, P0, proto, eps):
    sol = nonlinear(X0, P0, eps, proto, list(TAU)); x = [sol[t] for t in TAU]
    return {(a, b, c): k3(W, x[a], x[b], x[c]) / eps for (a, b, c) in IDX}
