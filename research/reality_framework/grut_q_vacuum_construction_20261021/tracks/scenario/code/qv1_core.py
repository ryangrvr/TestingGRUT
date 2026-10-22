"""QV1 core: the closed flat-FRW system of q-theory (Klinkhamer-Volovik action), in dimensionless form.

Track `scenario` (stage 2), 21 October 2026. NOTES.md section 1 states every equation used here.

Action (BRIEF.md; KV 0806.2805 eqs (2.1)-(2.8), credited):
    S = int d^4x sqrt(-g) [ K(q) R / 2 - eps(q) ] + S_m[psi, g],   q = F~/sqrt(-g)  (3-form A, F = dA)
Field equations (equations track sections 1.2-1.4; literature track section 2):
    (q)  eps'(q) - K'(q) R / 2 = mu            mu = global integration constant (the state)
    (g)  K G_mn = T_mn - (eps - mu q) g_mn - (g_mn Box - nabla_m nabla_n) K
    (m)  nabla^m T_mn = 0
State condition (P_q): mu = mu0[eps_inf], eps_inf = eps(q) + every q-independent vacuum energy of the final
T -> 0 vacuum. Under P_q the matter equation of state is normalised so that p(T=0) = 0 (its vacuum energy is part
of eps_inf), and the gravitating vacuum term is rho_V = eps_inf(q) - mu0 q, which vanishes at q0.

Dimensionless variables (units K0 = K(q0) = 1; the caller chooses the time unit):
    Y = (K(q) - K0)/K0 = K'(q0)(q - q0)/K0      (K linear in q, as in KV 0806.2805: K = q/3 with sigma = 1)
    V = dY/dt,  H = adot/a,  N = ln a,  J = int H^2 dt = int H dN (for N-space window averages)
eps_inf is quartic about q0 (exact for the KV form (4.7)):
    rho_V(Y) = (3/4) m^2 Y^2 (1 + b3 Y + b4 Y^2) + c (1 + beta Y)          [c = -(mu - mu0) q0: offset variants]
    R(Y, V)  = 3 m^2 Y (1 + 1.5 b3 Y + 2 b4 Y^2) + 2 beta c + 3 Gamma V/(1+Y)^2
    KV values (eps = (1/2chi)(-q^2/q0^2 + q^4/(3 q0^4)), K = q/3): b3 = 4/3, b4 = 1/3, beta = K0/(K'q0) = 1,
    m^2 = (2/3) K0 eps''(q0)/K'^2 = 2/(3 chi) in units M = 1.
First-order system (trace of (g), (q), FRW):
    dY/dt = V
    dV/dt = [S + 4 rho_V - (1+Y) R]/3 - 3 H V
          = S/3 + (c/3)[(4 - 2 beta) + 2 beta Y] - m^2 F(Y) - 3 H V - Gamma V/(1+Y)
            with F(Y) = Y + 1.5 b3 Y^2 + (b3/2 + 2 b4) Y^3 + b4 Y^4
    dH/dt = R/6 - 2 H^2
    matter: model dependent; S = rho_m - 3 p_m (minus the trace of T^m)
Constraint (Friedmann, 00):  C = 3 H^2 (1+Y) + 3 H V - rho_m - rho_V = 0;  dC/dt = -4 H C (Bianchi; verified in
s2 with sympy, including the dissipative term).
Dissipative extension (only in section 3, labelled QV1+Gamma, an effective source [Q/E]): the q-equation becomes
eps' - K'R/2 = mu - Gamma_q qdot with Gamma_q = (3/2) Gamma K0 K'^2/K^2, the metric equation keeps the constant mu,
and radiation receives Q = Gamma_q qdot^2 = (3/2) Gamma K0 V^2/(1+Y)^2. Gamma = N_s m^3/(192 pi M^2) is the scalaron
decay rate into N_s minimally coupled scalars (Starobinsky 1980; Vilenkin 1985; Gorbunov-Panin 1009.2448 eq (T4+)).
"""
import math

import numpy as np
from scipy.integrate import ode, solve_ivp

KV_SHAPE = dict(b3=4.0 / 3.0, b4=1.0 / 3.0, beta=1.0)


class QV1:
    """q-sector parameters in units K0 = 1."""

    def __init__(self, m, b3=4.0 / 3.0, b4=1.0 / 3.0, beta=1.0, c=0.0, Gamma=0.0):
        self.m = float(m)
        self.m2 = self.m ** 2
        self.b3, self.b4, self.beta, self.c, self.Gamma = float(b3), float(b4), float(beta), float(c), float(Gamma)

    def rhoV(self, Y):
        return 0.75 * self.m2 * Y * Y * (1 + self.b3 * Y + self.b4 * Y * Y) + self.c * (1 + self.beta * Y)

    def F(self, Y):
        b3, b4 = self.b3, self.b4
        return Y + 1.5 * b3 * Y ** 2 + (0.5 * b3 + 2 * b4) * Y ** 3 + b4 * Y ** 4

    def dF(self, Y):
        b3, b4 = self.b3, self.b4
        return 1 + 3 * b3 * Y + 3 * (0.5 * b3 + 2 * b4) * Y ** 2 + 4 * b4 * Y ** 3

    def R(self, Y, V):
        return (3 * self.m2 * Y * (1 + 1.5 * self.b3 * Y + 2 * self.b4 * Y * Y) + 2 * self.beta * self.c
                + 3 * self.Gamma * V / (1 + Y) ** 2)

    def Q(self, Y, V):
        return 1.5 * self.Gamma * V * V / (1 + Y) ** 2

    def Y_alg(self, S):
        """Algebraic (no-derivative) quasi-static root: m^2 F(Y) - (c/3)[(4-2b)+2bY] = S/3 (Newton)."""
        c, b = self.c, self.beta
        Y = (S / 3.0 + c * (4 - 2 * b) / 3.0) / self.m2
        for _ in range(60):
            g = self.m2 * self.F(Y) - c * ((4 - 2 * b) + 2 * b * Y) / 3.0 - S / 3.0
            dg = self.m2 * self.dF(Y) - 2 * b * c / 3.0
            dY = g / dg
            Y -= dY
            if abs(dY) <= 1e-17 * max(abs(Y), 1e-300):
                break
        return Y

    def H_from_constraint(self, Y, V, rho_m):
        A = 3.0 * (1 + Y)
        B = 3.0 * V
        Cc = -(rho_m + self.rhoV(Y))
        disc = B * B - 4 * A * Cc
        return (-B + math.sqrt(disc)) / (2 * A)

    def constraint(self, Y, V, H, rho_m):
        C = 3 * H * H * (1 + Y) + 3 * H * V - rho_m - self.rhoV(Y)
        scale = abs(rho_m) + abs(self.rhoV(Y)) + 3 * H * H * abs(1 + Y) + abs(3 * H * V)
        return C, C / scale


# ------------------------------------------------------------------------------------------------ matter models
class ToyMatter:
    """Radiation + dust + a supplied matter-vacuum step DeltaV(N) (tanh of width w in N = ln a, depth Delta).
    Under P_q the final matter vacuum is absorbed into eps_inf, so DeltaV -> 0 after the step. Latent heat goes
    into radiation: d rho_r/dN = -4 rho_r - dDeltaV/dN. Supplied history [P]."""
    nz = 2

    def __init__(self, Delta, Nc, w, extra_const=0.0):
        self.Delta, self.Nc, self.w, self.const = float(Delta), float(Nc), float(w), float(extra_const)

    def dV(self, N, k=0):
        """k-th N-derivative of DeltaV(N) = Delta (1 - tanh x)/2, x = (N - Nc)/w."""
        x = (N - self.Nc) / self.w
        if x > 300:
            return 0.0
        if x < -300:
            return self.Delta if k == 0 else 0.0
        t = math.tanh(x)
        s2 = 1 - t * t
        D, w = self.Delta, self.w
        if k == 0:
            return D * 0.5 * (1 - t)
        if k == 1:
            return -D * 0.5 * s2 / w
        if k == 2:
            return D * t * s2 / w ** 2
        if k == 3:
            return D * s2 * (1 - 3 * t * t) / w ** 3
        raise ValueError

    def rho(self, z, N):
        return z[0] + z[1] + self.dV(N) + self.const

    def S(self, z, N):
        return z[1] + 4 * self.dV(N) + 4 * self.const

    def S_derivs(self, z, N):
        """S, S_N, S_NN, S_NNN along the matter history (H-independent)."""
        rd = z[1]
        return (rd + 4 * self.dV(N) + 4 * self.const, -3 * rd + 4 * self.dV(N, 1), 9 * rd + 4 * self.dV(N, 2),
                -27 * rd + 4 * self.dV(N, 3))

    def rhom_N(self, z, N):
        # d rho_m/dN and d^2 rho_m/dN^2 (latent heat cancels in the first derivative)
        rr, rd = z
        r1 = -4 * rr - 3 * rd
        r2 = 16 * rr + 4 * self.dV(N, 1) + 9 * rd
        return r1, r2


class KVMatter:
    """KV 0806.2805 matter: one fluid with constant w_M, plus a radiation component fed by the decay source Q."""
    nz = 2

    def __init__(self, wM):
        self.wM = float(wM)

    def rho(self, z, N):
        return z[0] + z[1]

    def S(self, z, N):
        return (1 - 3 * self.wM) * z[0]


class SMMatter:
    """Standard Model plasma with a tabulated equation of state w(ln rho) (Saikawa-Shirai 1803.01038 tables,
    p(T = 0) = 0 normalisation), plus dust. State z = (x = ln rho_pl, rho_d)."""
    nz = 2

    def __init__(self, lnrho_grid, w_grid):
        self.xg = np.ascontiguousarray(lnrho_grid, dtype=float)
        self.wg = np.ascontiguousarray(w_grid, dtype=float)

    def w(self, x):
        return float(np.interp(x, self.xg, self.wg))

    def rho(self, z, N):
        return math.exp(z[0]) + z[1]

    def S(self, z, N):
        e = math.exp(z[0])
        return e * (1 - 3 * self.w(z[0])) + z[1]


# ------------------------------------------------------------------------------------------------ RHS builders
def make_rhs(q, mat):
    """Full system y = (Y, V, H, N, J, z...)."""
    m2, b3, b4, beta, c, G = q.m2, q.b3, q.b4, q.beta, q.c, q.Gamma
    k3 = 0.5 * b3 + 2 * b4

    if isinstance(mat, ToyMatter):
        def rhs(t, y):
            Y, V, H, N, J, rr, rd = y
            Fy = Y + 1.5 * b3 * Y * Y + k3 * Y ** 3 + b4 * Y ** 4
            Rr = 3 * m2 * Y * (1 + 1.5 * b3 * Y + 2 * b4 * Y * Y) + 2 * beta * c
            Qd = 0.0
            if G:
                Rr += 3 * G * V / (1 + Y) ** 2
                Qd = 1.5 * G * V * V / (1 + Y) ** 2
            S = rd + 4 * mat.dV(N) + 4 * mat.const
            Vd = S / 3.0 + (c / 3.0) * ((4 - 2 * beta) + 2 * beta * Y) - m2 * Fy - 3 * H * V
            if G:
                Vd -= G * V / (1 + Y)
            return [V, Vd, Rr / 6.0 - 2 * H * H, H, H * H,
                    -4 * H * rr - H * mat.dV(N, 1) + Qd, -3 * H * rd]
        return rhs

    if isinstance(mat, KVMatter):
        wM = mat.wM

        def rhs(t, y):
            Y, V, H, N, J, rM, rr = y
            Fy = Y + 1.5 * b3 * Y * Y + k3 * Y ** 3 + b4 * Y ** 4
            Rr = 3 * m2 * Y * (1 + 1.5 * b3 * Y + 2 * b4 * Y * Y) + 2 * beta * c + 3 * G * V / (1 + Y) ** 2
            Qd = 1.5 * G * V * V / (1 + Y) ** 2
            S = (1 - 3 * wM) * rM
            Vd = S / 3.0 + (c / 3.0) * ((4 - 2 * beta) + 2 * beta * Y) - m2 * Fy - 3 * H * V - G * V / (1 + Y)
            return [V, Vd, Rr / 6.0 - 2 * H * H, H, H * H, -3 * H * (1 + wM) * rM, -4 * H * rr + Qd]
        return rhs

    if isinstance(mat, SMMatter):
        xg, wg = mat.xg, mat.wg
        interp = np.interp

        def rhs(t, y):
            Y, V, H, N, J, x, rd = y
            wv = float(interp(x, xg, wg))
            e = math.exp(x)
            Fy = Y + 1.5 * b3 * Y * Y + k3 * Y ** 3 + b4 * Y ** 4
            Rr = 3 * m2 * Y * (1 + 1.5 * b3 * Y + 2 * b4 * Y * Y) + 2 * beta * c
            S = e * (1 - 3 * wv) + rd
            Vd = S / 3.0 + (c / 3.0) * ((4 - 2 * beta) + 2 * beta * Y) - m2 * Fy - 3 * H * V
            return [V, Vd, Rr / 6.0 - 2 * H * H, H, H * H, -3 * H * (1 + wv), -3 * H * rd]
        return rhs
    raise TypeError(mat)


def integrate_full(q, mat, y0, t0, t_out, rtol=1e-12, atol=1e-30, scale=None, nsteps=10 ** 9, counter=None):
    """Integrate with the Fortran DOP853 (scipy.integrate.ode) and return y at the sorted times t_out.
    The Fortran wrapper takes a scalar atol only, so components are rescaled by `scale` (u = y/scale) and the
    error control is atol + rtol |u_i| on the scaled variables."""
    rhs0 = make_rhs(q, mat)
    if scale is None:
        rhs = rhs0
        u0 = list(y0)
    else:
        sc = [float(s) for s in scale]
        n = len(sc)

        def rhs(t, u):
            f = rhs0(t, [u[i] * sc[i] for i in range(n)])
            return [f[i] / sc[i] for i in range(n)]
        u0 = [y0[i] / sc[i] for i in range(n)]
    if counter is not None:
        inner = rhs

        def rhs(t, u):
            counter[0] += 1
            return inner(t, u)
    r = ode(rhs)
    r.set_integrator('dop853', rtol=rtol, atol=atol, nsteps=nsteps, verbosity=-1)
    r.set_initial_value(u0, t0)
    out = np.empty((len(t_out), len(y0)))
    for i, t in enumerate(t_out):
        r.integrate(t)
        if not r.successful():
            raise RuntimeError('integration failed at t=%g (i=%d)' % (t, i))
        out[i] = r.y
    if scale is not None:
        out *= np.array(scale)[None, :]
    return out


# ------------------------------------------------------------------------------------------------ reductions
def reduced_toy(mat, q, z0, N0, N1, order=1, rtol=1e-13):
    """Adiabatic reduction for ToyMatter, integrated in N. order=0: GR + DeltaV(N) (the matter vacuum relative to
    c*(mu) = final vacuum under P_q).  order=1: + adiabatic q-correction
        3 H^2 (1 + (S + S_N)/(3 m^2)) = rho_m + S^2/(12 m^2)          (K0 = 1)
    State u = (rho_r, rho_d, t, J) with J = int H dN."""
    m2 = q.m2

    def Hof(N, rr, rd):
        z = (rr, rd)
        rho = mat.rho(z, N)
        if order == 0:
            return math.sqrt(rho / 3.0)
        S, SN, _, _ = mat.S_derivs(z, N)
        return math.sqrt((rho + S * S / (12 * m2)) / (3.0 * (1 + (S + SN) / (3 * m2))))

    def f(N, u):
        rr, rd, t, J = u
        H = Hof(N, rr, rd)
        return [-4 * rr - mat.dV(N, 1), -3 * rd, 1.0 / H, H]

    # three segments so that the step (width w) is resolved without forcing tiny steps elsewhere
    edges = [N0]
    for e in (mat.Nc - 12 * mat.w, mat.Nc + 12 * mat.w):
        if edges[-1] < e < N1:
            edges.append(e)
    edges.append(N1)
    u = [z0[0], z0[1], 0.0, 0.0]
    segs = []
    for a, b in zip(edges[:-1], edges[1:]):
        inside = (a >= mat.Nc - 12 * mat.w - 1e-12) and (b <= mat.Nc + 12 * mat.w + 1e-12)
        s = solve_ivp(f, (a, b), u, method='DOP853', rtol=rtol, atol=1e-30, dense_output=True,
                      max_step=(mat.w / 4 if inside else 0.02))
        segs.append((a, b, s.sol))
        u = list(s.y[:, -1])

    class _Sol:
        def sol(self, N):
            for a, b, sf in segs:
                if N <= b:
                    return sf(N)
            return segs[-1][2](N)
    return _Sol(), Hof


def quasi_static_toy(q, mat, N, rr, rd, H, Hdot=None):
    """Quasi-static Y_qs and V_qs to O(m^-4) on the toy background, with analytic N-derivatives.
    Y1 = Y_alg(S); Y2 = -(Y1dd + 3 H Y1d)/(m^2 F'(Y1)); background H derivatives from the GR relation."""
    z = (rr, rd)
    S, S1, S2, S3 = mat.S_derivs(z, N)
    m2 = q.m2
    Y1 = q.Y_alg(S)
    dF = q.dF(Y1)
    # derivatives of Y1 in N (leading order: Y1 ~ S/(3 m^2 F'))
    G1, G2, G3 = S1 / (3 * m2 * dF), S2 / (3 * m2 * dF), S3 / (3 * m2 * dF)
    rho = mat.rho(z, N)
    r1, r2 = mat.rhom_N(z, N)
    HN = r1 / (6 * H)
    HNN = (r2 - 6 * HN * HN) / (6 * H)
    # time derivatives: d/dt = H d/dN
    Y1d = H * G1
    Y1dd = H * (HN * G1 + H * G2)
    Y1ddd = H * (HN * (HN * G1 + H * G2) + H * (HNN * G1 + HN * G2 + HN * G2 + H * G3))
    Hd = H * HN
    Y2 = -(Y1dd + 3 * H * Y1d) / (m2 * dF)
    Y2d = -(Y1ddd + 3 * Hd * Y1d + 3 * H * Y1dd) / (m2 * dF)
    return Y1 + Y2, Y1d + Y2d, Y1, Y1d


def osc_energy(q, Y, V, Yqs, Vqs):
    """Einstein-frame energy density of the q-oscillation about the quasi-static track (units K0 = 1):
    sigma = sqrt(3/2) M ln(K/K0) ~ sqrt(3/2) Y  =>  E = (3/4) [ (V - Vqs)^2 + m^2 (Y - Yqs)^2 ]."""
    return 0.75 * ((V - Vqs) ** 2 + q.m2 * (Y - Yqs) ** 2)


def window_grid(t_start, t_end, period, n_windows, n_sub=16, log=True):
    """Window centres (log- or lin-spaced) each carrying n_sub+1 equally spaced points across one period."""
    lo, hi = t_start + period, t_end - period
    if log:
        centres = np.exp(np.linspace(math.log(lo), math.log(hi), n_windows))
    else:
        centres = np.linspace(lo, hi, n_windows)
    pts = []
    for c in centres:
        pts.append(c + period * (np.arange(n_sub + 1) / n_sub - 0.5))
    pts = np.array(pts)
    # remove overlaps by keeping windows disjoint
    keep = [0]
    for i in range(1, len(centres)):
        if pts[i, 0] > pts[keep[-1], -1]:
            keep.append(i)
    return centres[keep], pts[keep]
