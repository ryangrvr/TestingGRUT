"""Verifier v_vacua, part 1: q-sector claims of track `vacua` (NOTES.md sections 1.1-1.5, 3.4).

Independent code (sympy/numpy/scipy). Checks:
  1a  de Sitter reduction of the field equations (signs of E1, E2) from an explicit Einstein tensor
  1b  series Lambda(x) = x + a2 x^2 + a3 x^3 at a shared mu, q-shift, dK/K, and H^2 seen with K0 (no x^2 term?)
  1c  exact slope dLambda/dV along the static family and its pole
  1d  numerical static roots for KV-type eps (track: 0.010190 at K'=1, x=0.01; 0.0100177 at K'=0.3)
  1e  q-coupled matter (theta F~): Lambda_i = V_i - eps*(mu - theta_i); coexistence / Maxwell point
  1f  junction conditions: delta-function parts of G_ab and R for a codimension-1 wall, trace condition,
      and a linear thick-wall solve (kink -> -sigma for m w << 1, local tracking for m w >> 1)
  1g  size estimates (a2 x at Planck depth, xi bound, |H|^2-coupling shift)
Run: python3 v1_qsector.py  -> ../results/v1_qsector.json
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.optimize import fsolve
from scipy.integrate import solve_bvp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'v1_qsector.json')
res = {}

# ------------------------------------------------------------------ 1a: Einstein tensor of de Sitter, signs
t, x, y, z, H = sp.symbols('t x y z H', positive=True)
coords = [t, x, y, z]
a = sp.exp(H * t)
g = sp.diag(-1, a ** 2, a ** 2, a ** 2)


def einstein(g, coords):
    n = len(coords)
    gi = g.inv()
    Gam = [[[sp.simplify(sum(gi[l, s] * (sp.diff(g[s, m], coords[nn]) + sp.diff(g[s, nn], coords[m])
                                         - sp.diff(g[m, nn], coords[s])) for s in range(n)) / 2)
             for nn in range(n)] for m in range(n)] for l in range(n)]
    Ric = sp.zeros(n, n)
    for m in range(n):
        for nn in range(n):
            Ric[m, nn] = sp.simplify(sum(sp.diff(Gam[l][m][nn], coords[l]) for l in range(n))
                                     - sum(sp.diff(Gam[l][m][l], coords[nn]) for l in range(n))
                                     + sum(Gam[l][l][s] * Gam[s][m][nn] for l in range(n) for s in range(n))
                                     - sum(Gam[l][nn][s] * Gam[s][m][l] for l in range(n) for s in range(n)))
    Rs = sp.simplify(sum(gi[m, nn] * Ric[m, nn] for m in range(n) for nn in range(n)))
    G = sp.simplify(Ric - Rs * g / 2)
    return G, Rs


G_dS, R_dS = einstein(g, coords)
# (g): K G = T - (eps - mu q) g with T = -V g  =>  K G = -(V + eps - mu q) g =: -Lam g
# check G = -3 H^2 g  and R = 12 H^2
res['1a_dS'] = {'G_plus_3H2_g_is_zero': bool(sp.simplify(G_dS + 3 * H ** 2 * g) == sp.zeros(4, 4)),
                'R': str(R_dS),
                'reading': 'K G = -(V + eps - mu q) g with G = -3H^2 g gives 3 K H^2 = Lam = V + eps - mu q (E2); '
                           '(q) eps\' - K\' R/2 = mu with R = 12 H^2 = 4 Lam/K gives E1: eps\' - 2 K\' Lam/K = mu'}

# ------------------------------------------------------------------ 1b: series at a shared mu
X = sp.symbols('x')
K0, K1, K2, K3 = sp.symbols('K0 K1 K2 K3')
e2, e3, e4 = sp.symbols('e2 e3 e4')
mu, q0 = sp.symbols('mu q0')
d = sp.symbols('delta')
# Taylor data about q0, where eps'(q0) = mu (so q0 = q_mu is the K'=0 / Minkowski root for this mu)
eps0 = sp.symbols('eps0')
eps = eps0 + mu * d + e2 * d ** 2 / 2 + e3 * d ** 3 / 6 + e4 * d ** 4 / 24
epsp = sp.diff(eps, d)
Kq = K0 + K1 * d + K2 * d ** 2 / 2 + K3 * d ** 3 / 6
Kp = sp.diff(Kq, d)
# eps*(mu) = mu q0 - eps(q0) = mu q0 - eps0 ; x = V - eps*(mu) => V = x + mu q0 - eps0
V = X + mu * q0 - eps0
# unknowns: delta(x) = c1 x + c2 x^2 + c3 x^3, Lam(x) = x + b2 x^2 + b3 x^3
c1, c2, c3, b2, b3 = sp.symbols('c1 c2 c3 b2 b3')
dser = c1 * X + c2 * X ** 2 + c3 * X ** 3
Lser = X + b2 * X ** 2 + b3 * X ** 3
E1 = sp.expand(sp.series((epsp - 2 * Kp * Lser / Kq - mu).subs(d, dser), X, 0, 4).removeO())
# E2: Lam = V + eps(q) - mu q, with q = q0 + delta
E2 = sp.expand(sp.series((V + eps - mu * (q0 + d) - Lser).subs(d, dser), X, 0, 4).removeO())
eqs = []
for p in (1, 2, 3):
    eqs += [E1.coeff(X, p), E2.coeff(X, p)]
sol = sp.solve(eqs, [c1, c2, c3, b2, b3], dict=True)[0]
a2 = sp.simplify(sol[b2])
a3 = sp.simplify(sol[b3])
mq2 = sp.Rational(2, 3) * K0 * e2 / K1 ** 2
track_a2 = 2 * K1 ** 2 / (K0 ** 2 * e2)
track_a3 = 8 * K1 ** 2 * (3 * K2 * e2 - K1 * e3) / (3 * K0 ** 3 * e2 ** 3)
# H^2 = Lam / (3 K(q)) expanded, versus x/(3 K0)
Hs2 = sp.series((Lser / (3 * Kq.subs(d, dser))).subs(sol), X, 0, 4).removeO()
res['1b_series'] = {
    'a2': str(a2), 'a3': str(a3),
    'a2_equals_track': bool(sp.simplify(a2 - track_a2) == 0),
    'a2_equals_4/(3 K0 mq^2)': bool(sp.simplify(a2 - 4 / (3 * K0 * mq2)) == 0),
    'a3_equals_track': bool(sp.simplify(a3 - track_a3) == 0),
    'dq/dx': str(sp.simplify(sol[c1])),
    'dK/K per x (first order)': str(sp.simplify(K1 * sol[c1] / K0)),
    'dK/K equals a2': bool(sp.simplify(K1 * sol[c1] / K0 - a2) == 0),
    'H^2 = Lam/(3K(q)) series': str(sp.simplify(sp.expand(Hs2))),
    'H^2 x^2 coefficient': str(sp.simplify(sp.expand(Hs2).coeff(X, 2))),
}

# ------------------------------------------------------------------ 1c: exact slope along static family
qq, Lam, Vv, mm = sp.symbols('q Lam V mu')
ef = sp.Function('eps')(qq)
Kf = sp.Function('K')(qq)
F1 = sp.diff(ef, qq) - 2 * sp.diff(Kf, qq) * Lam / Kf - mm
F2 = Vv + ef - mm * qq - Lam
# implicit: d/dV of (F1, F2) = 0 with q(V), Lam(V)
J = sp.Matrix([[sp.diff(F1, qq), sp.diff(F1, Lam)], [sp.diff(F2, qq), sp.diff(F2, Lam)]])
rhs = -sp.Matrix([sp.diff(F1, Vv), sp.diff(F2, Vv)])
dsol = sp.simplify(J.LUsolve(rhs))
dLdV = sp.simplify(dsol[1])
Kp_, Kpp_, epp_ = sp.diff(Kf, qq), sp.diff(Kf, qq, 2), sp.diff(ef, qq, 2)
D = epp_ - 2 * Lam * (Kf * Kpp_ - Kp_ ** 2) / Kf ** 2
track_slope = Kf ** 2 * D / (Kf ** 2 * D - 4 * Lam * Kp_ ** 2)
dLdV_onshell = sp.simplify(dLdV.subs(mm, sp.diff(ef, qq) - 2 * sp.diff(Kf, qq) * Lam / Kf))
res['1c_slope'] = {'dLam/dV': str(dLdV), 'dLam/dV_on_shell(E1 used)': str(dLdV_onshell),
                   'equals_track_formula_on_shell': bool(sp.simplify(dLdV_onshell - track_slope) == 0),
                   'pole_at_Lam': 'K^2 D/(4 K\'^2) = (3/8) K m_q^2 at leading order (K\'\'=0 or Lam small)',
                   'check_pole_leading': bool(sp.simplify((K0 ** 2 * e2 / (4 * K1 ** 2)) - sp.Rational(3, 8) * K0 * mq2) == 0)}

# ------------------------------------------------------------------ 1d: numerical static roots, KV-type eps
def static_root(Kp1, xv):
    epsf = lambda q: 0.5 * (-q ** 2 + q ** 4 / 3)
    epsp_f = lambda q: -q + 2 * q ** 3 / 3
    Kf_ = lambda q: 1 + Kp1 * (q - 1)
    mu0 = epsp_f(1.0)              # -1/3 ; eps(1) - mu0 = 0 -> Minkowski equilibrium at q0 = 1
    Vx = xv                        # eps*(mu0) = mu0*1 - eps(1) = 0
    f = lambda u: [epsp_f(u[0]) - 2 * Kp1 * u[1] / Kf_(u[0]) - mu0, Vx + epsf(u[0]) - mu0 * u[0] - u[1]]
    u = fsolve(f, [1.0, xv], xtol=1e-13)
    return float(u[0] - 1), float(u[1]), [float(v) for v in f(u)]


rows = []
for Kp1 in (0.0, 0.3, 1.0):
    for xv in (-0.001, 0.001, 0.01, 0.05):
        dq, L, resid = static_root(Kp1, xv)
        # KV eps: e2 = 1, e3 = 4 at q = 1, K0 = 1, K2 = 0
        a2n, a3n = 2 * Kp1 ** 2, 8 * Kp1 ** 2 * (0 - Kp1 * 4) / 3
        rows.append({'Kprime': Kp1, 'x': xv, 'Lam_exact': L, 'two_term': xv + a2n * xv ** 2,
                     'three_term': xv + a2n * xv ** 2 + a3n * xv ** 3, 'q-q0': dq, 'residual': resid})
res['1d_numeric_KV'] = rows

# ------------------------------------------------------------------ 1e: q-coupled matter, coexistence
qv, th, e2s, qc, Vi, Vj, thi, thj, mus = sp.symbols('q theta e2 q_c V_i V_j theta_i theta_j mu', real=True)
# quadratic eps about q_c with eps'(q_c) = 0 (generic shift absorbed into mu): eps = e2 (q - q_c)^2 / 2
epsq = e2s * (qv - qc) ** 2 / 2


def Lam_coupled(Vk, thk):
    # A-equation with matter term -theta q: eps'(q) + theta = mu (K' = 0)
    qs = sp.solve(sp.Eq(sp.diff(epsq, qv) + thk, mus), qv)[0]
    # gravitating: V + eps - q eps'(q) = V + eps - q (mu - theta)
    return sp.simplify(Vk + epsq.subs(qv, qs) - qs * (mus - thk)), qs


Li, qi = Lam_coupled(Vi, thi)
Lj, qj = Lam_coupled(Vj, thj)
nu = sp.symbols('nu')
q_nu = sp.solve(sp.Eq(sp.diff(epsq, qv), nu), qv)[0]
eps_star = sp.simplify(nu * q_nu - epsq.subs(qv, q_nu))
res['1e_coupled'] = {
    'Lam_i - [V_i - eps*(mu - theta_i)]': str(sp.simplify(Li - (Vi - eps_star.subs(nu, mus - thi)))),
    'Lam_i - Lam_j (exact, quadratic eps)': str(sp.factor(sp.simplify(Li - Lj))),
    'first-order term check (should be 0)': str(sp.simplify(sp.diff(Li - Lj, thi).subs({thi: 0, thj: 0})
                                                             - q_nu.subs(nu, mus))),
    'q_i - q_j': str(sp.simplify(qi - qj)),
    'P_i - P_j at theta = 0': str(sp.simplify((-(Li - Lj)).subs({thi: 0, thj: 0}))),
    'mu_c (P_i = P_j) with theta': str(sp.solve(sp.Eq(Li, Lj), mus)),
}

# ------------------------------------------------------------------ 1f: junction conditions
n_ = sp.symbols('n', real=True)
A = sp.Function('A')(n_)
coords2 = [t, x, y, n_]
gw = sp.diag(-A ** 2, A ** 2, A ** 2, 1)        # ds^2 = A(n)^2 eta_ab dx^a dx^b + dn^2 (pure-tension symmetric wall)
Gw, Rw = einstein(gw, coords2)
App = sp.diff(A, n_, 2)
# distributional part: coefficient of A'' (A' jumps across the wall -> A'' contains [A'] delta(n))
cG_tt = sp.simplify(sp.diff(Gw[0, 0], App))    # G_tt coefficient of A''
cG_xx = sp.simplify(sp.diff(Gw[1, 1], App))
cR = sp.simplify(sp.diff(Rw, App))
cG_nn = sp.simplify(sp.diff(Gw[3, 3], App))
# expectation with Kx_ab = (1/2) d_n h_ab = (A'/A) h_ab, Kx = 3 A'/A:
#   G_ab |delta = -([Kx_ab] - h_ab [Kx]) = 2 [A'/A] h_ab  -> coefficient of A'' in G_ab = 2 h_ab / A
#   R |delta = -2 [Kx] = -6 [A'/A]                     -> coefficient of A'' in R   = -6 / A
res['1f_junction'] = {
    'coef_App_in_G_tt': str(cG_tt), 'expected_2 h_tt/A': str(sp.simplify(2 * (-A ** 2) / A)),
    'coef_App_in_G_xx': str(cG_xx), 'expected_2 h_xx/A': str(sp.simplify(2 * A ** 2 / A)),
    'coef_App_in_G_nn': str(cG_nn), 'coef_App_in_R': str(cR), 'expected_-6/A': str(-6 / A),
    'consequence': '(g) delta part: -K([Kx_ab] - h_ab [Kx]) + h_ab [d_n K] = S_ab; (q): R has no delta (q continuous, '
                   'K\' != 0) => [Kx] = 0; trace => [d_n K] = S/3; traceless => -K [Kx_ab] = S_ab - h_ab S/3. '
                   'Matches track 1.5 item 3.',
}
G_ok = (sp.simplify(cG_tt - 2 * (-A ** 2) / A) == 0) and (sp.simplify(cG_xx - 2 * A) == 0) and \
       (sp.simplify(cR + 6 / A) == 0) and (cG_nn == 0)
res['1f_junction']['distributional_parts_confirmed'] = bool(G_ok)


# linear thick wall: dK'' - m^2 dK = -rho(z), rho = sigma/(2w) sech^2(z/w) (unit integral sigma)
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


def wall(mw, sigma=1.0, w=1.0):
    m = mw / w
    L = max(40 * w, 40 / m)
    N = 400001
    zz = np.linspace(-L, L, N)
    h = zz[1] - zz[0]
    rho = sigma / (2 * w) / np.cosh(np.clip(zz / w, -300, 300)) ** 2
    Aop = diags([np.ones(N - 1) / h ** 2, (-2 / h ** 2 - m ** 2) * np.ones(N), np.ones(N - 1) / h ** 2], [-1, 0, 1],
                format='csc')
    dK = spsolve(Aop, -rho)              # dK = 0 at |z| = L (e^{-m L} negligible)
    dKp = np.gradient(dK, h)
    zo = 3 * w
    i1, i2 = np.argmin(np.abs(zz + zo)), np.argmin(np.abs(zz - zo))
    jump = dKp[i2] - dKp[i1]
    core = np.abs(zz) < 2 * w
    track_dev = np.max(np.abs(dK - rho / m ** 2)[core]) / np.max(rho / m ** 2)
    return {'m_w': mw, 'kink_jump_over_sigma(+-3w)': float(jump / sigma),
            'kink_corrected_exp(m zo)': float(jump * np.exp(m * zo) / sigma) if mw <= 0.1 else None,
            'max_rel_dev_from_local_tracking(|z|<2w)': float(track_dev),
            'analytic_center_dev 2/(m w)^2 (sech^2: rho\'\'/rho = -2/w^2)': 2 / mw ** 2}


res['1f_linear_wall_sech2'] = [wall(mw) for mw in (0.01, 0.1, 1.0, 10.0, 30.0)]

# ------------------------------------------------------------------ 1g: size estimates
M = 2.435e18
K = M ** 2
mq_nat = np.sqrt(2 / 3) * M
rhoL = 2.52e-47
mh = 125.2
a2x = lambda xv, mq=mq_nat: 4 / 3 * xv / (K * mq ** 2)
res['1g_sizes'] = {
    'm_q_natural': mq_nat,
    'a2 x at x=7.5e70 (track: ~4e-3)': a2x(7.5e70),
    'a2 x at x=1e71 (track says <0.5% down to 1e71)': a2x(1e71),
    'x at which a2 x = 0.5%': 0.005 * K * mq_nat ** 2 * 3 / 4,
    'a2 x at x=1e8 (track ~6e-66)': a2x(1e8),
    'xi_max = mh^2/M^2 (track ~3e-33)': mh ** 2 / M ** 2,
    'mh^2 h*^2 at h*=4.7e17 and 2.435e18 (track ~1e40)': [mh ** 2 * 4.7e17 ** 2, mh ** 2 * M ** 2],
    'mq w for EW wall w=1/mh (track 1.6e16)': mq_nat / mh,
    'mq R for R=1/(1e17 GeV) (track ~20)': mq_nat / 1e17,
    'tolerance delta lambda = 4 rho_L/h*^4 at h* = 4.7e17, 2.435e18, 1.22e19 (track ~1e-118)':
        [4 * rhoL / 4.7e17 ** 4, 4 * rhoL / M ** 4, 4 * rhoL / 1.2209e19 ** 4],
}

with open(OUT, 'w') as f:
    json.dump(res, f, indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
