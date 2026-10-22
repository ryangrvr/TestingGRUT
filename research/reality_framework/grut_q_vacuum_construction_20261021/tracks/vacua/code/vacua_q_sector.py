"""q-sector with several matter vacua (track `vacua`, NOTES.md sections 1-2).

Action (BRIEF.md): S = int sqrt(-g)[ K(q) R/2 - eps(q) ] + S_m[phi, g];  (q) eps' - K' R/2 = mu (global constant).
Matter vacuum i: constant phi_i with q-independent energy V_i (full T = 0 effective potential at its minimum).
Static maximally symmetric state: R = 4 Lam/K, H^2 = Lam/(3K), and
    E1: eps'(q) - 2 K'(q) Lam / K(q) - mu = 0,      E2: Lam - V - eps(q) + mu q = 0.
Checks (all printed and written to ../results/vacua_q_sector.json):
 A. series Lam_i = x + a2 x^2 + a3 x^3, x = V_i - eps*(mu) (Legendre transform), generic Taylor data;
    exact slope dLam/dV; numerical solution of E1, E2 for KV 0806.2805-type eps and linear K.
 B. q-dependent matter energy V(phi_i, q) = V_i + theta_i q: Lam_i = V_i - eps*(mu - theta_i) (K' = 0).
 C. coexistence / Maxwell construction at shared mu; multi-4-form rank argument (Volovik hep-ph/0309144 count).
 D. linearised planar wall: (d^2/dz^2 - m^2) dK = T/3 with T = -3 rho_w(z); kink [dK'] -> -sigma for m w << 1,
    local tracking dK -> rho_w/m^2 for m w >> 1.
"""
import json, os
import numpy as np
import sympy as sp
from scipy.optimize import fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'vacua_q_sector.json')
res = {}

# ------------------------------------------------------------------------------------------------ A. series
x, d = sp.symbols('x delta')
K0, K1, K2, K3, e2, e3, e4 = sp.symbols('K0 K1 K2 K3 e2 e3 e4', real=True)
# Taylor data at the Minkowski point q0 of this mu (eps'(q0) = mu, and V_star = eps*(mu) = mu q0 - eps(q0)).
# Write q = q0 + d, eps(q) - mu q = -V_star + e2 d^2/2 + e3 d^3/6 + e4 d^4/24; K = K0 + K1 d + K2 d^2/2 + K3 d^3/6.
U = e2 * d ** 2 / 2 + e3 * d ** 3 / 6 + e4 * d ** 4 / 24          # eps - mu q + V_star
Kq = K0 + K1 * d + K2 * d ** 2 / 2 + K3 * d ** 3 / 6
Lam_s = sp.symbols('Lam')
# E2: Lam = V + eps - mu q = x + U(d) with x = V - V_star ; E1: U'(d) = 2 K'(d) Lam / K(d)
c1, c2, c3 = sp.symbols('c1 c2 c3')
dser = c1 * x + c2 * x ** 2 + c3 * x ** 3
Lam_expr = x + U.subs(d, dser)
E1 = sp.diff(U, d).subs(d, dser) * Kq.subs(d, dser) - 2 * sp.diff(Kq, d).subs(d, dser) * Lam_expr
E1s = sp.series(sp.expand(E1), x, 0, 4).removeO()
sol = sp.solve([sp.expand(E1s).coeff(x, k) for k in (1, 2, 3)], [c1, c2, c3], dict=True)[0]
Lam_ser = sp.series(Lam_expr.subs(sol), x, 0, 4).removeO()
a2 = sp.simplify(Lam_ser.coeff(x, 2))
a3 = sp.simplify(Lam_ser.coeff(x, 3))
mq2 = sp.Rational(2, 3) * K0 * e2 / K1 ** 2
res['A_series'] = {
    'Lam = x + a2 x^2 + a3 x^3, x = V_i - eps*(mu)': True,
    'a2': str(a2), 'a2_minus_4/(3 K0 mq^2)': str(sp.simplify(a2 - sp.Rational(4, 3) / (K0 * mq2))),
    'a3': str(sp.factor(a3)),
    'dq/dx (= q_i - q0 per unit x)': str(sp.simplify(sol[c1])),
    'Delta Lam_ij': '(V_i - V_j) [1 + a2 (x_i + x_j) + O(x^2)]',
}
print('A: a2 =', a2, ' check a2 - 4/(3 K0 mq^2) =', sp.simplify(a2 - sp.Rational(4, 3) / (K0 * mq2)))
print('A: a3 =', sp.factor(a3))
print('A: dq/dx =', sp.simplify(sol[c1]))

# exact slope dLam/dV along the static family (implicit differentiation of E1, E2)
q, L, mu, V = sp.symbols('q Lam mu V')
eps_f, K_f = sp.Function('eps'), sp.Function('K')
F1 = sp.diff(eps_f(q), q) - 2 * sp.diff(K_f(q), q) * L / K_f(q) - mu
F2 = L - V - eps_f(q) + mu * q
J = sp.Matrix([[sp.diff(F1, q), sp.diff(F1, L)], [sp.diff(F2, q), sp.diff(F2, L)]])
dF_dV = sp.Matrix([sp.diff(F1, V), sp.diff(F2, V)])
dqL = -J.inv() * dF_dV
dLdV = sp.simplify(dqL[1].subs(mu, sp.diff(eps_f(q), q) - 2 * sp.diff(K_f(q), q) * L / K_f(q)))
res['A_exact_slope'] = {'dLam/dV': str(dLdV),
                        'reading': '= 1/(1 - 4 K\'^2 Lam/(K^2 D)), D = eps\'\' - 2 Lam (K K\'\' - K\'^2)/K^2: equals 1 '
                                   'at Lam = 0 and stays positive for Lam < K^2 D/(4 K\'^2) ~ (3/8) K m_q^2'}
print('A: dLam/dV =', dLdV)

# numerical: KV 0806.2805 eps (their eq. 4.7, units chi = q0 = 1 up to sign conventions) and linear K
def eps_n(qq):   # eps = (1/2)(-q^2 + q^4/3) + const, stable branch q ~ 1
    return 0.5 * (-qq ** 2 + qq ** 4 / 3.0)
def deps_n(qq):
    return 0.5 * (-2 * qq + 4 * qq ** 3 / 3.0)
def d2eps_n(qq):
    return 0.5 * (-2 + 4 * qq ** 2)
def K_n(qq, s):
    return 1.0 + s * (qq - 1.0)          # K(q0 = 1) = 1, K' = s
q0n = 1.0
mu_n = deps_n(q0n)                       # mu fixed by "Minkowski for V = V_star"
V_star = mu_n * q0n - eps_n(q0n)         # eps*(mu)
rows = []
for s in (0.0, 0.3, 1.0):
    for xv in (-1e-3, 1e-3, 1e-2, 5e-2):
        Vv = V_star + xv
        f = lambda z: [deps_n(z[0]) - 2 * s * z[1] / K_n(z[0], s) - mu_n, z[1] - Vv - eps_n(z[0]) + mu_n * z[0]]
        qs, Ls = fsolve(f, [q0n, xv], xtol=1e-14)
        mq2n = (2.0 / 3.0) * K_n(q0n, s) * d2eps_n(q0n) / s ** 2 if s else np.inf
        series2 = xv + (4.0 / 3.0) * xv ** 2 / mq2n
        rows.append({'Kprime': s, 'x': xv, 'Lam_exact': float(Ls), 'Lam_series2': float(series2),
                     'q_s - q0': float(qs - q0n), 'rel_dev_from_x': float((Ls - xv) / xv)})
res['A_numeric_KV_eps'] = {'eps': '(1/2)(-q^2 + q^4/3), q0 = 1, eps\'\'(q0) = 1, K = 1 + K\'(q-1)', 'rows': rows}
for r_ in rows:
    print('A num:', r_)

# ------------------------------------------------------------------------------------------------ B. q-coupled matter
th, e2s, qc, mus, Vi = sp.symbols('theta e2 q_c mu V_i', real=True)
eps_quad = e2s / 2 * (q - qc) ** 2
# E_i(q) = eps(q) + V_i + theta_i q ; Legendre: E_i*(mu) = sup_q [mu q - E_i(q)] ; Lam_i = -E_i*(mu)
qstar = sp.solve(sp.diff(mus * q - eps_quad - Vi - th * q, q), q)[0]
Estar = sp.simplify((mus * q - eps_quad - Vi - th * q).subs(q, qstar))
epsstar = sp.simplify((mus * q - eps_quad).subs(q, sp.solve(sp.diff(mus * q - eps_quad, q), q)[0]))
check = sp.simplify(-Estar - (Vi - epsstar.subs(mus, mus - th)))
res['B_q_coupled'] = {'coupling': 'S_m contains -int sqrt(-g) theta(phi) q  (= -int theta(phi) F~, metric independent)',
                      'Lam_i': 'V_i - eps*(mu - theta_i)', 'sympy_check_quadratic_eps (should be 0)': str(check),
                      'q_i': str(qstar),
                      'Delta Lam_ij first order': '(V_i - V_j) - q_mu (theta_j - theta_i) ... = (V_i - V_j) + q_mu (theta_i - theta_j)'}
print('B: Lam_i = V_i - eps*(mu - theta_i) check:', check, ' q_i =', qstar)

# ------------------------------------------------------------------------------------------------ C. coexistence
# no coupling, K' = 0: P_i(mu) = -Lam_i(mu) = eps*(mu) - V_i  ->  P_i - P_j = V_j - V_i for every mu.
Vj, thj = sp.symbols('V_j theta_j', real=True)
Pi = -(Vi - epsstar)
Pj = -(Vj - epsstar)
res['C_coexistence_no_coupling'] = {'P_i - P_j': str(sp.simplify(Pi - Pj)),
                                    'verdict': 'independent of mu: equal pressure (coexistence) iff V_i = V_j'}
# with coupling: Maxwell point
Pi_c = -(Vi - epsstar.subs(mus, mus - th))
Pj_c = -(Vj - epsstar.subs(mus, mus - thj))
mu_c = sp.solve(sp.Eq(Pi_c, Pj_c), mus)
res['C_coexistence_with_coupling'] = {'mu_c': [str(sp.simplify(m)) for m in mu_c],
                                      'P at mu_c': [str(sp.simplify(Pi_c.subs(mus, m))) for m in mu_c],
                                      'verdict': 'a coexistence value mu_c exists for any V_i, V_j if theta_i != '
                                                 'theta_j; requiring also P(mu_c) = 0 is one further relation among '
                                                 '(V_i, V_j, theta_i, theta_j, eps data)'}
print('C: mu_c =', mu_c)
# multi-4-form rank argument: P_i(mu_1..mu_k) = eps*(mu_1..mu_k) - V_i ; dP_i/dmu_a = q_a (same for all i)
res['C_multi_form'] = ('k 4-forms, matter vacua without q coupling: P_i(mu_a) = eps*(mu_a) - V_i, Jacobian rows '
                       'dP_i/dmu_a = q_a identical for all i (rank 1) -> at most one non-degenerate vacuum at P = 0; '
                       'Volovik\'s nu <= k needs the phases to differ in their charge densities, i.e. q-coupled matter')

# ------------------------------------------------------------------------------------------------ D. planar wall
def wall(mw, m=1.0, sigma=1.0, N=20001, Lbox=None):
    """Solve dK'' - m^2 dK = -rho(z) with rho = sigma/(2w) sech^2(z/w) (int rho = sigma)."""
    w = mw / m
    Lbox = Lbox or 40.0 * max(w, 1.0 / m)
    z = np.linspace(-Lbox, Lbox, N)
    ex = np.exp(-2 * np.abs(z / w))
    rho = sigma / (2 * w) * 4 * ex / (1 + ex) ** 2          # sigma/(2w) sech^2(z/w), overflow-safe
    # Green's function solution: dK(z) = int G(z-z') rho(z') dz', G = exp(-m|z|)/(2m)
    dz = z[1] - z[0]
    G = np.exp(-m * np.abs(z - z[:, None])) / (2 * m)
    dK = G @ rho * dz
    dKp = np.gradient(dK, z)
    # effective kink: dK' just outside the wall region (|z| = 3w) on each side
    zo = min(3 * w, 0.5 / m) if mw < 1 else 3 * w
    jump = np.interp(zo, z, dKp) - np.interp(-zo, z, dKp)
    jump_corr = float(jump * np.exp(m * zo)) if mw < 1 else None   # thin wall: undo the e^{-m zo} decay outside
    track_err = np.max(np.abs(dK - rho / m ** 2)) / np.max(rho / m ** 2)
    return {'m_w': mw, 'kink_jump_dKprime_at_+-zo': float(jump), 'jump_times_exp(m zo) [thin walls only]': jump_corr,
            'expected_thin(-sigma)': -sigma, 'zo': float(zo),
            'max_rel_dev_from_local_tracking_rho/m^2': float(track_err)}
rowsD = [wall(mw) for mw in (0.01, 0.1, 1.0, 10.0, 30.0)]
res['D_linear_wall'] = {'equation': "dK'' - m_q^2 dK = T/3 = -rho_w(z) (static planar wall, T = -3 rho_w)",
                        'rows': rowsD,
                        'reading': 'm_q w << 1: K continuous with kink [dK\'] = -sigma (thin-wall junction); '
                                   'm_q w >> 1: q tracks the local quasi-static value, deviation O(1/(m_q w)^2)'}
for r_ in rowsD:
    print('D:', r_)

with open(OUT, 'w') as f:
    json.dump(res, f, indent=1, default=str)
