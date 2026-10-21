"""Independent minisuperspace derivation of the flat-FRW q-theory equations.

Action (brief): S = int d^4x sqrt(-g) [ K(q) R / 2 - eps(q) ] + S_m,
q = Ftilde / sqrt(-g), Ftilde metric independent.

Minisuperspace: ds^2 = -N(t)^2 dt^2 + a(t)^2 dx^2,  A_{123} = A(t)  =>  Ftilde = dA/dt
(normalisation of the 3-form drops out), sqrt(-g) = N a^3, q = A'/(N a^3).
Matter: perfect fluid with rho(a) (minisuperspace Lagrangian -N a^3 rho(a); gives
T_00 = rho, p = -rho - a rho'(a)/3, i.e. the fluid continuity equation).

Checks performed (results written to ../results/frw_minisuperspace.json):
  1. Euler-Lagrange in A  -> d/dt[ K'(q) R/2 - eps'(q) ] = 0, i.e. eps' - K' R/2 = mu (const).
  2. Euler-Lagrange in N  -> Friedmann constraint 3H^2 K + 3 H Kdot = rho + eps - mu q.
  3. Euler-Lagrange in a  -> space-space equation, compared with the covariant (g) equation
       -(2 Hdot + 3H^2) K = p - (eps - mu q) + Kddot + 2 H Kdot.
  4. Covariant (g) equation reduced to FRW by an explicit tensor computation (independent route).
  5. Constraint propagation of the first-order system used for numerics.
"""
import json, os
import sympy as sp

t = sp.symbols('t', real=True)
mu = sp.symbols('mu', real=True)
N = sp.Function('N')(t)
a = sp.Function('a')(t)
A = sp.Function('A')(t)
Kf = sp.Function('K')
ef = sp.Function('eps')
rhof = sp.Function('rho')

q = sp.diff(A, t) / (N * a**3)
sqrtg = N * a**3
# Ricci scalar of -N^2 dt^2 + a^2 dx^2
R = 6 * (sp.diff(a, t, 2) / (a * N**2) + sp.diff(a, t)**2 / (a**2 * N**2)
         - sp.diff(a, t) * sp.diff(N, t) / (a * N**3))
L = sqrtg * (Kf(q) * R / 2 - ef(q)) - sqrtg * rhof(a)


def EL(L, f):
    return sp.diff(L, f) - sp.diff(sp.diff(L, sp.diff(f, t)), t) + sp.diff(sp.diff(L, sp.diff(f, t, 2)), t, 2)


out = {}
# ---------- 1. A-equation ----------
pA = sp.diff(L, sp.diff(A, t))  # conjugate momentum; L has no A-dependence and no A'' dependence
# A'' enters only through q? no: q contains A' only. Check:
assert sp.diff(L, sp.diff(A, t, 2)) == 0
qs = sp.symbols('q_s')
pA_simpl = sp.simplify(pA.subs(sp.diff(A, t), qs * N * a**3))
expected_pA = (sp.Subs(sp.diff(Kf(qs), qs), qs, qs).doit() * R / 2 - sp.diff(ef(qs), qs))
check1 = sp.simplify(pA_simpl - expected_pA)
out['A_equation_momentum_minus_(K_prime_R/2-eps_prime)'] = str(check1)
# => dpA/dt = 0  => K'(q) R/2 - eps'(q) = const =: -mu  => eps' - K' R/2 = mu

# ---------- gauge N=1 helpers ----------
H = sp.Function('H')(t)
Q = sp.Function('Q')(t)  # q(t) in N=1 gauge


def to_gauge(expr):
    """Set N=1, express A' via q: A' = q a^3, and a' = H a."""
    e = expr
    e = e.subs(sp.diff(A, t, 3), sp.diff(Q * a**3, t, 2))
    e = e.subs(sp.diff(A, t, 2), sp.diff(Q * a**3, t))
    e = e.subs(sp.diff(A, t), Q * a**3)
    e = e.subs({sp.diff(N, t, 2): 0, sp.diff(N, t): 0}).subs(N, 1)
    e = e.doit()
    e = e.subs(sp.diff(a, t, 3), sp.diff(H * a, t, 2))
    e = e.subs(sp.diff(a, t, 2), sp.diff(H * a, t))
    e = e.subs(sp.diff(a, t), H * a)
    e = e.doit()
    e = e.subs(sp.diff(a, t, 2), sp.diff(H * a, t)).subs(sp.diff(a, t), H * a).doit()
    return sp.simplify(e)


# ---------- 2. N-equation (Friedmann constraint) ----------
EN = EL(L, N)
EN_g = to_gauge(EN)
# mu defined by the A-equation: eps'(Q) - K'(Q) R/2 = mu.  R in gauge:
Rg = 6 * (sp.diff(H, t) + 2 * H**2)
Kq = Kf(Q)
Kd = sp.diff(Kq, t)
rho_a = rhof(a)
# The constraint from EL_N should read  -(3H^2 K + 3 H Kdot) + rho + eps - Q eps' + Q K' R/2 = 0
# i.e. using the A-eq:  3H^2K + 3HKdot = rho + eps - mu Q
pred_N = -(3 * H**2 * Kq + 3 * H * Kd) + rho_a + ef(Q) - Q * sp.diff(ef(Q), Q) + Q * sp.diff(Kf(Q), Q) * Rg / 2
pred_N = sp.simplify(pred_N.doit())
ratio_N = sp.simplify(EN_g / (a**3))
out['N_equation_over_prediction'] = str(sp.simplify(ratio_N / pred_N))  # expect -1 (EL_N = -a^3 * constraint)

# ---------- 3. a-equation (space-space) ----------
Ea = EL(L, a)
Ea_g = to_gauge(Ea)
p_expr = -rho_a - a * sp.diff(rho_a, a) / 3  # fluid pressure implied by rho(a)
# covariant prediction (g)_ij with g_ij = a^2 delta_ij, divided by a^2:
#   -(2Hdot+3H^2)K = p - (eps - mu Q) + Kddot + 2 H Kdot,   mu = eps'(Q) - K'(Q) R/2
mu_onshell = sp.diff(ef(Q), Q) - sp.diff(Kf(Q), Q) * Rg / 2
pred_a = (-(2 * sp.diff(H, t) + 3 * H**2) * Kq) - (p_expr - (ef(Q) - mu_onshell * Q) + sp.diff(Kq, t, 2) + 2 * H * Kd)
pred_a = sp.simplify(pred_a.doit())
# EL_a = dL/da - ... ; relation: EL_a = -3 a^2 N * [space-space eq] up to sign; fit coefficient
coef = sp.simplify(Ea_g / pred_a)
out['a_equation_over_covariant_prediction'] = str(coef)

# ---------- 4. covariant tensor reduction (independent route) ----------
# Build FRW metric, compute G_mu_nu, Box K, nabla nabla K explicitly from Christoffels.
x, y, z = sp.symbols('x y z', real=True)
X = [t, x, y, z]
af = sp.Function('a')(t)
qf = sp.Function('Q')(t)
g = sp.diag(-1, af**2, af**2, af**2)
gi = g.inv()
Gam = [[[sum(gi[l, s] * (sp.diff(g[s, m], X[n]) + sp.diff(g[s, n], X[m]) - sp.diff(g[m, n], X[s])) for s in range(4)) / 2
         for n in range(4)] for m in range(4)] for l in range(4)]


def Riem(r, s, m, n):
    e = sp.diff(Gam[r][n][s], X[m]) - sp.diff(Gam[r][m][s], X[n])
    e += sum(Gam[r][m][l] * Gam[l][n][s] - Gam[r][n][l] * Gam[l][m][s] for l in range(4))
    return e


Ric = sp.Matrix(4, 4, lambda s, n: sp.simplify(sum(Riem(r, s, r, n) for r in range(4))))
Rs = sp.simplify(sum(gi[i, j] * Ric[i, j] for i in range(4) for j in range(4)))
Gmn = sp.simplify(Ric - Rs * g / 2)
Kc = Kf(qf)
dK = [sp.diff(Kc, X[m]) for m in range(4)]
DDK = sp.Matrix(4, 4, lambda m, n: sp.diff(Kc, X[m], X[n]) - sum(Gam[l][m][n] * dK[l] for l in range(4)))
BoxK = sp.simplify(sum(gi[m, n] * DDK[m, n] for m in range(4) for n in range(4)))
rho, p = sp.symbols('rho p')
T = sp.diag(rho, af**2 * p, af**2 * p, af**2 * p)
mu_c = sp.diff(ef(qf), qf) - sp.diff(Kf(qf), qf) * Rs / 2
E = Kc * Gmn - (T - (ef(qf) - mu_c * qf) * g - (g * BoxK - DDK))
Hs = sp.diff(af, t) / af
fr00 = sp.simplify(E[0, 0] - (3 * Hs**2 * Kc + 3 * Hs * sp.diff(Kc, t) - rho - (ef(qf) - mu_c * qf)))
frii = sp.simplify(E[1, 1] / af**2 - (-(2 * sp.diff(Hs, t) + 3 * Hs**2) * Kc - p + (ef(qf) - mu_c * qf) - sp.diff(Kc, t, 2) - 2 * Hs * sp.diff(Kc, t)))
out['covariant_FRW_00_residual'] = str(fr00)
out['covariant_FRW_ii_residual'] = str(frii)
out['covariant_FRW_R'] = str(sp.simplify(Rs - 6 * (sp.diff(af, t, 2) / af + Hs**2)))
out['covariant_FRW_BoxK'] = str(sp.simplify(BoxK + sp.diff(Kc, t, 2) + 3 * Hs * sp.diff(Kc, t)))

# FRW Bianchi (symbolic, generic K, eps, a, q):
#   div of [K G - (T_vac + (gBox-DD)K)] = ? ; equivalently the identity
#   d/dt of Friedmann constraint consistent with the ii equation + continuity iff mu' = 0.
# Direct covariant divergence of E_mn with mu -> mu(t) := eps' - K' R/2 (off shell):
#   nabla^m E_mn  =  -nabla^m T_mn - q d_n mu  (identity) -> check nu=0 component
Emat = Kc * Gmn + (ef(qf) - mu_c * qf) * g + (g * BoxK - DDK)  # "required T" for given geometry & q
# covariant divergence of symmetric tensor
def div(Tm):
    res = []
    for n in range(4):
        e = 0
        for m in range(4):
            for r in range(4):
                e += gi[m, r] * (sp.diff(Tm[r, n], X[m]) - sum(Gam[l][m][r] * Tm[l, n] + Gam[l][m][n] * Tm[r, l] for l in range(4)))
        res.append(sp.simplify(e))
    return res


dv = div(Emat)
ident0 = sp.simplify(dv[0] + qf * sp.diff(mu_c, t))
out['FRW_offshell_identity_div(E)_nu=0_plus_q_dmu'] = str(ident0)
out['FRW_offshell_identity_spatial'] = [str(sp.simplify(dv[i])) for i in (1, 2, 3)]

# ---------- 5. first-order system and constraint propagation ----------
qv, v, Hv, rd, rr, Vm, Vmdot, muv = sp.symbols('q v H rho_d rho_r V_m Vdot_m mu', real=True)
K = Kf(qv); K1 = sp.diff(K, qv); K2 = sp.diff(K, qv, 2)
e0 = ef(qv); e1 = sp.diff(e0, qv)
rhoV = e0 - muv * qv
rho_m = rd + rr + Vm
p_m = rr / 3 - Vm
Rq = 2 * (e1 - muv) / K1                     # from the q-equation
Kdot = K1 * v
# trace eq: 3 Box K = T - 4 rhoV + K R, Box K = -(Kddot + 3 H Kdot), T = -rho + 3p
Tm = -rho_m + 3 * p_m
Kddot = -3 * Hv * Kdot - (Tm - 4 * rhoV + K * Rq) / 3
vdot = (Kddot - K2 * v**2) / K1
Hdot = Rq / 6 - 2 * Hv**2
rddot = -3 * Hv * rd
rrdot = -4 * Hv * rr - Vmdot
C = 3 * Hv**2 * K + 3 * Hv * Kdot - rho_m - rhoV
state = [qv, v, Hv, rd, rr, Vm]
rates = [v, vdot, Hdot, rddot, rrdot, Vmdot]
Cdot = sum(sp.diff(C, s) * r for s, r in zip(state, rates))
lam = sp.simplify(Cdot / C)
out['constraint_propagation_Cdot_over_C'] = str(lam)
out['constraint_propagation_Cdot_plus_H_C'] = str(sp.simplify(Cdot + Hv * C))

for k_, v_ in out.items():
    print(k_, ':', v_)

os.makedirs(os.path.join(os.path.dirname(__file__), '..', 'results'), exist_ok=True)
with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'frw_minisuperspace.json'), 'w') as f:
    json.dump(out, f, indent=1)
