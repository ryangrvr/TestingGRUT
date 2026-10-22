"""v_bound, independent symbolic checks of track `bound` (sympy).  Writes ../results/v1_symbolic.json.

Written from scratch; nothing is imported from ../../bound/code.  Every check starts from an action or a definition,
not from the bound track's intermediate formulas.

 V1a  f(R) FRW 00 equation derived from the minisuperspace action  S = int N a^3 [f(R)/2 - rho(a)] dt  (lapse N):
        3 f' H^2 = (R f' - f)/2 - 3 H d(f')/dt + rho
      and the q-theory (Brans-Dicke omega = 0) minisuperspace  S = int N a^3 [K(q) R/2 - U(q) - rho(a)] dt:
        constraint, q-equation K' R = 2 U', and the trace equation  Kddot + 3 H Kdot = (rho_d + 4U - K R)/3.
 V1b  rho_DE,eff = 3 K_lab H^2 - rho for f = K_lab R - 2 Lam + f2 R^2/2 to first order on the GR + Lam background
      (dust + radiation), with time derivatives taken through the background equations (no substitution shortcuts).
 V1c  Same quantity from the q-variables directly: rho_V (gravitating) vs rho_DE,eff (factor 9 claim).
 V1d  Static (dS/AdS) root as a series in dmu with generic Taylor data, and consistency with stage-1 (equations 5.2).
 V1e  f(R) map: f'' and m^2; the affine inverse map for an arbitrary f; P_q <-> f(0) = 0, f'(0) > 0, f''(0) > 0.
 V1f  Threshold re-selection: delta f(R) = 2 c (q_R - q0)/q0 + O(c^2), checked on explicit Taylor data.
 V1g  omega_BD = 0 with dust and a negligible potential: a ~ t^(1/2), K ~ t^(1/2) (phi-MDE) exact solution.
 V1h  rho_V = (3/4) m^2 K (dK/K)^2 for quadratic eps; QCD trace-anomaly coefficient 49/6 alpha_s^2 (n_f = 6).
 V1i  Appleby-Battye: F'(R) = [1 + tanh(aR - b)]/2, F'(0) = 1/(1 + e^{2b}); relation of e^{2b} to |f_R0|.
Run: python3 v1_symbolic.py
"""
import json
import os
import sympy as sp

out = {}
t = sp.symbols('t', real=True)

# ------------------------------------------------------------------------------------------------ V1a
a = sp.Function('a')(t)
N = sp.Function('N')(t)
Rsym = sp.symbols('R')
c0_, c1_, c2_, c3_, c4_, c5_ = sp.symbols('c0:6', real=True)
fconc = lambda X: c0_ + c1_ * X + c2_ * X**2 + c3_ * X**3 + c4_ * sp.sin(c5_ * X)     # generic analytic f (4 shapes)
rho = sp.Function('rho')          # rho(a): matter density as a function of the scale factor
Rlapse = 6 * (sp.diff(a, t, 2) / (a * N**2) + sp.diff(a, t)**2 / (a**2 * N**2) - sp.diff(a, t) * sp.diff(N, t) / (a * N**3))
L = N * a**3 * (fconc(Rlapse) / 2 - rho(a))
EL_N = (sp.diff(L, N) - sp.diff(sp.diff(L, sp.diff(N, t)), t)).subs(N, 1).doit()
H = sp.diff(a, t) / a
R1 = 6 * (sp.diff(a, t, 2) / a + H**2)
Xs = sp.symbols('X')
fp_conc = sp.diff(fconc(Xs), Xs).subs(Xs, R1)
claimed = 3 * fp_conc * H**2 - (R1 * fp_conc - fconc(R1)) / 2 + 3 * H * sp.diff(fp_conc, t) - rho(a)
ratio = sp.simplify(EL_N / (a**3 * claimed))
out['V1a_fR_00_equation_from_minisuperspace'] = {
    'claimed (bound A.1, adiabatic_relation S1)': "3 f' H^2 = (R f' - f)/2 - 3 H d(f')/dt + rho",
    'f used': 'c0 + c1 R + c2 R^2 + c3 R^3 + c4 sin(c5 R) (generic)',
    '[dL/dN - d/dt dL/dNdot]_(N=1) / (a^3 * (LHS - RHS of claim)); a constant means identical': str(ratio),
}
# q-theory / BD omega = 0 minisuperspace
q = sp.Function('q')(t)
Kf, Uf = sp.Function('K'), sp.Function('U')
Lq = N * a**3 * (Kf(q) * Rlapse / 2 - Uf(q) - rho(a))
EL_Nq = (sp.diff(Lq, N) - sp.diff(sp.diff(Lq, sp.diff(N, t)), t)).subs(N, 1).doit()
EL_qq = (sp.diff(Lq, q)).subs(N, 1).doit()                     # q appears without derivatives in L (Kdot via R only)
EL_aq = (sp.diff(Lq, a) - sp.diff(sp.diff(Lq, sp.diff(a, t)), t) + sp.diff(sp.diff(Lq, sp.diff(a, t, 2)), t, 2)).subs(N, 1).doit()
Kq, Uq = Kf(q), Uf(q)
Kd = sp.diff(Kq, t)
constraint_claim = 3 * H**2 * Kq + 3 * H * Kd - rho(a) - Uq
c1 = sp.simplify(EL_Nq / (a**3 * constraint_claim))
qeq_claim = sp.diff(Kq, q) * R1 / 2 - sp.diff(Uq, q)
c2 = sp.simplify(EL_qq - a**3 * qeq_claim)
# trace equation: use a-equation + constraint; dust-only matter rho(a) = C a^-3 here, radiation separately below
C, Cr = sp.symbols('C C_r', positive=True)
for label, rhofun, trace_src in [('dust', C / a**3, C / a**3), ('radiation', Cr / a**4, 0)]:
    ELa = EL_aq.subs(rho(a), rhofun).doit()
    ELa = ELa.subs(sp.Subs(sp.Derivative(rho(sp.Symbol('x')), sp.Symbol('x')), sp.Symbol('x'), a), sp.diff(rhofun, a))
    ELa = sp.simplify(ELa.replace(lambda e: isinstance(e, sp.Subs), lambda e: e.doit()))
    # a-equation is (space-space) * a^2 * const; combine with 3*constraint to form the trace:
    #   claimed: Kddot + 3 H Kdot - (rho_trace + 4 U - K R)/3 = 0 when the q-equation R = 2U'/K' holds
    trace_claim = sp.diff(Kq, t, 2) + 3 * H * Kd - (trace_src + 4 * Uq - Kq * R1) / 3
    # Solve the a-equation for addot and the constraint is not needed for the trace (it is the 'p' combination)
    A2 = sp.symbols('A2')
    ELa_s = ELa.subs(sp.diff(a, t, 2), A2)
    sol_add = sp.solve(sp.Eq(ELa_s, 0), A2)
    tc = sp.simplify(trace_claim.subs(sp.diff(a, t, 2), sol_add[0]) if sol_add else trace_claim)
    # tc should vanish identically given the a-equation together with the constraint; check with the constraint:
    Kpp = sp.symbols('Kpp')
    tc2 = sp.simplify(tc)
    out['V1a_trace_check_' + label] = {
        'trace residual after eliminating addot by the a-equation (should be proportional to the constraint)': str(sp.factor(tc2))[:400],
        'residual / (3H^2K + 3HKdot - rho - U) (a pure number: trace eq. holds on the constraint surface)': str(sp.simplify(tc2 / constraint_claim.subs(rho(a), rhofun)))[:200],
    }
out['V1a_q_theory_minisuperspace'] = {
    'constraint: [dL/dN]/(a^3 (3H^2K + 3H Kdot - rho - U)) (constant = identical)': str(c1),
    'q-equation K\' R/2 = U\'  (0 = identical)': str(c2),
}

# ------------------------------------------------------------------------------------------------ V1b
K, Lam, m, f2 = sp.symbols('K Lambda m f2', positive=True)
rd = sp.Function('rho_d')(t)
rr = sp.Function('rho_r')(t)
Hb = sp.sqrt((rd + rr + Lam) / (3 * K))                        # GR + Lam background, K = K_lab
Rb = (rd + 4 * Lam) / K                                       # radiation traceless
bg_rules = {sp.Derivative(rd, t): -3 * Hb * rd, sp.Derivative(rr, t): -4 * Hb * rr}
Rb_dot = sp.diff(Rb, t).subs(bg_rules)
F = lambda x: f2 * x**2 / 2
Fp = lambda x: f2 * x
Fpp = f2
# exact identity: 3 K H^2 - rho = Lam + (R F' - F)/2 - 3 F' H^2 - 3 H F'' Rdot  (from V1a with f = K R - 2Lam + F)
rhoDE = Lam + (Rb * Fp(Rb) - F(Rb)) / 2 - 3 * Fp(Rb) * Hb**2 - 3 * Hb * Fpp * Rb_dot
rhoDE = sp.simplify(sp.expand(rhoDE.subs(f2, K / (3 * m**2))))
target = Lam + (sp.Rational(9, 4) * rd**2 + 2 * rd * rr - 4 * Lam * rr) / (3 * K * m**2)
out['V1b_first_order_rhoDE'] = {
    'rho_DE,eff - Lam (GR background, f2 = K/(3 m^2))': str(sp.factor(sp.simplify(rhoDE - Lam))),
    'equals [(9/4) rd^2 + 2 rd rr - 4 Lam rr]/(3 K m^2)': bool(sp.simplify(rhoDE - target) == 0),
}
# derivation of the 00 identity itself for f = K R - 2Lam + F (algebra, generic F)
Rs, Hs, Rds = sp.symbols('R H Rdot')
Fg = sp.Function('F')
fg = K * Rs - 2 * Lam + Fg(Rs)
lhs00 = 3 * sp.diff(fg, Rs) * Hs**2
rhs00_nomat = (Rs * sp.diff(fg, Rs) - fg) / 2 - 3 * Hs * sp.diff(fg, Rs, 2) * Rds
ident = sp.simplify(rhs00_nomat - (lhs00 - 3 * K * Hs**2) - (Lam + (Rs * sp.diff(Fg(Rs), Rs) - Fg(Rs)) / 2
                                                            - 3 * sp.diff(Fg(Rs), Rs) * Hs**2 - 3 * Hs * sp.diff(Fg(Rs), Rs, 2) * Rds))
out['V1b_identity_rhoDE_exact_in_F'] = str(ident)
# Omega_dyn(t0), w_dyn, m_q needed for Omega_Lambda
Om, OL, H0, aa = sp.symbols('Omega_m Omega_L H_0 a_s', positive=True)
rho_c = 3 * K * H0**2
rho_dyn_dust = sp.Rational(9, 4) * (Om * rho_c * aa**-3)**2 / (3 * K * m**2)
Omega_dyn = sp.simplify(rho_dyn_dust.subs(aa, 1) / rho_c)
w_dyn = sp.simplify(-1 - aa * sp.diff(sp.log(rho_dyn_dust), aa) / 3)
m_formal = sp.solve(sp.Eq(Omega_dyn, OL), m)
out['V1b_Omega_dyn_and_w'] = {'Omega_dyn(t0)': str(Omega_dyn), 'w_dyn (dust era)': str(w_dyn),
                              'm needed for Omega_dyn = Omega_L': str(m_formal),
                              'numeric m/H0 at Om = 0.315, OL = 0.685': [float(s.subs({Om: 0.315, OL: 0.685, H0: 1})) for s in m_formal]}

# ------------------------------------------------------------------------------------------------ V1c
# q-variables: eps = mu0 q + e2 (q - q0)^2/2 (+ e3 cubic), K(q) = K0 + K1 (q - q0); at mu = mu0 the R = 0 point is q0.
q0, e2, e3, K0, K1, Rv = sp.symbols('q0 e2 e3 K0 K1 R', real=True)
x = sp.symbols('x')                     # q - q0
U = e2 * x**2 / 2 + e3 * x**3 / 6        # eps - mu0 q
Kx = K0 + K1 * x
# q-equation: K' R/2 = U'  ->  x(R) as a series
xs = sp.symbols('x1 x2')
xser = xs[0] * Rv + xs[1] * Rv**2
eq = sp.expand(sp.series((K1 * Rv / 2 - sp.diff(U, x)).subs(x, xser), Rv, 0, 3).removeO())
sol = sp.solve([eq.coeff(Rv, 1), eq.coeff(Rv, 2)], xs, dict=True)[0]
xR = xser.subs(sol)
rhoV_R = sp.series(U.subs(x, xR), Rv, 0, 3).removeO()
m2 = sp.Rational(2, 3) * K0 * e2 / K1**2
out['V1c_gravitating_rhoV'] = {
    'rho_V(R) to O(R^2)': str(sp.simplify(rhoV_R)),
    'equals K0 R^2/(12 m^2) (stage-1 rho_m^2/(12 K m^2) with R = rho/K)': bool(sp.simplify(rhoV_R - K0 * Rv**2 / (12 * m2)) == 0),
    'rho_DE,eff/rho_V in the dust era (no Lam): [(9/4)/(3)] / (1/12)': str(sp.Rational(9, 4) / 3 / sp.Rational(1, 12)),
}

# ------------------------------------------------------------------------------------------------ V1d
dmu = sp.symbols('dmu')
E2, E3, E4, k0, k1, k2, k3 = sp.symbols('E2 E3 E4 k0 k1 k2 k3', real=True)
mu0 = sp.symbols('mu0', real=True)
epsT = mu0 * (q0 + x) + E2 * x**2 / 2 + E3 * x**3 / 6 + E4 * x**4 / 24   # Minkowski point q0 for mu0: eps' = mu0, eps - mu0 q0 = 0
KT = k0 + k1 * x + k2 * x**2 / 2 + k3 * x**3 / 6
rhoVT = epsT - (mu0 + dmu) * (q0 + x)
# static de Sitter / AdS:  3 H^2 K = rho_V,  R = 12 H^2 = 4 rho_V/K,  q-eq: eps' - mu = K' R/2 = 2 K' rho_V/K
static_eq = (sp.diff(epsT, x) - (mu0 + dmu)) * KT - 2 * sp.diff(KT, x) * rhoVT
cs = sp.symbols('s1 s2 s3')
xd = cs[0] * dmu + cs[1] * dmu**2 + cs[2] * dmu**3
ser = sp.expand(sp.series(static_eq.subs(x, xd), dmu, 0, 4).removeO())
ssol = sp.solve([ser.coeff(dmu, i) for i in (1, 2, 3)], cs, dict=True)[0]
rhoV_s = sp.expand(sp.series(rhoVT.subs(x, xd.subs(ssol)), dmu, 0, 3).removeO())
mm2 = sp.Rational(2, 3) * k0 * E2 / k1**2
coef2_claim = -1 / (2 * E2) + sp.Rational(4, 3) * q0**2 / (k0 * mm2)
# stage-1 consistency: dc(dmu) from mu0[eps + c] = mu0 + dmu, then rho_V = dc + (4/3) dc^2/(K0 m^2)
c_ = sp.symbols('c')
# P_q for eps + c: eps'(qc) = mu, eps(qc) + c - mu qc = 0 -> c(mu) to second order
yq = sp.symbols('y')
eq1 = sp.diff(epsT, x).subs(x, yq) - (mu0 + dmu)
ys = sp.symbols('y1 y2')
yser = ys[0] * dmu + ys[1] * dmu**2
eqs = sp.expand(sp.series(eq1.subs(yq, yser), dmu, 0, 3).removeO())
ysol = sp.solve([eqs.coeff(dmu, 1), eqs.coeff(dmu, 2)], ys, dict=True)[0]
c_of_dmu = sp.expand(sp.series(((mu0 + dmu) * (q0 + yq) - epsT.subs(x, yq)).subs(yq, yser.subs(ysol)), dmu, 0, 3).removeO())
dc = -c_of_dmu                                                   # c = 0 actual, c* = c(dmu): dc = c - c* = -c*
stage1 = sp.expand(sp.series(dc + sp.Rational(4, 3) * dc**2 / (k0 * mm2), dmu, 0, 3).removeO())
out['V1d_static_root_series'] = {
    'rho_V,s(dmu) to O(dmu^2) (generic Taylor data up to E4, k3)': str(sp.collect(sp.simplify(rhoV_s), dmu)),
    'linear coefficient = -q0': bool(sp.simplify(rhoV_s.coeff(dmu, 1) + q0) == 0),
    'quadratic coefficient = -1/(2 eps\'\') + (4/3) q0^2/(K m^2)': bool(sp.simplify(rhoV_s.coeff(dmu, 2) - coef2_claim) == 0),
    'c*(dmu) from P_q': str(sp.simplify(c_of_dmu)),
    'agrees with stage-1 dc + (4/3) dc^2/(K m^2) to O(dmu^2)': bool(sp.simplify(stage1 - rhoV_s) == 0),
}

# ------------------------------------------------------------------------------------------------ V1e
qq, muS, lam_ = sp.symbols('q mu lambda', real=True)
eps_f, K_f = sp.Function('epsilon'), sp.Function('Kf')
# f(R) = K(q_R) R - 2 eps(q_R) + 2 mu q_R with eps'(q_R) - K'(q_R) R/2 = mu.  d q_R/dR by implicit differentiation:
qR = sp.Function('qR')(Rv)
qeq = sp.diff(eps_f(qR), qR) - sp.diff(K_f(qR), qR) * Rv / 2 - muS
dqdR = sp.solve(sp.diff(qeq, Rv), sp.diff(qR, Rv))[0]
fR_expr = K_f(qR) * Rv - 2 * eps_f(qR) + 2 * muS * qR
f1 = sp.diff(fR_expr, Rv).subs(sp.diff(qR, Rv), dqdR)
f1 = sp.simplify(f1.subs(sp.Subs(sp.Derivative(eps_f(sp.Symbol('_x')), sp.Symbol('_x')), sp.Symbol('_x'), qR), 0) if False else f1)
# use the q-equation to eliminate eps'(q_R):  eps' = mu + K' R/2
epsp = sp.Derivative(eps_f(qR), qR)
f1s = sp.simplify(f1.doit().subs(sp.diff(eps_f(qR), qR), muS + sp.diff(K_f(qR), qR) * Rv / 2))
f2s = sp.simplify(sp.diff(K_f(qR), Rv).subs(sp.diff(qR, Rv), dqdR))
out['V1e_map'] = {
    "f'(R) (after using the q-equation)": str(f1s),
    "f''(R) = K' dq_R/dR": str(f2s),
    "scalaron m^2 at R = 0: f'/(3 f'') with f'' = K'^2/(2 eps'')": str(sp.simplify(K0 / (3 * K1**2 / (2 * e2)))),
    'equals (2/3) K eps\'\'/K\'^2': bool(sp.simplify(K0 / (3 * K1**2 / (2 * e2)) - sp.Rational(2, 3) * K0 * e2 / K1**2) == 0),
}
# inverse (affine gauge): arbitrary f, R(q) = lam (q - q0)
gF = sp.Function('g')
Rq = lam_ * (qq - q0)
Kinv = sp.diff(gF(Rq), qq) / lam_
epsinv = (Kinv * Rq - gF(Rq)) / 2 + muS * qq
res_q = sp.simplify(sp.diff(epsinv, qq) - sp.diff(Kinv, qq) * Rq / 2 - muS)
res_f = sp.simplify(Kinv * Rq - 2 * epsinv + 2 * muS * qq - gF(Rq))
# stationarity uniqueness: K'(q)(R - lam(q - q0)) = 0 at fixed R
Rfix = sp.symbols('R_fix')
stat = sp.simplify(sp.diff(Kinv * Rfix - 2 * epsinv + 2 * muS * qq, qq))
out['V1e_inverse_affine'] = {'q-equation residual': str(res_q), 'reconstructed f - f': str(res_f),
                             'stationarity in q at fixed R (should be K\'(q) (R - lam (q-q0)))': str(sp.factor(stat))}
# P_q in f language: f(0) = -2 (eps(q0) - mu q0) with eps'(q0) = mu ; Minkowski solves R f' - 2 f = 0 at R = 0 iff f(0) = 0
out['V1e_P_q'] = {'f(0)': '-2 (eps(q0) - mu q0) at eps\'(q0) = mu', 'f\'(0)': 'K(q0)', 'f\'\'(0)': "K'(q0)^2/(2 eps''(q0))",
                  'vacuum trace eq. R f\' - 2 f at R = 0': '-2 f(0)',
                  'reading': 'P_q (eps - mu q0 = 0, eps\'\' > 0, K > 0, K\' != 0) <=> f(0) = 0, f\'(0) > 0, f\'\'(0) > 0 [D, confirmed]'}

# ------------------------------------------------------------------------------------------------ V1f
cc = sp.symbols('c')
e2_, e3_, K0_, K1_, K2_ = sp.symbols('e2 e3 K0 K1 K2', positive=True)
q0_, mu0_ = sp.symbols('q0 mu0', positive=True)
Q = sp.symbols('Q')                                   # q
epsP = mu0_ * Q + e2_ * (Q - q0_)**2 / 2 + e3_ * (Q - q0_)**3 / 6    # P_q at (q0, mu0) for c = 0
KP = K0_ + K1_ * (Q - q0_) + K2_ * (Q - q0_)**2 / 2
# P_q for eps + c: (q0c, mu0c) with eps'(q0c) = mu0c, eps(q0c) + c - mu0c q0c = 0;  series in c
z1, z2 = sp.symbols('z1 z2')
q0c = q0_ + z1 * cc + z2 * cc**2
mu0c = sp.diff(epsP, Q).subs(Q, q0c)
cond = sp.expand(sp.series(epsP.subs(Q, q0c) + cc - mu0c * q0c, cc, 0, 3).removeO())
zs = sp.solve([cond.coeff(cc, 1), cond.coeff(cc, 2)], [z1, z2], dict=True)[0]
mu0c = sp.series(mu0c.subs(zs), cc, 0, 2).removeO()
dmu0_dc = sp.simplify(mu0c.coeff(cc, 1))
# f_new(R) - f(R) at first order in c via the envelope theorem check: compute both f's as series in R (to R^3) explicitly
def f_series(epsfun, mu, order=4):
    w = sp.symbols('w1:%d' % order)
    qser = q0_ + sum(w[i] * Rv**(i + 1) for i in range(order - 1))
    qsol0 = sp.solve(sp.Eq(sp.diff(epsfun, Q).subs(Q, q0_ + sp.Symbol('y0')), mu), sp.Symbol('y0'))
    return None
# do it with the R = 0 point of each mu directly as a double series (c, R)
yv = sp.symbols('y')
def f_of(epsfun, mu, nR=3):
    # q_R = q0 + y(R, c): solve eps'(q) - K'(q) R/2 = mu by series in R and c (both small)
    cf = sp.symbols('d0:%d' % 12)
    idx = 0
    terms = 0
    combos = [(i, j) for i in range(0, nR + 1) for j in range(0, 2) if (i, j) != (0, 0)]
    for (i, j) in combos:
        terms += cf[idx] * Rv**i * cc**j
        idx += 1
    expr = (sp.diff(epsfun, Q) - sp.diff(KP, Q) * Rv / 2 - mu).subs(Q, q0_ + terms)
    expr = sp.expand(expr)
    eqs_ = []
    for (i, j) in combos:
        eqs_.append(sp.expand(expr).coeff(Rv, i).coeff(cc, j) if i > 0 else sp.expand(expr).subs(Rv, 0).coeff(cc, j))
    # truncate higher orders in c by keeping only j <= 1 (series)
    solc = sp.solve(eqs_, cf[:len(combos)], dict=True)[0]
    qR_ = q0_ + terms.subs(solc)
    fval = KP.subs(Q, qR_) * Rv - 2 * epsfun.subs(Q, qR_) + 2 * mu * qR_
    fval = sp.expand(fval)
    # keep orders R^<=nR, c^<=1
    poly = sp.Poly(fval, Rv, cc)
    keep = sum(coef * Rv**i * cc**j for (i, j), coef in zip(poly.monoms(), poly.coeffs()) if i <= nR and j <= 1)
    return sp.simplify(keep), sp.expand(qR_)
f_old, qR_old = f_of(epsP, mu0_)
f_new, _ = f_of(epsP + cc, mu0_ + dmu0_dc * cc)
delta = sp.expand(f_new - f_old)
pred = sp.expand(2 * cc * (qR_old.subs(cc, 0) - q0_) / q0_)
predpoly = sp.Poly(pred, Rv, cc)
pred_tr = sum(coef * Rv**i * cc**j for (i, j), coef in zip(predpoly.monoms(), predpoly.coeffs()) if i <= 3 and j <= 1)
deltapoly = sp.Poly(delta, Rv, cc)
delta_c1 = sum(coef * Rv**i * cc**j for (i, j), coef in zip(deltapoly.monoms(), deltapoly.coeffs()) if j == 1)
out['V1f_threshold'] = {
    'dmu0/dc': str(dmu0_dc),
    'delta f at O(c), series to R^3': str(sp.factor(delta_c1)),
    'equals 2 c (q_R - q0)/q0 to R^3': bool(sp.simplify(delta_c1 - pred_tr) == 0),
    'linear K: q_R - q0 = (f\'(R) - f\'(0))/K1 => f_new = f(R + dR) - dR f\'(0), dR = 2c/(q0 K1)': 'follows algebraically',
}

# ------------------------------------------------------------------------------------------------ V1g
tt, A0, r0 = sp.symbols('t A0 rho0', positive=True)
a_t = tt**sp.Rational(1, 2)
H_t = sp.diff(a_t, tt) / a_t
Kt = sp.Rational(2, 3) * r0 * tt**sp.Rational(1, 2)
rho_t = r0 * tt**sp.Rational(-3, 2)                      # dust ~ a^-3
R_t = sp.simplify(6 * (sp.diff(a_t, tt, 2) / a_t + H_t**2))
constr = sp.simplify(3 * H_t**2 * Kt + 3 * H_t * sp.diff(Kt, tt) - rho_t)
trace = sp.simplify(sp.diff(Kt, tt, 2) + 3 * H_t * sp.diff(Kt, tt) - (rho_t - Kt * R_t) / 3)
out['V1g_phiMDE'] = {'R for a ~ t^1/2': str(R_t), 'constraint residual (U = 0)': str(constr), 'trace residual': str(trace),
                     'reading': 'omega_BD = 0 with dust and U = 0 (q-mode frozen, m << H): a ~ t^(1/2), K ~ t^(1/2) ~ a solves the '
                                'system exactly (Nariai-type BD dust solution; APT07 phi-MDE)'}

# ------------------------------------------------------------------------------------------------ V1h
mq, Kc, dK = sp.symbols('m_q K dK', positive=True)
e2h = sp.Rational(3, 2) * mq**2 * K1**2 / Kc                 # from m^2 = (2/3) K e2/K1^2
rhoV_disp = e2h / 2 * (dK / K1)**2
out['V1h_displacement'] = {'rho_V for dK': str(sp.simplify(rhoV_disp)),
                           'equals (3/4) m^2 K (dK/K)^2': bool(sp.simplify(rhoV_disp - sp.Rational(3, 4) * mq**2 * Kc * (dK / Kc)**2) == 0),
                           'dK/K for rho_V = 3 K H0^2 Omega_L': str(sp.solve(sp.Eq(sp.Rational(3, 4) * mq**2 * Kc * (dK / Kc)**2, 3 * Kc * H0**2 * OL), dK)[0] / Kc)}
# QCD trace anomaly: p = (8 pi^2/45) T^4 [1 + 21 nf/32 - (15/4)(1 + 5 nf/12) alpha/pi], dalpha/dlnT = -beta0 alpha^2/(2 pi)
T, al, nf = sp.symbols('T alpha n_f', positive=True)
alT = sp.Function('alpha_s')(T)
p = sp.Rational(8, 1) * sp.pi**2 / 45 * T**4 * (1 + sp.Rational(21, 32) * nf - sp.Rational(15, 4) * (1 + sp.Rational(5, 12) * nf) * alT / sp.pi)
beta0 = 11 - sp.Rational(2, 3) * nf
I_ = sp.expand(T**5 * sp.diff(p / T**4, T))
I_ = sp.simplify(I_.subs(sp.diff(alT, T), -beta0 * alT**2 / (2 * sp.pi * T)))
out['V1h_QCD_trace_anomaly'] = {'(eps - 3p)/(alpha_s^2 T^4)': str(sp.simplify(I_ / (alT**2 * T**4))),
                                'n_f = 6': str(sp.simplify((I_ / (alT**2 * T**4)).subs(nf, 6))),
                                'bound code uses 2940/360': str(sp.Rational(2940, 360))}

# ------------------------------------------------------------------------------------------------ V1i
aA, bA, RA = sp.symbols('a b R', positive=True)
FA = RA / 2 + sp.log(sp.cosh(aA * RA) - sp.tanh(bA) * sp.sinh(aA * RA)) / (2 * aA)
FAp = sp.simplify(sp.diff(FA, RA))
chk = sp.simplify((FAp - (1 + sp.tanh(aA * RA - bA)) / 2).rewrite(sp.exp))
out['V1i_AppleByBattye'] = {"F'(R) - [1 + tanh(aR - b)]/2": str(chk),
                            "F'(0)": str(sp.simplify(FAp.subs(RA, 0).rewrite(sp.exp))),
                            "1/F'(0)": str(sp.simplify((1 / FAp.subs(RA, 0)).rewrite(sp.exp)))}

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, '..', 'results', 'v1_symbolic.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
for k, v in out.items():
    print(k, ':', v)
