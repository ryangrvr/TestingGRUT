"""v_scenario / v1: independent symbolic checks of the scenario track's closed system (QV1).

Written by the verifier without using the scenario track's code. Every check is derived here from the action
    S = int sqrt(-g) [ K(q) R/2 - eps(q) ] + S_m,   q-equation eps' - K' R/2 = mu (constant)
via an FRW minisuperspace Lagrangian (dual form, q varied as an independent scalar with + mu q), then compared
with the statements in tracks/scenario/NOTES.md sections 1.1-1.6, 2.4, 5.3.

Checks
 A. minisuperspace -> FRW constraint, space equation, q-equation (generic K, eps)
 B. Bianchi: dC/dt = -4 H C with dissipative q-friction Gamma_q qdot provided matter gains Q = Gamma_q qdot^2
 C. KV forms: q0, mu0, m_q^2, rho_V(Y), R(Y), F(Y) polynomials, beta, Einstein-frame damping rate = Gamma
 D. R1 reduction: 3H^2(1+Y)+3HYdot = rho + rho_V with Y = O(1/m^2)  ->  delta rho_q formula
 E. rho_V = S^2/(12 K0 m^2) == KV (9/8) kappa^2 H^4 K'^2/eps''  (kappa = S/rho, 3K0H^2 = rho)
 F. oscillation form factor: Fourier transform of the tanh-step derivative -> (x/sinh x), x = pi m tau/2;
    sudden-step energy (4/3) Delta^2/(K0 m^2)
 G. P_q re-absorption with the KV quartic: (q0'/q0)^2 = (1+sqrt(1+8 chi c))/2, delta K/K, m^2 shifts
 H. fixed-mu static roots (mpmath, 120 digits) vs c(1 + (4/3)c/(K0 m^2) + c3 ...)
Output: ../results/v1_symbolic.json
"""
import json
import os

import mpmath as mp
import sympy as sp

out = {}
t = sp.symbols('t', real=True)
mu, Gq = sp.symbols('mu Gamma_q', real=True)
a = sp.Function('a')(t)
Nl = sp.Function('N')(t)
q = sp.Function('q')(t)
Kf = sp.Function('K')
ef = sp.Function('eps')

# ------------------------------------------------------------------ A. minisuperspace
# metric ds^2 = -N^2 dt^2 + a^2 dx^2 ; R = 6 ( addot/(a N^2) + adot^2/(a^2 N^2) - adot Ndot/(a N^3) )
R_ms = 6 * (a.diff(t, 2) / (a * Nl ** 2) + a.diff(t) ** 2 / (a ** 2 * Nl ** 2) - a.diff(t) * Nl.diff(t) / (a * Nl ** 3))
# matter: barotropic fluid L_m = -N a^3 rho(a)  (gives p = -d(rho a^3)/d(a^3))
rho_a = sp.Function('rho')
L = Nl * a ** 3 * (Kf(q) * R_ms / 2 - ef(q) + mu * q) - Nl * a ** 3 * rho_a(a)


def EL(Lag, f):
    return sp.diff(Lag, f) - sp.diff(sp.diff(Lag, f.diff(t)), t) + sp.diff(sp.diff(Lag, f.diff(t, 2)), t, 2)


E_N = sp.simplify(EL(L, Nl).subs(Nl, 1).doit())
E_a = sp.simplify(EL(L, a).subs(Nl, 1).doit())
E_q = sp.simplify(EL(L, q).subs(Nl, 1).doit())
H = a.diff(t) / a
K = Kf(q)
Kd = K.diff(t)
rhoV = ef(q) - mu * q
rho = rho_a(a)
p = -sp.diff(rho * a ** 3, a) / (3 * a ** 2)
C_mine = 3 * H ** 2 * K + 3 * H * Kd - rho - rhoV                      # Friedmann (00), scenario 1.4 form
space_mine = -(2 * H.diff(t) + 3 * H ** 2) * K - (p - rhoV + K.diff(t, 2) + 2 * H * Kd)   # equations-track 3.1 form
R_frw = 6 * (H.diff(t) + 2 * H ** 2)
qeq_mine = ef(q).diff(q) - Kf(q).diff(q) * R_frw / 2 - mu
r1 = sp.simplify(E_N - a ** 3 * C_mine)            # E_N = +a^3 C
r2 = sp.simplify(E_a - 3 * a ** 2 * space_mine)    # E_a = 3 a^2 * (space)
r2b = sp.simplify(E_a + 3 * a ** 2 * space_mine)
r3 = sp.simplify(E_q + a ** 3 * qeq_mine)
out['A_minisuperspace'] = dict(
    constraint_residual=str(r1), space_residual=str(r2) if r2 == 0 else str(r2b), q_eq_residual=str(r3),
    pass_=bool(r1 == 0 and (r2 == 0 or r2b == 0) and r3 == 0),
    note='E_N = +a^3 C, E_a = +-3a^2 (space eq), E_q = -a^3 (eps\' - K\'R/2 - mu): the scenario FRW equations follow '
         'from the action')

# ------------------------------------------------------------------ B. Bianchi with dissipation
# Treat q(t), H(t) as unknown functions; matter rho(t), p(t) generic; impose
#   (q') eps'(q) - K'(q) R/2 = mu - Gamma_q qdot          -> Hdot
#   (g_ii) -(2Hdot+3H^2)K = p - rhoV + Kddot + 2H Kdot    -> qddot (metric eq keeps constant mu)
#   (m) rhodot = -3H(rho+p) + Q
Hs, qs, vs, rs, ps, Qs = sp.symbols('H q v rho p Q', real=True)
Kq = sp.Function('K')(qs)
eq = sp.Function('eps')(qs)
K1, K2 = Kq.diff(qs), Kq.diff(qs, 2)
e1 = eq.diff(qs)
Hdot = ((e1 - mu + Gq * vs) * 2 / K1) / 6 - 2 * Hs ** 2
rV = eq - mu * qs
# Kddot = K1 qddot + K2 v^2 = -(2Hdot+3H^2)K - p + rV - 2 H K1 v
qddot = (-(2 * Hdot + 3 * Hs ** 2) * Kq - ps + rV - 2 * Hs * K1 * vs - K2 * vs ** 2) / K1
rhodot = -3 * Hs * (rs + ps) + Qs
C = 3 * Hs ** 2 * Kq + 3 * Hs * K1 * vs - rs - rV
dC = (sp.diff(C, Hs) * Hdot + sp.diff(C, qs) * vs + sp.diff(C, vs) * qddot + sp.diff(C, rs) * rhodot)
res_space = sp.simplify(sp.expand(dC + 3 * Hs * C))
# trace form (scenario 1.4 / equations 3.2): Kddot + 3H Kdot = (S + 4 rhoV - K R)/3, S = rho - 3p, R = 6(Hdot+2H^2)
Rq = 2 * (e1 - mu + Gq * vs) / K1
qddot_tr = ((rs - 3 * ps + 4 * rV - Kq * Rq) / 3 - 3 * Hs * K1 * vs - K2 * vs ** 2) / K1
dC_tr = (sp.diff(C, Hs) * Hdot + sp.diff(C, qs) * vs + sp.diff(C, vs) * qddot_tr + sp.diff(C, rs) * rhodot)
res = sp.simplify(sp.expand(dC_tr + 4 * Hs * C))
out['B_bianchi'] = dict(trace_form_dCdt_plus_4HC=str(res),
                        trace_form_with_Q_equal_Gamma_q_v2=str(sp.simplify(res.subs(Qs, Gq * vs ** 2))),
                        pass_=bool(sp.simplify(res.subs(Qs, Gq * vs ** 2)) == 0),
                        space_form_dCdt_plus_3HC_with_Q=str(sp.simplify(res_space.subs(Qs, Gq * vs ** 2))),
                        note='trace-form system (the one the scenario integrates): dC/dt = -4HC + (Gamma_q qdot^2 - Q); '
                             'with the space-space equation instead of the trace, dC/dt = -3HC + (same). Either way '
                             'the constraint surface is invariant iff Q = Gamma_q qdot^2')

# ------------------------------------------------------------------ C. KV forms
qq, q0, chi, M = sp.symbols('q q0 chi M', positive=True)
Y, c = sp.symbols('Y c', real=True)
epsKV = (1 / (2 * chi)) * (-qq ** 2 / q0 ** 2 + qq ** 4 / (3 * q0 ** 4))
KKV = qq / 3
mu0 = sp.diff(epsKV, qq).subs(qq, q0)
eq_resid = sp.simplify(epsKV.subs(qq, q0) - mu0 * q0)
K0 = KKV.subs(qq, q0)
m2 = sp.Rational(2, 3) * K0 * sp.diff(epsKV, qq, 2).subs(qq, q0) / sp.diff(KKV, qq) ** 2
subsnat = {q0: 3 * M ** 2, chi: 1 / M ** 4}
m2nat = sp.simplify(m2.subs(subsnat))
# rho_V(Y) and R(Y) in units K0 = 1 (M = 1), mu = mu0 - c/q0
nat1 = {q0: 3, chi: 1}
muc = mu0 - c / q0
rhoV_Y = sp.expand(sp.simplify((epsKV - muc * qq).subs(qq, q0 * (1 + Y)).subs(nat1)))
R_Y = sp.expand(sp.simplify((2 * (sp.diff(epsKV, qq) - muc) / sp.diff(KKV, qq)).subs(qq, q0 * (1 + Y)).subs(nat1)))
m2v = m2nat.subs(M, 1)
b3, b4, beta = sp.Rational(4, 3), sp.Rational(1, 3), 1
rhoV_claim = sp.Rational(3, 4) * m2v * Y ** 2 * (1 + b3 * Y + b4 * Y ** 2) + c * (1 + beta * Y)
R_claim = 3 * m2v * Y * (1 + sp.Rational(3, 2) * b3 * Y + 2 * b4 * Y ** 2) + 2 * beta * c
# trace equation Ydd + 3H Yd = [S + 4 rho_V - (1+Y) R]/3
S_ = sp.symbols('S')
force = sp.expand((S_ + 4 * rhoV_Y - (1 + Y) * R_Y) / 3)
F_claim = Y + sp.Rational(3, 2) * b3 * Y ** 2 + (b3 / 2 + 2 * b4) * Y ** 3 + b4 * Y ** 4
force_claim = S_ / 3 + (c / 3) * ((4 - 2 * beta) + 2 * beta * Y) - m2v * F_claim
# Einstein-frame damping with Gamma_q = (3/2) Gamma K0 K'^2/K^2: qdd gets -(2/3) K Gamma_q qd / K'^2
Gam = sp.symbols('Gamma', positive=True)
Kc, K1c, K0c = sp.symbols('K Kp K0', positive=True)
Gamma_q = sp.Rational(3, 2) * Gam * K0c * K1c ** 2 / Kc ** 2
damp = sp.simplify(sp.Rational(2, 3) * Kc * Gamma_q / K1c ** 2)
out['C_KV_forms'] = dict(
    mu0=str(sp.simplify(mu0)), eps_minus_mu0q_at_q0=str(eq_resid), m_q2_natural=str(m2nat),
    m_q_GeV=float(sp.sqrt(m2nat).subs(M, 2.435e18)),
    rhoV_Y=str(rhoV_Y), rhoV_matches=bool(sp.simplify(rhoV_Y - sp.expand(rhoV_claim)) == 0),
    R_Y=str(R_Y), R_matches=bool(sp.simplify(R_Y - sp.expand(R_claim)) == 0),
    trace_force_matches=bool(sp.simplify(force - sp.expand(force_claim)) == 0),
    damping_coefficient=str(damp), damping_is_Gamma_over_1pY=bool(sp.simplify(damp - Gam * K0c / Kc) == 0),
    note='energy of a damped oscillator qdd + gamma qd decays at rate gamma: Gamma at Y -> 0')

# ------------------------------------------------------------------ D. R1 reduction (formal small parameter e = 1/m^2)
e = sp.symbols('e', positive=True)
Sf = sp.Function('S')(t)
Hf = sp.Function('H')(t)
rhof = sp.Function('rho')(t)
y1, y2 = sp.Function('y1')(t), sp.Function('y2')(t)
Yt = e * y1 + e ** 2 * y2
Fy = Yt + 2 * Yt ** 2 + sp.Rational(4, 3) * Yt ** 3 + Yt ** 4 / 3
qeq = sp.expand(Yt.diff(t, 2) + 3 * Hf * Yt.diff(t) + Fy / e - Sf / 3)
o0 = qeq.coeff(e, 0)
sol_y1 = sp.solve(o0, y1)[0]
o1 = sp.expand(qeq.coeff(e, 1).subs(y1, sol_y1).doit())
sol_y2 = sp.solve(o1, y2)[0]
m2e = 1 / e
rhoV_e = sp.Rational(3, 4) * m2e * Yt ** 2 * (1 + sp.Rational(4, 3) * Yt + Yt ** 2 / 3)
Cexpr = 3 * Hf ** 2 * (1 + Yt) + 3 * Hf * Yt.diff(t) - rhof - rhoV_e
Cser = sp.expand(sp.series(Cexpr.subs({y1: sol_y1, y2: sol_y2}).doit(), e, 0, 2).removeO())
c0 = sp.simplify(Cser.coeff(e, 0))
c1 = sp.simplify(Cser.coeff(e, 1))
claim_c1 = -(Sf ** 2 / 12 - Hf ** 2 * Sf - Hf * Sf.diff(t))   # C = 3H^2 - rho - e*[S^2/12 - H^2 S - H Sdot]
out['D_R1'] = dict(y1=str(sol_y1), y2=str(sp.simplify(sol_y2)), C_order0=str(c0), C_order1=str(c1),
                   R1_matches=bool(sp.simplify(c1 - claim_c1) == 0),
                   note='3K0H^2 = rho + [S^2/(12K0) - H^2 S - H Sdot]/m^2 + O(m^-4); shape coefficients enter '
                        'only at O(m^-4)')

# ------------------------------------------------------------------ E. rho_V = S^2/(12 K0 m^2) vs KV kick
Sx, rx, K0x, Kp, e2 = sp.symbols('S rho K0 Kp e2', positive=True)
m2g = sp.Rational(2, 3) * K0x * e2 / Kp ** 2
H2 = rx / (3 * K0x)
kv = sp.Rational(9, 8) * (Sx / rx) ** 2 * H2 ** 2 * Kp ** 2 / e2
mine = Sx ** 2 / (12 * K0x * m2g)
out['E_kick_equivalence'] = dict(ratio=str(sp.simplify(kv / mine)), pass_=bool(sp.simplify(kv / mine) == 1))

# ------------------------------------------------------------------ F. oscillation form factor
# Linear oscillator Ydd + m^2 Y = S/3 (K0=1); S jumps by -4 Delta (DeltaV -> 0) with profile (1 - tanh(t/tau))/2.
# Amplitude after the step: |int dt e^{i m t} d/dt [S/3]| / m ; energy (3/4)[Yd^2 + m^2 Y^2] = (3/4) m^2 A^2.
mpm = mp.mp
mpm.dps = 30


def form_factor(m, tau):
    f = lambda s: mp.cos(m * s) * (1 / (2 * tau)) / mp.cosh(s / tau) ** 2
    return mp.quad(f, [-mp.inf, 0, mp.inf])


ff = {}
for xval in [0.1, 1.0, 2.0, 4.0, 8.0]:
    m_ = mp.mpf(1)
    tau = 2 * xval / (mp.pi * m_)
    num = form_factor(m_, tau)
    ff[str(xval)] = dict(numeric=float(num), x_over_sinhx=float(xval / mp.sinh(xval)))
Dl, mm = sp.symbols('Delta m', positive=True)
A_sudden = 4 * Dl / (3 * mm ** 2)
E_sudden = sp.simplify(sp.Rational(3, 4) * mm ** 2 * A_sudden ** 2)
out['F_form_factor'] = dict(table=ff, E_sudden=str(E_sudden),
                            pass_=all(abs(v['numeric'] / v['x_over_sinhx'] - 1) < 1e-12 for v in ff.values()),
                            note='energy = (4/3) Delta^2/(K0 m^2) (x/sinh x)^2, x = pi m tau/2 (tau: tanh width in t)')

# ------------------------------------------------------------------ G. P_q re-absorption (KV quartic)
u = sp.symbols('u', positive=True)
cc = sp.symbols('c', real=True)
P = sp.expand((qq * sp.diff(epsKV, qq) - epsKV).subs(qq, q0 * sp.sqrt(u)))
roots = sp.solve(sp.Eq(P, cc), u)
root = [r for r in roots if sp.limit(r, cc, 0) == 1][0]
dK = sp.series(sp.sqrt(root) - 1, cc, 0, 3).removeO()
e2new = sp.diff(epsKV, qq, 2).subs(qq, q0 * sp.sqrt(root))
m2new_over = sp.series((sp.sqrt(root) * e2new / sp.diff(epsKV, qq, 2).subs(qq, q0)), cc, 0, 2).removeO()
m2K_over = sp.series(e2new / sp.diff(epsKV, qq, 2).subs(qq, q0), cc, 0, 2).removeO()
out['G_Pq_reabsorption'] = dict(root_u=str(sp.simplify(root)), deltaK_over_K=str(sp.simplify(dK)),
                                m2_ratio=str(sp.simplify(m2new_over)), m2_over_K_ratio=str(sp.simplify(m2K_over)),
                                existence='1 + 8 chi c >= 0  <=> c >= -1/(8 chi) = -M^4/8')

# ------------------------------------------------------------------ H. fixed-mu static roots
mp.mp.dps = 400


def static_root(cM4):
    """KV quartic, M = 1 (chi = 1, q0 = 3, K = q/3), mu = mu0 fixed, eps -> eps + c. Unknown Y (q = 3(1+Y)), H^2.
    q-eq: eps'(q) - (1/3)(1/2)(12 H^2) = mu0 ; constraint 3 H^2 K = eps + c - mu0 q."""
    cM4 = mp.mpf(cM4)
    epsf = lambda qv: (-(qv ** 2) / 9 + qv ** 4 / 243) / 2
    deps = lambda qv: (-2 * qv / 9 + 4 * qv ** 3 / 243) / 2
    mu0v = deps(mp.mpf(3))
    rv = lambda Yv: epsf(3 * (1 + Yv)) + cM4 - mu0v * 3 * (1 + Yv)
    # H^2 = rho_V/(3K), K = 1+Y ; q-eq: deps - 2 H^2 = mu0
    g = lambda Yv: deps(3 * (1 + Yv)) - 2 * rv(Yv) / (3 * (1 + Yv)) - mu0v
    Yv = mp.findroot(g, cM4 * 2)       # first-order guess Y ~ 2 K1 c/(K0 e2) /q0 ...
    # cancellation-free cross-check with the polynomial identity of check C
    rp = lambda Z: Z * Z * (1 + 4 * Z / 3 + Z * Z / 3) / 2 + cM4
    gp = lambda Z: (Z + 2 * Z * Z + 2 * Z ** 3 / 3) / 3 - 2 * rp(Z) / (3 * (1 + Z))
    Zv = mp.findroot(gp, cM4 * 2)
    assert abs(rp(Zv) / rv(Yv) - 1) < mp.mpf(10) ** -60, (rp(Zv), rv(Yv))
    return Yv, rv(Yv)


Mred = mp.mpf('2.435e18')
M4 = Mred ** 4
cases = {'rho_Lambda': mp.mpf('2.52e-47'), 'Delta_EW': mp.mpf(125.25) ** 2 * mp.mpf(246.22) ** 2 / 8,
         'heavy_10TeV_abs': mp.mpf('1.5') * mp.mpf(1e4) ** 4 / (64 * mp.pi ** 2), 'toy_0.03': mp.mpf('0.03') * M4}
Hres = {}
for k, cv in cases.items():
    Yv, rv = static_root(cv / M4)
    ratio_m1 = rv / (cv / M4) - 1
    Hres[k] = dict(c_GeV4=mp.nstr(cv, 6), rhoV_over_c_minus_1=mp.nstr(ratio_m1, 8),
                   series_2c_over_M4=mp.nstr(2 * cv / M4, 8),
                   series_with_c3=mp.nstr(2 * cv / M4 - mp.mpf(32) / 3 * (cv / M4) ** 2, 8))


def static_root_mu_offset(cM4):
    """type ii (QV1-cal form): mu = mu0 - c/q0, eps unchanged. rho_V = eps - mu q = rho_V0(Y) + c(1+Y);
    q-eq static: eps'(q) - mu = (1/3)(1/2) 4 rho_V/(1+Y)  (K' = 1/3, R = 4 rho_V/K)."""
    cM4 = mp.mpf(cM4)
    rp = lambda Z: Z * Z * (1 + 4 * Z / 3 + Z * Z / 3) / 2 + cM4 * (1 + Z)
    gp = lambda Z: (Z + 2 * Z * Z + 2 * Z ** 3 / 3) / 3 + cM4 / 3 - 2 * rp(Z) / (3 * (1 + Z))
    Zv = mp.findroot(gp, cM4)
    return Zv, rp(Zv)


Hii = {}
for k, cv in cases.items():
    Zv, rv = static_root_mu_offset(cv / M4)
    Hii[k] = dict(rhoV_over_c_minus_1=mp.nstr(rv / (cv / M4) - 1, 8),
                  c_over_K0m2=mp.nstr(cv / M4 * mp.mpf(3) / 2, 8),
                  four_thirds_c_over_K0m2=mp.nstr(cv / M4 * 2, 8))
out['H_mu_offset_static_typeii'] = dict(cases=Hii, note='QV1-cal (mu offset): rho_V = c + c^2/(K0 m^2) + ..., '
                                                        'coefficient 1, not 4/3 (the 4/3 belongs to eps -> eps + c '
                                                        'at fixed mu)')
out['H_fixed_mu_static'] = dict(cases=Hres, note='(4/3) c/(K0 m^2) = 2 c/M^4 for natural KV scales; '
                                                 'c3 (equations-track formula, KV data) = -32/3 in M = 1 units')
# chi c values (P_q G renormalisation)
out['G_Pq_reabsorption']['chi_c'] = {k: mp.nstr(v / M4, 6) for k, v in cases.items()}
out['G_Pq_reabsorption']['relative_mu_offset_QV1cal'] = dict(
    c_over_M4=mp.nstr(cases['rho_Lambda'] / M4, 6),
    delta_mu_over_mu0=mp.nstr(3 * cases['rho_Lambda'] / M4, 6),
    note='|mu0| q0 = M^4/3 for the KV quartic, so delta mu/mu0 = 3 rho_L/M^4 (the "7.2e-121" is c/M^4)')

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, '..', 'results', 'v1_symbolic.json'), 'w') as f:
    json.dump(out, f, indent=1, default=str)
for k, v in out.items():
    print(k, {kk: vv for kk, vv in v.items() if kk in ('pass_', 'R1_matches', 'rhoV_matches', 'R_matches',
                                                        'trace_force_matches', 'damping_is_Gamma_over_1pY',
                                                        'm_q2_natural', 'm_q_GeV', 'deltaK_over_K', 'm2_ratio',
                                                        'm2_over_K_ratio', 'root_u', 'ratio', 'cases', 'chi_c',
                                                        'relative_mu_offset_QV1cal', 'with_Q_equal_Gamma_q_v2',
                                                        'trace_form_with_Q_equal_Gamma_q_v2', 'space_form_dCdt_plus_3HC_with_Q', 'constraint_residual', 'space_residual',
                                                        'q_eq_residual', 'table')})
