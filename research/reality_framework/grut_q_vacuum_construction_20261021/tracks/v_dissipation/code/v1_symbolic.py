"""v_dissipation, script 1: independent symbolic checks of the derivations in ../../dissipation/NOTES.md.

Nothing is imported from ../../dissipation/code.  Formulations differ from the track's s1_symbolic.py:
  A  the FRW components of the BRIEF (g) operator are built here from an explicit Christoffel computation of
     G_mn, Box K and nabla_m nabla_n K on ds^2 = -dt^2 + a^2 dx^2 (the track used pre-reduced FRW formulas);
  B  KV16a unit map: solve for the energy unit E that makes the printed (6.1) equal to the dimensionful (5.6);
  C  KV16a App. B solution substituted into the three first-order ODEs, and the leading-order asymptote;
  D  per-e-fold drain law, radiation closed form, trace-drain integrals;
  E  KSV16 (17a-c), (18)-(21) re-derived from the BRIEF equations with K = f/3, eps generic; constraint propagation;
  F  K22 endpoint violation and the ad hoc asymptote coefficients;
  G  QV1-diss slow-manifold (quasi-static) drain law with K' != 0: closed form and its expansion in df = f - f0
     (new here: explains the track's '(H/m_q)^6' deviation table).
  H  2207.02826: graviton-equivalent coefficient of the Gamma = g~ |rdot/r| r_vac^2 source.
Output: ../results/v1_symbolic.json.  Run: python3 -I v1_symbolic.py
"""
import json, os
import sympy as sp

out = {}
ok = {}

# =============================================================================================== A
t = sp.symbols('t', real=True)
x, y, z = sp.symbols('x y z', real=True)
a = sp.Function('a', positive=True)(t)
q = sp.Function('q')(t)
Kf = sp.Function('K')
ef = sp.Function('eps')
X = [t, x, y, z]
g = sp.diag(-1, a**2, a**2, a**2)
gi = g.inv()
n = 4
Gam = [[[sp.simplify(sum(gi[l, s] * (sp.diff(g[s, m], X[nn]) + sp.diff(g[s, nn], X[m]) - sp.diff(g[m, nn], X[s]))
                         for s in range(n)) / 2) for nn in range(n)] for m in range(n)] for l in range(n)]


def ricci():
    Ric = sp.zeros(n, n)
    for m in range(n):
        for nn in range(n):
            expr = 0
            for l in range(n):
                expr += sp.diff(Gam[l][m][nn], X[l]) - sp.diff(Gam[l][m][l], X[nn])
                for s in range(n):
                    expr += Gam[l][l][s] * Gam[s][m][nn] - Gam[l][nn][s] * Gam[s][m][l]
            Ric[m, nn] = sp.simplify(expr)
    return Ric


Ric = ricci()
Rs = sp.simplify(sum(gi[m, nn] * Ric[m, nn] for m in range(n) for nn in range(n)))
Gmn = sp.simplify(Ric - Rs * g / 2)
H = sp.diff(a, t) / a
Rs_expected = 6 * (sp.diff(H, t) + 2 * H**2)
ok['A_R_MTW_sign'] = sp.simplify(Rs - Rs_expected) == 0

Kq = Kf(q)
grad = [sp.diff(Kq, X[m]) for m in range(n)]
hess = sp.Matrix(n, n, lambda m, nn: sp.diff(Kq, X[m], X[nn]) - sum(Gam[l][m][nn] * grad[l] for l in range(n)))
boxK = sp.simplify(sum(gi[m, nn] * hess[m, nn] for m in range(n) for nn in range(n)))
mu_x = sp.diff(ef(q), q) - sp.diff(Kf(q), q) * Rs / 2          # off-shell mu(t) := eps' - K' R/2
rhoV = ef(q) - mu_x * q
E = sp.simplify(Kq * Gmn + rhoV * g + (g * boxK - hess))      # BRIEF (g): K G = T - rhoV g - (g Box - nn)K => T = E
rho_E = sp.simplify(E[0, 0])                                   # T_00 = rho (u_0 = -1)
p_E = sp.simplify(E[1, 1] / a**2)
lhs = sp.diff(rho_E, t) + 3 * H * (rho_E + p_E)
rhs = q * sp.diff(mu_x, t)
resA = sp.simplify(sp.expand(lhs - rhs))
ok['A_continuity_equals_q_mudot'] = resA == 0
out['A'] = {'rho_E (Friedmann side)': str(sp.simplify(rho_E)),
            'identity': 'rho_E-dot + 3H(rho_E + p_E) - q * d(mu)/dt',
            'residual': str(resA),
            'meaning': 'with E = T^m on shell: rho_m-dot + 3H(rho_m+p_m) = q mu_eff-dot = q S (m-diss)'}
# rho_V-dot = (1/2) R K-dot - q S
resA2 = sp.simplify(sp.diff(rhoV, t) - (Rs / 2 * sp.diff(Kq, t) - q * sp.diff(mu_x, t)))
ok['A_rhoVdot'] = resA2 == 0
# total energy: rho_m + [rhoV - 3 H K-dot ...] obeys standard continuity -> Friedmann side checked by A
# the 00 component: 3 K H^2 + 3 H Kdot = rho_m + rhoV
res00 = sp.simplify(rho_E - (3 * Kq * H**2 + 3 * H * sp.diff(Kq, t) - rhoV))
ok['A_Friedmann_00'] = res00 == 0

# =============================================================================================== B
Gs, gam, E_, r, hh, M = sp.symbols('G gamma E r_V h M', positive=True)
# dimensionful (5.6) at constant G, radiation only: Hdot + 2H^2 = (16 pi G/3) rho_V ; rhoV-dot = -gamma |H| (Hdot+2H^2)^2
# scaled: rho = r E^4, H = h E, t = tau/E  ->  hdot + 2 h^2 = (16 pi G/3) E^2 r ;  r-dot = -gamma h (hdot+2h^2)^2
coef = sp.Rational(16, 3) * sp.pi * Gs * E_**2
Esol = sp.solve(sp.Eq(coef, 2), E_)
out['B'] = {'hdot+2h^2 in E units': str(coef) + ' * r_V',
            'E making it 2 r_V': str(Esol), 'E in M (G = 1/(8 pi M^2))': str([sp.simplify(e.subs(Gs, 1 / (8 * sp.pi * M**2))) for e in Esol])}
ok['B_Estar_sqrt3M'] = sp.simplify(Esol[0].subs(Gs, 1 / (8 * sp.pi * M**2)) - sp.sqrt(3) * M) == 0
# E_P units: G E_P^2 = 1 -> hdot+2h^2 = (16 pi/3) r.  Rescale r' = (8 pi/3) r so that hdot+2h^2 = 2 r';
# r'-dot = (8 pi/3) r-dot = -(8 pi/3) gamma h (2 r')^2  -> gamma_printed = (8 pi/3) gamma_phys
gp = sp.Rational(8, 3) * sp.pi
out['B']['gamma_printed/gamma_phys in E_P reading'] = str(gp)
out['B']['gamma_phys for printed 1/100 (E_P reading)'] = float(sp.Rational(1, 100) / gp)
# check drain in E* units: r-dot = -gamma h (2r)^2 = -4 gamma h r^2
out['B']['rV-dot (E* units)'] = '-4*gamma*h*r_V**2'

# =============================================================================================== C
xi, xi0, C1 = sp.symbols('xi xi0 C1', real=True)
rV = 1 / (4 * gam * (xi - xi0))
kk = 1 / (2 * gam * (xi - xi0))
hx = sp.sqrt((C1 + sp.exp(4 * xi0) * sp.Ei(4 * xi - 4 * xi0)) / (gam * sp.exp(4 * xi)))
c1 = sp.simplify(sp.diff(rV, xi) + 4 * gam * rV**2)
c2 = sp.simplify(sp.diff(kk, xi) + 2 * gam * kk**2)
c3 = sp.simplify(sp.diff(hx**2, xi) / 2 + 2 * hx**2 - kk)   # h h' = (h^2)'/2
c4 = sp.simplify(kk - 2 * rV)
ok['C_AppB_rV'] = c1 == 0
ok['C_AppB_k'] = c2 == 0
ok['C_AppB_h'] = c3 == 0
ok['C_AppB_k_eq_2rV'] = c4 == 0
# r_M/r_V = h^2/r_V - 1 at large xi with C1 = 0: Ei(y) ~ e^y/y (1 + 1/y + 2/y^2 + ...)
yv = sp.symbols('y', positive=True)
ratio_series = sp.series((sp.exp(yv) / yv * (1 + 1 / yv + 2 / yv**2)) * sp.exp(-yv) / (1 / yv) - 1, yv, sp.oo, 3)
out['C'] = {'h^2/r_V - 1 at large y=4(xi-xi0)': str(ratio_series), 'so': 'r_M/r_V = 1/(4(xi-xi0)) + O(xi^-2)'}
# vacuum dominated leading order: h = (4 gamma xi)^(-1/2), xi(tau) from d xi/d tau = h
tau = sp.symbols('tau', positive=True)
xit = (3 * tau / (4 * sp.sqrt(gam)))**sp.Rational(2, 3)
c5 = sp.simplify(sp.diff(xit, tau) - (4 * gam * xit)**sp.Rational(-1, 2))
c6 = sp.simplify((1 / (4 * gam * xit)) / (6 * gam * tau)**sp.Rational(-2, 3))
c7 = sp.simplify((4 * gam * xit)**sp.Rational(-1, 2) / (6 * gam * tau)**sp.Rational(-1, 3))
ok['C_asym_xi(tau)'] = c5 == 0
ok['C_asym_rV=(6 gamma tau)^-2/3'] = sp.simplify(c6 - 1) == 0
ok['C_asym_h=(6 gamma tau)^-1/3'] = sp.simplify(c7 - 1) == 0

# =============================================================================================== D
rho_m, p_m, rhoVs, Ms, N, kap, H1, rVi = sp.symbols('rho_m p_m rho_V M N kappa H_1 rho_Vi', positive=True)
Rexpr = (rho_m - 3 * p_m + 4 * rhoVs) / Ms**2                  # from trace of (g) at K' = 0
# trace check: -K R = T - 4 rhoV  with T = -rho + 3p
ok['D_trace_R'] = sp.simplify(-Ms**2 * Rexpr - ((-rho_m + 3 * p_m) - 4 * rhoVs)) == 0
drain_perN = gam / 36 * Rexpr**2                               # |H| dt = dN for H > 0
rad = drain_perN.subs({p_m: rho_m / 3})
ok['D_radiation_drain_indep_of_rho_r'] = sp.simplify(rad - sp.Rational(4, 9) * gam * rhoVs**2 / Ms**4) == 0
f_ = sp.Function('f')
sol = sp.dsolve(sp.Eq(f_(N).diff(N), -sp.Rational(4, 9) * gam * f_(N)**2 / Ms**4), f_(N), ics={f_(0): rVi})
closed = sp.simplify(sol.rhs)
ok['D_closed_form'] = sp.simplify(closed - rVi / (1 + sp.Rational(4, 9) * gam * rVi * N / Ms**4)) == 0
# plasma with trace fraction kappa during radiation era: R = kappa rho/M^2, rho = 3 M^2 H^2, H = H1 e^{-2N}
Hh = sp.symbols('H', positive=True)
perN_trace = sp.simplify(gam / 36 * (kap * 3 * Ms**2 * Hh**2 / Ms**2)**2)
ok['D_trace_perN=gamma kappa^2 H^4/4'] = sp.simplify(perN_trace - gam * kap**2 * Hh**4 / 4) == 0
Dtrace = sp.integrate(perN_trace.subs(Hh, H1 * sp.exp(-2 * N)), (N, 0, sp.oo))
Ddust = sp.integrate((gam / 36 * (3 * Ms**2 * Hh**2 / Ms**2)**2).subs(Hh, H1 * sp.exp(-sp.Rational(3, 2) * N)), (N, 0, sp.oo))
ok['D_trace_total=gamma kappa^2 H1^4/32'] = sp.simplify(Dtrace - gam * kap**2 * H1**4 / 32) == 0
ok['D_dust_total=gamma H1^4/24'] = sp.simplify(Ddust - gam * H1**4 / 24) == 0
# in E* units with dust: d r_V/dN = -(gamma/4)(r_d e^{-3N} + 4 r_V)^2 ; leading dust boundary gamma r_d^2/24
rd = sp.symbols('r_d', positive=True)
Rstar = sp.simplify(Rexpr.subs({p_m: 0, rho_m: 9 * Ms**4 * rd, rhoVs: 9 * Ms**4 * r}) / (3 * Ms**2))  # R / E*^2
perN_star = sp.simplify(gam / 36 * (Rstar * 3 * Ms**2)**2 / (9 * Ms**4))
ok['D_Estar_dust_drain'] = sp.simplify(perN_star - gam / 4 * (rd + 4 * r)**2) == 0
lead = sp.integrate(gam / 4 * (rd * sp.exp(-3 * N))**2, (N, 0, sp.oo))
ok['D_dust_boundary_leading=gamma rd^2/24'] = sp.simplify(lead - gam * rd**2 / 24) == 0
# inflaton plateau: d rhoV/dN = -(4 gamma/9) V^2/M^4
Vv = sp.symbols('V', positive=True)
ok['D_inflaton'] = sp.simplify(drain_perN.subs({rho_m: 0, p_m: 0, rhoVs: Vv}) - sp.Rational(4, 9) * gam * Vv**2 / Ms**4) == 0
# attractor in physical units
tt = sp.symbols('t', positive=True)
rho_att = 9 * Ms**4 * (6 * gam * sp.sqrt(3) * Ms * tt)**sp.Rational(-2, 3)
H_att = sp.sqrt(3) * Ms * (6 * gam * sp.sqrt(3) * Ms * tt)**sp.Rational(-1, 3)
ok['D_attractor_Friedmann'] = sp.simplify(3 * Ms**2 * H_att**2 - rho_att) == 0
ok['D_attractor_drain'] = sp.simplify(sp.diff(rho_att, tt) + gam * H_att * (2 * rho_att / (3 * Ms**2))**2) == 0
t_of_H = sp.solve(sp.Eq(H_att, Hh), tt)
ok['D_t(H)=M^2/(2 gamma H^3)'] = sp.simplify(t_of_H[0] - Ms**2 / (2 * gam * Hh**3)) == 0
Natt = sp.integrate(H_att, (tt, 0, tt))
ok['D_N=(3 sqrt3 M t/(4 sqrt gamma))^(2/3)'] = sp.simplify(Natt - (3 * sp.sqrt(3) * Ms * tt / (4 * sp.sqrt(gam)))**sp.Rational(2, 3)) == 0
ok['D_Ht=(2/3)N'] = sp.simplify(H_att * tt - sp.Rational(2, 3) * Natt) == 0

# =============================================================================================== E
# KSV16: BRIEF equations with K = f/3 in units q0 = chi0 = 1, generic eps(f), u(tau) = mu_eff, matter w
ta = sp.symbols('tau', real=True)
f = sp.Function('f')(ta); h = sp.Function('h')(ta); u = sp.Function('u')(ta); rM = sp.Function('r_M')(ta)
w = sp.symbols('w')
eps = sp.Function('eps')
K = f / 3
Rf = 6 * (h.diff(ta) + 2 * h**2)
rVt = eps(f) - u * f
qeq = sp.Eq(u, eps(f).diff(f) - sp.Rational(1, 2) * sp.Rational(1, 3) * Rf)       # (q) with mu -> u(tau)
fried = 3 * K * h**2 + 3 * h * K.diff(ta) - (rM + rVt)                          # 00 component
# space-space component (derived from (g) in FRW): -K(2 hdot + 3 h^2) = p - rhoV + Kddot + 2 h Kdot
ii = -K * (2 * h.diff(ta) + 3 * h**2) - (w * rM - rVt + K.diff(ta, 2) + 2 * h * K.diff(ta))
# check (17b): h fdot + f h^2 = rV~ + r_M  with rV~ = eps - f eps' + f(hdot + 2h^2)
rVtilde_ksv = eps(f) - f * eps(f).diff(f) + f * (h.diff(ta) + 2 * h**2)
ok['E_rVtilde_ksv_eq_eps-uf'] = sp.simplify((rVtilde_ksv - rVt).subs(u, qeq.rhs)) == 0
ok['E_17b'] = sp.simplify((fried - (h * f.diff(ta) + f * h**2 - rVtilde_ksv - rM)).subs(u, qeq.rhs)) == 0
# (17a): hdd + 4 h hdot = fdot eps'' - s  where s = du/dtau
s_ = sp.symbols('s')
e17a = sp.simplify(qeq.rhs.diff(ta) - (f.diff(ta) * eps(f).diff(f, 2) - h.diff(ta, 2) - 4 * h * h.diff(ta)))
ok['E_17a_(u-dot = s)'] = e17a == 0
# constraint propagation: C := fried; use ii-eq for fddot, u-dot = s, rM-dot = -3h(1+w) rM + f s
C = fried
dC = C.diff(ta)
fdd = sp.solve(ii, f.diff(ta, 2))[0]
sub = {f.diff(ta, 2): fdd, u.diff(ta): s_, rM.diff(ta): -3 * h * (1 + w) * rM + f * s_}
dC = dC.subs(sub)
# eliminate u via (q)?  No: u is independent (state); hdot from (q): hdot = eps' - u - 2 h^2
hdot_q = eps(f).diff(f) - u - 2 * h**2
dC2 = sp.simplify(dC.subs(h.diff(ta), hdot_q))
C2 = sp.simplify(C.subs(h.diff(ta), hdot_q))
ratio = sp.simplify(dC2 + 3 * h * C2)
ok['E_constraint_dC=-3hC'] = ratio == 0
# without the matter source f s, the constraint is violated at rate f s
dC_nosrc = sp.simplify(C.diff(ta).subs({f.diff(ta, 2): fdd, u.diff(ta): s_, rM.diff(ta): -3 * h * (1 + w) * rM}).subs(h.diff(ta), hdot_q))
ok['E_no_source_violation=f s'] = sp.simplify(dC_nosrc + 3 * h * C2 - f * s_) == 0
out['E'] = {'dC + 3hC without matter source': str(sp.simplify(dC_nosrc + 3 * h * C2))}
# de Sitter family and KSV u0
fs_ = sp.symbols('f', positive=True)
epsK = sp.Rational(1, 2) * (-fs_**2 + fs_**4 / 3)
u0 = sp.diff(epsK, fs_).subs(fs_, 1)
ok['E_u0=-1/3'] = sp.simplify(u0 + sp.Rational(1, 3)) == 0
out['E']['u0'] = str(u0)
# fdd derived here (to compare with any implementation)
out['E']['fddot (from (g)_ii, K=f/3)'] = str(sp.simplify(fdd))

# =============================================================================================== F
lam, mu_, xs = sp.symbols('lambda mu xi', real=True)
rV22 = lam + xs**2 - mu_ * xs
xb = (mu_ + sp.sqrt(mu_**2 - 4 * lam)) / 2
ok['F_rV(xibar)=0'] = sp.simplify(rV22.subs(xs, xb)) == 0
viol = sp.simplify(sp.diff(lam + xs**2, xs).subs(xs, xb) - mu_)
out['F'] = {'eps\'(xibar) - mu': str(viol),
            'mu=4, lambda=1': float(viol.subs({mu_: 4, lam: 1})), 'mu=4, lambda=-1': float(viol.subs({mu_: 4, lam: -1})),
            'xibar(mu=4, lambda=1)': float(xb.subs({mu_: 4, lam: 1})), 'xibar(mu=4, lambda=-1)': float(xb.subs({mu_: 4, lam: -1}))}
# ad hoc asymptote: r-dot = -g h^5, 3 h^2 = 8 pi r
gt = sp.symbols('g', positive=True)
rr = sp.Function('r')
solF = sp.dsolve(sp.Eq(rr(tau).diff(tau), -gt * (8 * sp.pi * rr(tau) / 3)**sp.Rational(5, 2)), rr(tau))
cR = (sp.Rational(3, 2))**sp.Rational(-2, 3) * (8 * sp.pi / 3)**sp.Rational(-5, 3)
test = cR * (gt * tau)**sp.Rational(-2, 3)
ok['F_adhoc_asymptote'] = sp.simplify(test.diff(tau) + gt * (8 * sp.pi * test / 3)**sp.Rational(5, 2)) == 0
out['F']['r_V coefficient'] = float(cR)
out['F']['h coefficient'] = float(sp.sqrt(8 * sp.pi / 3 * cR))
# quasi-static r_M on the ad hoc attractor: rM-dot + 4 h rM = g h^5 -> rM ~ g h^4/4 when h tau >> 1
# =============================================================================================== G
# QV1-diss quasi-static drain law (KSV units, K = f/3).  Static de Sitter branch: f h^2 = f eps' - eps, u = 2 eps/f - eps'.
# First order in the drift s: fdot = s / U'(f), U = 2 eps/f - eps'.  rV~-dot = -f s + k fdot, k = 2 h^2, f s = gamma h k^2.
fq, gg = sp.symbols('f gamma', positive=True)
E0 = sp.Function('eps')
h2 = sp.diff(E0(fq), fq) - E0(fq) / fq
U = 2 * E0(fq) / fq - sp.diff(E0(fq), fq)
k2 = 2 * h2
fs_drain = gg * sp.sqrt(h2) * k2**2                 # = f s
fdot_ = (fs_drain / fq) / sp.diff(U, fq)
rate_tau = -fs_drain + k2 * fdot_
rate_N = sp.simplify(rate_tau / sp.sqrt(h2))
rV_branch = fq * h2                                 # gravitating rV on the static branch (r_M = O(s))
closed_G = sp.simplify(rate_N / (-4 * gg * rV_branch**2))
out['G'] = {'d rV/dN / (-4 gamma rV^2)  (quasi-static, first order in drift)': str(closed_G)}
# compare with eps''/(f^2 eps'' - 2 rV)
target = sp.diff(E0(fq), fq, 2) / (fq**2 * sp.diff(E0(fq), fq, 2) - 2 * rV_branch)
ok['G_closed_form'] = sp.simplify(closed_G - target) == 0
# expansion about f0 where rV(f0) = 0, i.e. f0 eps'(f0) = eps(f0):  generic Taylor data e1 = eps'(f0) etc.
f0, d, e2, e3, e4 = sp.symbols('f0 delta e2 e3 e4')
e1 = sp.symbols('e1')
epsT = e1 * f0 + e1 * d + e2 * d**2 / 2 + e3 * d**3 / 6 + e4 * d**4 / 24   # eps(f0) = f0 eps'(f0) = f0 e1
F = f0 + d
epsF = epsT
eps1 = sp.diff(epsT, d); eps2 = sp.diff(epsT, d, 2)
rVd = F * eps1 - epsF
expr = eps2 / (F**2 * eps2 - 2 * rVd) * f0**2 - 1
ser = sp.series(expr, d, 0, 3).removeO()
ser = sp.simplify(sp.expand(ser))
out['G']['f0^2 * rate/(-4 gamma rV^2) - 1  (series in delta = f - f0)'] = str(ser)
ok['G_first_order_cancels'] = sp.simplify(ser.coeff(d, 1)) == 0
c2 = sp.simplify(ser.coeff(d, 2))
out['G']['coefficient of delta^2'] = str(c2)
# for KSV eps = lambda + (1/2)(-f^2 + f^4/3) at lambda = 0: f0 = 1, e2 = 1, e3 = 4  -> -4 delta^2; delta ~ h^2/e2 ...
c2_ksv0 = c2.subs({f0: 1, e2: 1, e3: 4})
out['G']['lambda=0 (f0=1,e2=1,e3=4,m_q^2=2)'] = str(c2_ksv0) + ' * delta^2, delta = rV/(f0^2 e2) + ... = h^2 + ...'
# also lambda = 1: f0 = sqrt2, eps'' = -1 + 2 f^2 = 3, eps''' = 4 f = 4 sqrt2
c2_ksv1 = sp.nsimplify(c2.subs({f0: sp.sqrt(2), e2: 3, e3: 4 * sp.sqrt(2)}))
out['G']['lambda=1 (f0=sqrt2,e2=3,e3=4 sqrt2)'] = str(sp.simplify(c2_ksv1)) + ' * delta^2'

# =============================================================================================== H
# Graviton-equivalent coefficient of 2207.02826's Gamma = g~ |rdot/r| r_vac^2 (E_P-type units, 3 h^2 = 8 pi (...)):
# with radiation: hdot + 2 h^2 = (16 pi/3) r_vac, so gamma |h| (hdot+2h^2)^2 = gamma (16 pi/3)^2 |h| r_vac^2.
ge = sp.simplify((sp.Rational(16, 3) * sp.pi)**2 / (32 * sp.pi**2))
out['H'] = {'g~ equivalent to gamma = 1/(32 pi^2)': str(ge), 'value': float(ge),
            'ratio 2e11 / g~_grav': float(2e11 / ge), 'log10': float(sp.log(2e11 / ge, 10))}

out['checks'] = {k_: bool(v_) for k_, v_ in ok.items()}
out['all_pass'] = all(out['checks'].values())
os.makedirs(os.path.join(os.path.dirname(__file__), '..', 'results'), exist_ok=True)
with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'v1_symbolic.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out['checks'], indent=1))
print('ALL PASS' if out['all_pass'] else 'SOME FAIL')
print(json.dumps(out['G'], indent=1))
print(json.dumps(out['H'], indent=1))
