"""Track `dissipation`, script 1: symbolic checks of the dissipative q-theory equations.

Conventions (BRIEF.md): signature (-,+,+,+), hbar = c = 1, K = M^2 = 1/(8 pi G) for constant G,
R = 6 (Hdot + 2 H^2) in flat FRW (MTW sign; KV use Weinberg's sign, R_W = -R, irrelevant for R^2).

Checks (all written to ../results/symbolic_checks.json):
  A. FRW form of the off-shell identity  nabla^m E_mn = -q d_n mu(x):  with mu(t) := eps'(q) - K'(q) R/2,
     d/dt E_rho + 3H (E_rho + E_p) = q dmu/dt.  So a dissipative q-equation dmu_eff/dt = S is consistent
     with the unchanged Einstein equation iff matter obeys rho_m' + 3H(rho_m + p_m) = q S.
  B. 1601.00601 (constant G): dimensionless system (6.1) <-> dimensionful (5.3), (5.6); unit E_* = sqrt(3) M.
     App. B exact solution and the asymptote (6 gamma tau)^(-2/3).
  C. 1601.04676: eqs (17a-c), (18a), (19)-(21), (23) from the brief's equations with K = f q0/3 (sigma = 1 units).
  D. 2207.03453: Gamma = 0 limit of (15b) admits every constant xi; endpoint violates the action's
     n-equation by sqrt(mu^2 - 4 lambda); exact asymptotic coefficients of (18) and r_M.
  E. Relaxation laws per e-fold in brief conventions; radiation closed form; dust fixed point.
Run:  python3 -I s1_symbolic.py   (from this directory; < 1 min)
"""
import json, os
import sympy as sp

out = {}
t = sp.symbols('t', real=True)
a = sp.Function('a')(t)
qf = sp.Function('q')(t)
Kf = sp.Function('K')
ef = sp.Function('eps')
H = sp.diff(a, t) / a
R = 6 * (sp.diff(H, t) + 2 * H**2)

# ---------------- A. FRW identity with mu(t) defined off shell ----------------
K = Kf(qf)
Kd = sp.diff(K, t)
Kdd = sp.diff(K, t, 2)
mu_x = sp.diff(ef(qf), qf) - sp.diff(K, qf) * R / 2          # off-shell mu(x)
rhoV = ef(qf) - mu_x * qf
E_rho = 3 * H**2 * K + 3 * H * Kd - rhoV                      # = rho_m on shell
E_p = -(2 * sp.diff(H, t) + 3 * H**2) * K - Kdd - 2 * H * Kd + rhoV   # = p_m on shell
lhs = sp.diff(E_rho, t) + 3 * H * (E_rho + E_p)
rhs = qf * sp.diff(mu_x, t)
resA = sp.simplify(sp.expand(lhs - rhs))
out['A_FRW_identity_dErho_plus_3H_Erho_plus_Ep_equals_q_mudot'] = {
    'residual': str(resA), 'holds': resA == 0,
    'consequence': 'rho_m-dot + 3H(rho_m+p_m) = q * d(mu_eff)/dt: dissipated vacuum energy goes into matter'}
# with a constant-mu E (mu treated as a number): identity becomes (mu(x)-mu) qdot
mu_c = sp.symbols('mu_c')
rhoV_c = ef(qf) - mu_c * qf
E_rho_c = 3 * H**2 * K + 3 * H * Kd - rhoV_c
E_p_c = -(2 * sp.diff(H, t) + 3 * H**2) * K - Kdd - 2 * H * Kd + rhoV_c
lhs_c = sp.diff(E_rho_c, t) + 3 * H * (E_rho_c + E_p_c)
resA2 = sp.simplify(sp.expand(lhs_c + (mu_x - mu_c) * sp.diff(qf, t)))
out['A_constant_mu_form'] = {'residual_of_[dErho+3H(Erho+Ep)] = -(mu(x)-mu) qdot': str(resA2), 'holds': resA2 == 0}

# ---------------- B. 1601.00601 constant-G dimensionless system ----------------
M, gam, Es = sp.symbols('M gamma E_star', positive=True)
Hs, Hd, rhoV_s, rhor_s = sp.symbols('H Hdot rho_V rho_r', real=True)
# Friedmann + Raychaudhuri for radiation + vacuum, constant K = M^2
fr = sp.Eq(3 * Hs**2 * M**2, rhoV_s + rhor_s)
ray = sp.Eq(Hd, -(rhor_s * sp.Rational(4, 3)) / (2 * M**2))
sol = sp.solve([fr, ray], [rhor_s, Hd], dict=True)[0]
k_dimful = sp.simplify(sol[Hd] + 2 * Hs**2)        # Hdot + 2H^2 = R/6
out['B_Hdot_plus_2H2_radiation'] = str(k_dimful)  # expect 2 rho_V/(3 M^2)
# dimensionless map: r = rho/E*^4, h = H/E*, tau = E* t with E*^2 = 3 M^2
rV, rr, h, hd = sp.symbols('r_V r_r h hdot', real=True)
subsB = {rhoV_s: rV * Es**4, Hs: h * Es}
k_dimless = sp.simplify((k_dimful.subs(subsB) / Es**2).subs(Es, sp.sqrt(3) * M))
out['B_hdot_plus_2h2_in_Estar_units'] = str(k_dimless)   # expect 2 r_V   (KV 5.3b)
# drain: rho_V-dot = -gamma |H| (Hdot+2H^2)^2 (KV 5.6) -> r_V-dot = -gamma |h| (hdot+2h^2)^2 in E* units
drain_dimful = -gam * Hs * k_dimful**2
drain_dimless = sp.simplify((drain_dimful.subs(subsB) / Es**5).subs(Es, sp.sqrt(3) * M))
out['B_rV_dot_in_Estar_units'] = str(drain_dimless)  # expect -4 gamma h r_V^2
# same in E_P units (G = 1, i.e. M^2 = 1/(8 pi), r = rho/E_P^4): hdot+2h^2 = (16 pi/3) r_V, as 2207.03453 eq (9)
k_EP = sp.simplify(k_dimful.subs({rhoV_s: rV, Hs: h, M: 1 / sp.sqrt(8 * sp.pi)}))
out['B_hdot_plus_2h2_in_EP_units'] = str(k_EP)
out['B_note'] = ('KV 1601.00601 (6.1) [hdot+2h^2 = 2 r_V and -rV-dot = gamma|h|(hdot+2h^2)^2] equals (5.6) '
                 'with the same gamma only in units E_* = (3/(8 pi G))^(1/2) = sqrt(3) M = 4.22e18 GeV; '
                 'in E_P (G=1) units, as in 2207.03453 eq (7c), gamma_(6.1) = (8 pi/3) gamma_phys.')
# App. B exact solution: dr/dxi = -4 gamma r^2, k = 2 r, h dh/dxi + 2 h^2 = k
xi, xi0, C1 = sp.symbols('xi xi0 C1', real=True)
r_sol = 1 / (4 * gam * (xi - xi0))
chk1 = sp.simplify(sp.diff(r_sol, xi) + 4 * gam * r_sol**2)
h2_sol = (C1 + sp.exp(4 * xi0) * sp.Ei(4 * xi - 4 * xi0)) / (gam * sp.exp(4 * xi))
# h dh/dxi = (1/2) d(h^2)/dxi
chk2 = sp.simplify(sp.Rational(1, 2) * sp.diff(h2_sol, xi) + 2 * h2_sol - 2 * r_sol)
out['B_AppB_rV_solution_residual'] = str(chk1)
out['B_AppB_h_solution_residual'] = str(chk2)
# asymptote: h^2 ~ r_V (1 + 1/(4 xi) + ...): Ei(y) e^-y ~ 1/y + 1/y^2 + ...
y = sp.symbols('y', positive=True)
import mpmath as mp
ei_chk = {str(Y): float(mp.e**(-Y) * mp.ei(Y) / (1 / mp.mpf(Y) + 1 / mp.mpf(Y)**2 + 2 / mp.mpf(Y)**3)) for Y in (50, 200, 1000)}
out['B_Ei_asymptotic_check_ratio_to_1/y+1/y^2+2/y^3'] = ei_chk   # -> 1
# tau(xi) on the leading asymptote: dxi/dtau = h = (4 gamma xi)^(-1/2) => xi = (3 tau/(4 sqrt(gamma)))^(2/3)
tau = sp.symbols('tau', positive=True)
xi_tau = (3 * tau / (4 * sp.sqrt(gam)))**sp.Rational(2, 3)
chk3 = sp.simplify(sp.diff(xi_tau, tau) - (4 * gam * xi_tau)**sp.Rational(-1, 2))
rV_tau = sp.simplify(1 / (4 * gam * xi_tau))
ratio = sp.simplify(rV_tau / (6 * gam * tau)**sp.Rational(-2, 3))
out['B_asymptote'] = {'dxi/dtau - h residual': str(chk3), 'r_V(tau)': str(rV_tau),
                      'r_V / (6 gamma tau)^(-2/3)': str(ratio)}
out['B_radiation_fraction_on_attractor'] = 'r_M/r_V = 1/(4 xi) + O(xi^-2)  (from the Ei series): Omega_M -> 0'

# ---------------- C. 1601.04676 from the brief's equations ----------------
f = sp.Function('f')(t)
hh = sp.Function('h')(t)
u = sp.Function('u')(t)
rM = sp.Function('r_M')(t)
w, gf, gh = sp.symbols('w gamma_f gamma_h', real=True)
eps = sp.Function('epsilon')
# units: q = f q0, K = f q0/3 (sigma = 1), rho = r/chi0, H = h/sqrt(chi0 q0), t = tau sqrt(chi0 q0)
# => (dimensionless) Friedmann: f h^2 + h fdot = r_V + r_M ; q-eq: eps'(f) - (hdot + 2h^2) = u (K' = 1/3, R/2*K' = hdot+2h^2)
Kdl = f / 3
friedmann = 3 * hh**2 * Kdl + 3 * hh * sp.diff(Kdl, t) - (eps(f) - u * f) - rM
qeq_u = sp.diff(eps(f), f) - (sp.diff(hh, t) + 2 * hh**2)      # = u on shell
rVt = eps(f) - sp.diff(eps(f), f) * f + f * (sp.diff(hh, t) + 2 * hh**2)
out['C_rVtilde_equals_eps_minus_ueff_f'] = str(sp.simplify(rVt - (eps(f) - qeq_u * f)))
# (17a)/(21): d u_eff/dtau = s  <=>  hdd + 4 h hdot = fdot eps'' - s
s = sp.symbols('s')
eq8a = sp.diff(qeq_u, t) - s
eq8a_form = sp.diff(hh, t, 2) + 4 * hh * sp.diff(hh, t) - (sp.diff(f, t) * sp.diff(eps(f), f, 2) - s)
out['C_17a_equivalence_residual'] = str(sp.simplify(eq8a + eq8a_form))
# (17c) from the identity: r_M-dot + 3h(1+w) r_M = f s ; check constraint propagation of (17b) with (21),(17c),
# and the space-space equation used to evolve f.
fd, fdd_, hd_, hdd_ = sp.symbols('fd fdd hd hdd')
rMs, us_, fs_, hs_ = sp.symbols('rM u f h')
eps_s = sp.Function('epsilon')
# space-space: -(2 hdot + 3 h^2) K = p_M - r_V + Kdd + 2 h Kdot   (K = f/3, r_V = eps - u f)
ss = -(2 * hd_ + 3 * hs_**2) * fs_ / 3 - (w * rMs - (eps_s(fs_) - us_ * fs_) + fdd_ / 3 + 2 * hs_ * fd / 3)
fdd_sol = sp.solve(ss, fdd_)[0]
hd_sol = sp.diff(eps_s(fs_), fs_) - us_ - 2 * hs_**2           # q-eq
ud = s                                                         # dissipative q-eq
rMd = -3 * hs_ * (1 + w) * rMs + fs_ * s                      # matter with source f s
Cexpr = fs_ * hs_**2 + hs_ * fd - (eps_s(fs_) - us_ * fs_) - rMs
# total derivative of C along the flow
dC = (sp.diff(Cexpr, fs_) * fd + sp.diff(Cexpr, hs_) * hd_sol + sp.diff(Cexpr, fd) * fdd_sol.subs(hd_, hd_sol)
      + sp.diff(Cexpr, us_) * ud + sp.diff(Cexpr, rMs) * rMd)
lam_C = sp.simplify(dC / Cexpr)
out['C_constraint_propagation_dC_over_C'] = str(sp.simplify(lam_C))
# without the matter source the constraint is violated at rate f s:
dC_nosrc = dC - sp.diff(Cexpr, rMs) * fs_ * s
out['C_constraint_without_matter_source_dC_minus_kC'] = str(sp.simplify(dC_nosrc - lam_C * Cexpr))
# de Sitter family (10): f, h const, r_M = 0
hbar2 = sp.solve(sp.Eq(fs_ * hs_**2, eps_s(fs_) - (sp.diff(eps_s(fs_), fs_) - 2 * hs_**2) * fs_), hs_**2)
out['C_deSitter_h2'] = str(hbar2)  # expect eps' - eps/f
# KSV eps and u0
ff = sp.symbols('f', positive=True)
lamb = sp.symbols('lambda', real=True)
eps_KSV = sp.Rational(1, 2) * (-ff**2 + ff**4 / 3)
u0 = sp.diff(eps_KSV, ff).subs(ff, 1)
out['C_KSV_u0'] = str(u0)  # -1/3
eps_l = lamb + eps_KSV
f0sq = sp.solve(sp.Eq(eps_l - ff * sp.diff(eps_l, ff), 0), ff**2) if False else None
X = sp.symbols('X', positive=True)
f0sq = sp.solve(sp.Eq((eps_l - ff * sp.diff(eps_l, ff)).subs(ff, sp.sqrt(X)), 0), X)
out['C_KSV_with_lambda_f0^2_roots'] = [str(r) for r in f0sq]
out['C_note_eq13_vs_18b'] = ('1601.04676 eq (13) writes q S = Gamma_q qdot^2 + Gamma_H Hdot^2, but the integrated '
                            'eqs (18b), (21) put s~ = gamma_f fdot^2 + gamma_h hdot^2 = du_eff/dtau and f s~ into '
                            'matter (17c): the two differ by a factor f (=q/q0). Notational; the numerics follow (17)-(21).')

# ---------------- D. 2207.03453 brane ----------------
xs, mus, lam_s, gt = sp.symbols('xi mu lambda gamma_tilde', real=True)
rV_brane = lam_s + xs**2 - mus * xs
xibar = (mus + sp.sqrt(mus**2 - 4 * lam_s)) / 2
out['D_rV_at_xibar'] = str(sp.simplify(rV_brane.subs(xs, xibar)))
viol = sp.simplify(sp.diff(lam_s + xs**2, xs).subs(xs, xibar) - mus)
out['D_action_n_equation_violation_eps_prime_minus_mu_at_endpoint'] = str(viol)  # sqrt(mu^2-4 lambda)
out['D_note'] = ('Action (1) gives eps\'(q) = mu algebraically (eq 2a). Eq (15b) xi-dot r_V\'(xi) = -Gamma '
                 'at Gamma = 0 admits every constant xi, i.e. it deforms the weaker dr_V/dt = 0 (eq 6), '
                 'not eq (2a). The late-time endpoint xi-bar violates eq (2a) by sqrt(mu^2 - 4 lambda) '
                 '(= 2 sqrt(3) for mu=4, lambda=1), unless lambda = mu^2/4.')
# exact asymptote for Gamma = g~ |h| h^4 with 3 h^2 = 8 pi r_V (r_M negligible)
c = 8 * sp.pi / 3
rr_ = sp.symbols('r', positive=True)
# dr/dtau = -g~ (c r)^(5/2)  => r^(-3/2) = (3/2) g~ c^(5/2) tau
r_asym = ((sp.Rational(3, 2) * gt * c**sp.Rational(5, 2) * tau))**sp.Rational(-2, 3)
coef_r = sp.simplify(r_asym * (gt * tau)**sp.Rational(2, 3))
h_asym = sp.sqrt(c * r_asym)
coef_h = sp.simplify(h_asym * (gt * tau)**sp.Rational(1, 3))
out['D_asymptote_rV_coefficient_of_(gt tau)^(-2/3)'] = [str(coef_r), float(coef_r.subs(gt, 1))]
out['D_asymptote_h_coefficient_of_(gt tau)^(-1/3)'] = [str(coef_h), float(coef_h.subs(gt, 1))]
out['D_asymptote_rM'] = 'r_M ~ g~ h^4 / 4  (quasi-static balance 4 h r_M = g~ h^5, valid while h tau >> 1)'
# with the particle-production Gamma (8): r_V(chi) = 1/((16pi/3)^2 gamma (chi - chi0))  [2207 eq (12)]
chi = sp.symbols('chi')
rch = 1 / ((16 * sp.pi / 3)**2 * gam * (chi - xi0))
out['D_eq12_residual'] = str(sp.simplify(sp.diff(rch, chi) + (16 * sp.pi / 3)**2 * gam * rch**2))

# ---------------- E. relaxation laws per e-fold (brief conventions, constant K = M^2) ----------------
rho_m, p_m, rv, N = sp.symbols('rho_m p_m rho_V N', real=True)
Rsc = (rho_m - 3 * p_m + 4 * rv) / M**2      # trace of K G = T - rho_V g
drain = -gam / 36 * Rsc**2                 # d rho_V / dN  (H > 0)
out['E_drho_V_dN'] = str(sp.factor(drain))
rv_rad = sp.Function('rv')(N)
ode = sp.Eq(sp.diff(rv_rad, N), -gam / 36 * (4 * rv_rad / M**2)**2)
solE = sp.dsolve(ode, rv_rad, ics={rv_rad.subs(N, 0): sp.Symbol('rho_Vi', positive=True)})
out['E_radiation_closed_form'] = str(sp.simplify(solE.rhs))   # rho_Vi/(1 + (4 gamma/9) rho_Vi N / M^4)
out['E_dust_fixed_point'] = 'drain vanishes at R = 0, i.e. rho_V = -(rho_m - 3 p_m)/4 (= -rho_d/4 for dust) < 0'
# matter-era drain from H1 to infinity with rho_V ~ 0: d rho_V/dt = -(gamma/36) H R^2, R = 3 H^2 (dust dominated)
t1 = sp.symbols('t1', positive=True)
Hm = 2 / (3 * t)
integral = sp.integrate((gam / 36) * Hm * (3 * Hm**2)**2, (t, t1, sp.oo))
H1 = 2 / (3 * t1)
out['E_matter_era_drain_in_units_of_H1^4'] = str(sp.simplify(integral / H1**4))   # gamma/24
# radiation era with trace fraction kappa: (gamma/4) kappa^2 int H^4 dN, H ~ a^-2 => H_1^4/8
kap = sp.symbols('kappa', positive=True)
out['E_radiation_trace_drain_in_units_of_H1^4'] = str(sp.simplify(gam / 36 * 9 * kap**2 / 8))  # gamma kappa^2/32

checks = {
    'A_identity': resA == 0,
    'A_const_mu_form': resA2 == 0,
    'B_k_equals_2rV': sp.simplify(k_dimless - 2 * rV) == 0,
    'B_drain_Estar': sp.simplify(drain_dimless + 4 * gam * h * rV**2) == 0,
    'B_AppB_rV': chk1 == 0, 'B_AppB_h': chk2 == 0, 'B_asym_ratio_1': sp.simplify(ratio - 1) == 0,
    'C_rVtilde': sp.simplify(rVt - (eps(f) - qeq_u * f)) == 0,
    'C_17a': sp.simplify(eq8a + eq8a_form) == 0,
    'C_constraint_dC=-3hC': sp.simplify(lam_C + 3 * hs_) == 0,
    'C_no_source_violates_by_fs': sp.simplify(dC_nosrc - lam_C * Cexpr - fs_ * s) == 0,
    'B_Ei_asym': all(abs(v - 1) < 1e-3 for v in ei_chk.values()),
    'D_endpoint_zero': sp.simplify(rV_brane.subs(xs, xibar)) == 0,
    'D_eq12': sp.simplify(sp.diff(rch, chi) + (16 * sp.pi / 3)**2 * gam * rch**2) == 0,
}
out['checks'] = {k: bool(v) for k, v in checks.items()}
out['all_checks_pass'] = all(out['checks'].values())
os.makedirs('../results', exist_ok=True)
with open('../results/symbolic_checks.json', 'w') as fh:
    json.dump(out, fh, indent=1, default=str)
for k, v in out.items():
    print(k, ':', v)
