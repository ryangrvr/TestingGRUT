"""Track `bound`, symbolic part (sympy). Writes results/adiabatic_relation.json.

Derives, for the action S = int sqrt(-g)[K(q) R/2 - eps(q)] + S_m at fixed mu (= f(R) gravity, equations track 2.2):

 S1  First-order (in 1/m_q^2) effective dark-energy density of ANY f(R), written as
        f(R) = K_lab (R - 2 Lam/K_lab) + F(R),   F small,
        rho_DE,eff = 3 K_lab H^2 - rho_m - rho_r = Lam + (R F' - F)/2 - 3 F' H^2 - 3 H F'' Rdot   (on the GR background)
     and its value for the unscreened q-mode (F = f2 R^2/2, f2 = K/(3 m^2)):
        rho_dyn = [ (9/4) rho_d^2 + 2 rho_d rho_r - 4 Lam rho_r ] / (3 K m^2)
     i.e. Omega_dyn(t0) = (9/4) Omega_m^2 (H0/m)^2 (radiation neglected), independent of Lam; w_dyn = +1 in the dust era.
 S2  The same from the q-variables directly (generic Taylor data of eps, K about the R = 0 point of the given mu).
 S3  Lam_cal(mu) series: static (dS/AdS) root as a power series in dmu = mu - mu0.
 S4  The map (eps, K, mu) -> f(R) and its inverse:  f' = K(q_R), f'' = K'^2/(2D), d f/d mu = 2 q_R, f(0) = 2 eps*(mu);
     the canonical inverse  K(q) = f'(lam (q - q0)),  eps(q) = [K(q) lam (q-q0) - f(lam(q-q0))]/2 + mu q  solves the q-equation
     identically for ANY f; P_q  <=>  f(0) = 0, f'(0) > 0, f''(0) > 0 (+ final-state reading: R -> 0 is the attractor).
 S5  Threshold re-selection: eps -> eps + c with mu -> mu0[eps + c] changes f by  delta f(R) = 2 c (q_R - q0)/q0 + O(c^2);
     for linear K this is a shift of the curvature structure by delta R = 2 c/(q0 K'), for the affine gauge q - q0 = R/lam a
     renormalization of K by 2 c/(lam q0).
Run:  python3 adiabatic_relation.py  (from this directory; writes ../results/adiabatic_relation.json)
"""
import json
import os
import sympy as sp

out = {}

# ---------------------------------------------------------------- S1: general first-order formula in f(R) language
R, H, Rdot, Klab, Lam, rd, rr, f2, m = sp.symbols('R H Rdot K_lab Lambda rho_d rho_r f_2 m', real=True)
Fs = sp.Function('F')
# FRW 00 equation of f(R) gravity (S = int f/2):  3 f' H^2 = (R f' - f)/2 - 3 H d(f')/dt + rho
# with f = Klab R - 2 Lam + F(R):  solve for 3 Klab H^2 - rho  to first order in F.
fR = Klab * R - 2 * Lam + Fs(R)
f1 = sp.diff(fR, R)
f2expr = sp.diff(fR, R, 2)
lhs = 3 * f1 * H**2
rhs_nomatter = (R * f1 - fR) / 2 - 3 * H * f2expr * Rdot
rhoDE_eff = sp.simplify(rhs_nomatter - 3 * (f1 - Klab) * H**2)   # = 3 Klab H^2 - rho
general = sp.expand(rhoDE_eff)
out['S1_general_first_order'] = {
    'formula': 'rho_DE,eff = 3 K_lab H^2 - rho = ' + str(general),
    'reading': 'Lam + (R F\' - F)/2 - 3 F\' H^2 - 3 H F\'\' Rdot, exact in F for the 00 equation; first order when the '
               'background (R, H, Rdot) is taken from GR + Lam',
}
# quadratic F (unscreened q-mode): F = f2 R^2/2, with GR+Lam background including radiation (traceless) and dust
Fq = f2 * R**2 / 2
rho_dyn_q = (R * sp.diff(Fq, R) - Fq) / 2 - 3 * sp.diff(Fq, R) * H**2 - 3 * H * sp.diff(Fq, R, 2) * Rdot
H2_GR = (rd + rr + Lam) / (3 * Klab)
R_GR = (rd + 4 * Lam) / Klab
Rdot_over_H = -3 * rd / Klab                       # d/dt (rho_d/K) = -3 H rho_d/K
rho_dyn_sub = sp.expand(rho_dyn_q.subs(Rdot, Rdot_over_H * H).subs(H**2, H2_GR).subs(R, R_GR))
rho_dyn_sub = sp.expand(sp.simplify(rho_dyn_sub.subs(f2, Klab / (3 * m**2))))
target = ((sp.Rational(9, 4) * rd**2 + 2 * rd * rr - 4 * Lam * rr) / (3 * Klab * m**2))
out['S1_quadratic'] = {
    'rho_dyn': str(sp.factor(rho_dyn_sub)),
    'equals [(9/4) rho_d^2 + 2 rho_d rho_r - 4 Lam rho_r]/(3 K m^2)': bool(sp.simplify(rho_dyn_sub - sp.expand(target)) == 0),
    'Lam_dependence_without_radiation': str(sp.simplify(sp.diff(rho_dyn_sub.subs(rr, 0), Lam))),
}
# Omega_dyn today and w_dyn in the dust era
Om, H0 = sp.symbols('Omega_m H_0', positive=True)
rho_c = 3 * Klab * H0**2
Omega_dyn = sp.simplify((sp.Rational(3, 4) * (Om * rho_c)**2 / (Klab * m**2)) / rho_c)
out['S1_Omega_dyn_today'] = str(Omega_dyn)
a = sp.symbols('a', positive=True)
rho_dyn_a = sp.Rational(3, 4) * (Om * rho_c * a**-3)**2 / (Klab * m**2)
w_dyn = sp.simplify(-1 - sp.Rational(1, 3) * a * sp.diff(sp.log(rho_dyn_a), a))
out['S1_w_dyn_dust_era'] = str(w_dyn)

# ---------------------------------------------------------------- S2: the same from q-variables (Taylor data)
# Expand about q_R0 (the R = 0 point of this mu: eps'(q_R0) = mu); K0, K1, K2 and e2, e3 Taylor data there.
K0, K1, K2, e2, e3, e1 = sp.symbols('K0 K1 K2 e2 e3 e1', real=True)
dq = sp.symbols('dq')
# eps'(q) - mu = e2 dq + e3 dq^2/2 ;  K'(q) = K1 + K2 dq ;  q-equation  eps' - mu = K' R/2  -> dq(R)
dq_ser = sp.symbols('d1 d2')
d1, d2 = dq_ser
eqq = sp.expand((e2 * (d1 * R + d2 * R**2) + e3 * (d1 * R)**2 / 2) - (K1 + K2 * d1 * R) * R / 2)
sol = sp.solve([eqq.coeff(R, 1), eqq.coeff(R, 2)], [d1, d2], dict=True)[0]
dqR = sol[d1] * R + sol[d2] * R**2
KofR = K0 + K1 * dqR + K2 * dqR**2 / 2
f_series = sp.series(sp.integrate(sp.series(KofR, R, 0, 3).removeO(), R), R, 0, 4).removeO()   # f' = K(q_R)
f2_from_q = sp.simplify(sp.diff(f_series, R, 2).subs(R, 0))
out['S2_fpp_at_0_from_q'] = {'f\'\'(0)': str(f2_from_q), 'expected K1^2/(2 e2)': bool(sp.simplify(f2_from_q - K1**2 / (2 * e2)) == 0),
                             'm^2 = K0/(3 f\'\'(0)) = (2/3) K0 e2/K1^2': str(sp.simplify(K0 / (3 * f2_from_q)))}

# ---------------------------------------------------------------- S3: Lam_cal(mu) as a series in dmu (static root)
dmu = sp.symbols('dmu')
q0, E2, E3, k0, k1, k2 = sp.symbols('q0 E2 E3 k0 k1 k2', real=True)   # Taylor data at the Minkowski point (q0, mu0)
x = sp.symbols('x')       # q - q0
mu0 = sp.symbols('mu0')
epsT = E2 * x**2 / 2 + E3 * x**3 / 6 + mu0 * (q0 + x)          # eps(q0)-mu0 q0 = 0, eps'(q0) = mu0
KT = k0 + k1 * x + k2 * x**2 / 2
c1s, c2s = sp.symbols('c1s c2s')
xs = c1s * dmu + c2s * dmu**2
rhoV = epsT - (mu0 + dmu) * (q0 + x)
static = sp.diff(epsT, x) - (mu0 + dmu) - 2 * sp.diff(KT, x) * rhoV / KT
ser = sp.expand(sp.series(static.subs(x, xs), dmu, 0, 3).removeO())
s_sol = sp.solve([ser.coeff(dmu, 1), ser.coeff(dmu, 2)], [c1s, c2s], dict=True)[0]
rhoV_s = sp.expand(sp.series(rhoV.subs(x, xs).subs(s_sol), dmu, 0, 3).removeO())
out['S3_rhoV_static_series'] = {
    'rho_V,s(dmu)': str(sp.collect(sp.simplify(rhoV_s), dmu)),
    'leading': str(sp.simplify(rhoV_s.coeff(dmu, 1))),
    'note': 'rho_V,s = -q0 dmu + O(dmu^2): the calibrated part. With m^2 = (2/3) k0 E2/k1^2 the dmu^2 coefficient equals '
            '-1/(2 E2) + (4/3) q0^2/(k0 m^2) (checked below).',
    'dmu2_coefficient_check': bool(sp.simplify(rhoV_s.coeff(dmu, 2) - (-1 / (2 * E2) + sp.Rational(4, 3) * q0**2 * k1**2 / (k0 * sp.Rational(2, 3) * k0 * E2))) == 0),
}

# ---------------------------------------------------------------- S4: map and inverse
q, muS, lam = sp.symbols('q mu lambda', real=True)
fF = sp.Function('f')
Kq = sp.Function('K')(q)
epsq = sp.Function('epsilon')(q)
# envelope: f(R) = K(q_R) R - 2 eps(q_R) + 2 mu q_R with q_R stationary -> df/dR = K, df/dmu = 2 q_R
Lq = Kq * R - 2 * epsq + 2 * muS * q
stationarity = sp.diff(Lq, q)            # = 2*(K'R/2 - eps' + mu) = 0  <=> q-equation
out['S4_map'] = {
    'f(R) = K(q_R) R - 2 eps(q_R) + 2 mu q_R; stationarity in q is the q-equation': str(sp.simplify(stationarity / 2)),
    'df/dR (envelope)': 'K(q_R)', 'df/dmu (envelope)': '2 q_R',
    'f\'\'(R)': 'K\'^2 / (2 (eps\'\' - K\'\' R/2))  [from dq_R/dR = K\'/(2 D)]',
    'f(0)': '-2 (eps(q_0) - mu q_0) at eps\'(q_0) = mu  =  2 eps*(mu)  (Legendre transform)',
}
# canonical inverse: q - q0 = R/lam  (q tracks R affinely)
Rq = lam * (q - q0)
Kinv = sp.diff(fF(Rq), q) / lam                     # = f'(lam (q-q0))
eps_inv = (Kinv * Rq - fF(Rq)) / 2 + muS * q
q_eq_residual = sp.simplify(sp.diff(eps_inv, q) - sp.diff(Kinv, q) * Rq / 2 - muS)
f_back = sp.simplify(Kinv * Rq - 2 * eps_inv + 2 * muS * q - fF(Rq))
out['S4_inverse_affine_gauge'] = {
    'K(q)': 'f\'(lam (q - q0))', 'eps(q)': '[K(q) lam (q-q0) - f(lam (q-q0))]/2 + mu q',
    'q-equation residual at R = lam (q - q0) (should be 0)': str(q_eq_residual),
    'reconstructed f minus f (should be 0)': str(f_back),
    'reading': 'every C^2 function f(R) is the f of a q-theory; the realization is not unique (any reparametrization q(R) '
               'and the shift eps -> eps + d q, mu -> mu + d give the same f at this mu)',
}
# P_q in f-language
out['S4_P_q'] = {
    'statement': 'P_q (eps_inf\'(q0) = mu, eps_inf(q0) - mu q0 = 0, eps_inf\'\' > 0, K > 0)  <=>  f(0) = 0 with f\'(0) > 0 and '
                 'f\'\'(0) > 0, f being the full (final, renormalized) f of the selected mu. f(0) = 0 <=> R = 0 solves the '
                 'vacuum trace equation R f\' - 2 f = 0 (Minkowski is a vacuum solution). The final-state reading needs in '
                 'addition that the cosmological trajectory ends at R = 0; f(0) = 0 alone does not imply it.',
}

# ---------------------------------------------------------------- S5: threshold re-selection
cc = sp.symbols('c')
# first order: f_new(R) - f(R) = -2c + 2 dmu0 q_R + O(c^2) with dmu0 = c/q0  (equations track 5.1)
out['S5_threshold'] = {
    'delta f(R)': '2 c (q_R - q0)/q0 + O(c^2)',
    'linear K': 'K = K0 + K1 (q - q0): q_R - q0 = (f\'(R) - f\'(0))/K1  =>  f_new(R) = f(R + dR) - dR f\'(0) + O(c^2), '
                'dR = 2c/(q0 K1)  [first-order identity checked below; exact numerics in hu_sawicki.py B4]',
    'affine gauge': 'q - q0 = R/lam  =>  f_new(R) = f(R) + (2 c/(lam q0)) R + O(c^2), a renormalization of K only '
                    '[numerics in hu_sawicki.py B4]',
    'local statement': 'at a point of the curvature structure where the local slope is K\'(q_R), the structure moves by '
                       'dR = 2 c/(q0 K\'(q_R)); keeping a feature at R ~ R_dS ~ 4 rho_Lam/K requires q0 K\'/K >~ c/rho_Lam '
                       'across the feature',
}
# check the linear-K shift statement to first order with a generic f
fG = sp.Function('g')
K1s, eps_c = sp.symbols('K1s epsilon_c')
dR = 2 * eps_c / (q0 * K1s)
lhs_shift = sp.series(fG(R + dR) - dR * sp.Subs(sp.diff(fG(x), x), x, 0), eps_c, 0, 2).removeO()
rhs_shift = fG(R) + 2 * eps_c * (sp.diff(fG(R), R) - sp.Subs(sp.diff(fG(x), x), x, 0)) / (q0 * K1s)
out['S5_linearK_first_order_identity'] = bool(sp.simplify(sp.expand(lhs_shift.doit() - rhs_shift.doit())) == 0)

here = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(here, '..', 'results', 'adiabatic_relation.json')
with open(path, 'w') as fh:
    json.dump(out, fh, indent=1)
for k, v in out.items():
    print(k, ':', v)
