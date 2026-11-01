"""m1_vacuum SKEPTIC, part 1: exact/high-precision algebra (sympy + mpmath), written to challenge
tracks/m1_vacuum items 1, 3, 4 with a different toolchain than the track's float-numpy code.
Independent code (Claude skeptic, 1 Nov 2026). Nothing from Codex's code is used.

Checks: (A) V' factorisation, discriminant, thresholds, positivity bracket -- as exact rationals;
        (B) bulk quartic derived symbolically, root count by Sturm (not Descartes), 30-digit root;
        (C) dispersion cubic derived from the 6x6 first-order linear system (not from the 3x3 determinant);
        (D) long-wave coefficient vs thermodynamic sound speed -- symbolic identity along the family;
        (E) spinodal density by solving a D = 2 S^2 at 30 digits;
        (F) hidden tree bound: exact numbers AND a brute-force minimisation of the actual rp.py potential
            (lambda_X|X|^4 + M^2|Y|^2 + lambda_Y|Y|^4 + lambda_S(|S|^2-h^2/2)^2 + 2Re[kappa Y* X^3 + mu S* Y^2])
            at fixed |Sigma| to test V_min(s) >= s^2 * margin, and a search for any V < 0;
        (G) wall numbers: V(-v), V(0), kink tension (exact), roll time by mpmath quad.
"""
import json, sys
import sympy as sp
import mpmath as mp
import numpy as np
from scipy.optimize import minimize

mp.mp.dps = 30
OUT = sys.argv[1] if len(sys.argv) > 1 else 'sk_algebra.json'
res = {}

lam, eps = sp.Rational(3, 100), sp.Rational(2, 1000)
a, b, c = sp.Rational(1, 10), sp.Integer(3), sp.Rational(1, 5)      # lambda_X/g^2, lambda_phi/g^2, eps/g^2
u, S, F2, z, k2 = sp.symbols('u S F2 z k2', real=True)

# ---------------------------------------------------------------- (A) one-vacuum algebra, exact
Vn = lam / 4 * (u**2 - 1)**2 - eps * (u - 1)**3            # V/v^4
Vp = sp.diff(Vn, u)                                          # V'/v^3
quad = sp.cancel(Vp / (u - 1))
disc = sp.discriminant(quad, u)
eps_s = sp.symbols('eps_s', positive=True)
quad_s = lam * u**2 + (lam - 3 * eps_s) * u + 3 * eps_s
thr = sp.solve(sp.discriminant(quad_s, u), eps_s)
bracket = sp.cancel(Vn / (u - 1)**2)
bdisc = sp.discriminant(bracket, u)
Vppp = sp.diff(Vn, u, 3)
res['A_one_vacuum'] = dict(
    Vprime_over_v3_factor=str(sp.factor(Vp)),
    quadratic=str(sp.expand(quad)),
    discriminant_exact=str(disc), discriminant_float=float(disc), codex=-0.000144,
    thresholds_exact=[str(t) for t in thr], thresholds_float=[float(t) for t in thr],
    codex_threshold=0.001715728752538096,
    bracket=str(sp.expand(bracket)), bracket_discriminant_exact=str(bdisc), bracket_discriminant_float=float(bdisc),
    bracket_leading_coefficient=str(sp.Poly(bracket, u).LC()),
    V_minus_v_over_v4=str(Vn.subs(u, -1)), V_zero_over_v4=str(Vn.subs(u, 0)),
    Vppp_at_v_over_v=str(Vppp.subs(u, 1)), Vppp_symmetric=str(6 * lam), cubic_reduction=str(sp.simplify(1 - Vppp.subs(u, 1) / (6 * lam))),
    inflection_u=[float(r) for r in sp.solve(sp.diff(Vn, u, 2), u)],
    # a second stationary point appears for eps below the lower threshold: where is it at eps = 0.001 (Codex grid point)?
    stationary_points_eps_0p001=[float(r) for r in sp.Poly(quad_s.subs(eps_s, sp.Rational(1, 1000)), u).nroots()],
    V_at_metastable_min_eps_0p001_over_v4=float(min(Vn.subs(eps, sp.Rational(1, 1000)).subs(u, r) for r in sp.Poly(quad_s.subs(eps_s, sp.Rational(1, 1000)), u).nroots())),
)

# ---------------------------------------------------------------- (B) bulk quartic, exact, root count by Sturm
Uphi = b / 4 * (S**2 - 1)**2 - c * (S - 1)**3
Uphi_p = sp.diff(Uphi, S)
# Omega^2(S) = S^2 + 2 sqrt(a Uphi); stationarity 2 S sqrt(Uphi) = - sqrt(a) Uphi'; square and divide by (S-1)^2
sq = sp.expand(4 * S**2 * Uphi - a * Uphi_p**2)
quart = sp.cancel(sq / (S - 1)**2)
Pq = sp.Poly(sp.expand(quart), S)
coeffs = [sp.nsimplify(cf) for cf in Pq.all_coeffs()]
n_pos_roots_sturm = Pq.count_roots(0, sp.oo)
n_roots_01 = Pq.count_roots(0, 1)
Sb = mp.findroot(lambda s: sum(mp.mpf(str(cf)) * s**(4 - i) for i, cf in enumerate(coeffs)), mp.mpf('0.15'))
Uphi_f = lambda s: mp.mpf(3) / 4 * (s * s - 1)**2 - mp.mpf('0.2') * (s - 1)**3
Uphi_pf = lambda s: 3 * (s * s - 1) * s - 3 * mp.mpf('0.2') * (s - 1)**2
F2b = 2 * mp.sqrt(Uphi_f(Sb) / mp.mpf('0.1'))
Omb = mp.sqrt(Sb**2 + mp.mpf('0.1') * F2b)
# unsquared equation check (sign): 2 S sqrt(Uphi) + sqrt(a) Uphi' = 0 at the root
unsq = 2 * Sb * mp.sqrt(Uphi_f(Sb)) + mp.sqrt(mp.mpf('0.1')) * Uphi_pf(Sb)
# also check the other real root (-0.0747) is spurious (sign of unsquared equation)
Sneg = mp.findroot(lambda s: sum(mp.mpf(str(cf)) * s**(4 - i) for i, cf in enumerate(coeffs)), mp.mpf('-0.07'))
unsq_neg = 2 * Sneg * mp.sqrt(Uphi_f(Sneg)) + mp.sqrt(mp.mpf('0.1')) * Uphi_pf(Sneg)
res['B_bulk'] = dict(
    quartic_coefficients_exact=[str(cf) for cf in coeffs], quartic_float=[float(cf) for cf in coeffs],
    codex_quartic=[2.1, 3.76, 2.864, -0.288, -0.036],
    positive_real_roots_by_Sturm=int(n_pos_roots_sturm), real_roots_in_0_1_by_Sturm=int(n_roots_01),
    S_bulk_30digits=mp.nstr(Sb, 25), F_bulk=mp.nstr(mp.sqrt(F2b), 25), Omega_bulk=mp.nstr(Omb, 25),
    n_bulk_over_mX3=mp.nstr(Omb * F2b * 100, 25),
    unsquared_equation_residual_at_root=mp.nstr(unsq, 5),
    other_real_root=mp.nstr(Sneg, 12), unsquared_residual_at_other_root=mp.nstr(unsq_neg, 8),
    codex=dict(S=0.1511721078446304, F=2.406476756777611, Omega=0.7758647074595397, n=449.3133678845478),
    Uphi_at_bulk_mX4=mp.nstr(Uphi_f(Sb) * 100, 12),
    symmetric_Omega_bulk=mp.nstr(mp.sqrt(2 * mp.sqrt(mp.mpf('0.1') * 3 / 4)), 20),
)

# ---------------------------------------------------------------- (C) dispersion cubic from the 6x6 first-order system
Fs, Ss, Om, D, nu, kk = sp.symbols('F S Omega D nu k', real=True)
# state (psi1, psi2, s, psi1', psi2', s'), d/dt = nu.  Equations (my derivation from the EOM; see SKEPTIC.md):
#  psi1'' = -(k^2 + 2aF^2) psi1 - 2 Om psi2' - 2 F S s
#  psi2'' = -k^2 psi2 + 2 Om psi1'
#  s''    = -(k^2 + D) s - 2 F S psi1
M = sp.zeros(6, 6)
M[0, 3] = 1; M[1, 4] = 1; M[2, 5] = 1
M[3, 0] = -(kk**2 + 2 * a * Fs**2); M[3, 4] = -2 * Om; M[3, 2] = -2 * Fs * Ss
M[4, 1] = -kk**2; M[4, 3] = 2 * Om
M[5, 2] = -(kk**2 + D); M[5, 0] = -2 * Fs * Ss
cp_expr = sp.expand((M - nu * sp.eye(6)).det())          # det(M - nu) = 0 <=> nu is a growth exponent
# the char poly must be even in nu; compare with Codex's cubic in z = nu^2
odd_part = sp.simplify((cp_expr - cp_expr.subs(nu, -nu)) / 2)
cubic_codex_nu = sp.expand((((nu**2 + kk**2 + 2 * a * Fs**2) * (nu**2 + kk**2) + 4 * Om**2 * nu**2) * (nu**2 + kk**2 + D) - 4 * Fs**2 * Ss**2 * (nu**2 + kk**2)))
res['C_dispersion_cubic'] = dict(odd_part_of_charpoly=str(odd_part),
                                 mine_minus_codex=str(sp.simplify(cp_expr - cubic_codex_nu)),
                                 mine_charpoly_6x6=str(cp_expr))

# ---------------------------------------------------------------- (D) long-wave coefficient = -c_s^2 along the family (symbolic)
# family: Uphi'(S) + F2 S = 0, Omega^2 = S^2 + a F2, n = Omega F2.  Parametrise by S: F2(S) = -Uphi'(S)/S.
F2S = sp.simplify(-Uphi_p / S)
Om2S = sp.simplify(S**2 + a * F2S)
nS = sp.sqrt(Om2S) * F2S
dOm_dS = sp.diff(sp.sqrt(Om2S), S)
dn_dS = sp.diff(nS, S)
cs2 = sp.simplify(nS * (dOm_dS / dn_dS) / sp.sqrt(Om2S))         # n dOmega/dn / Omega
DS = sp.diff(Uphi, S, 2) + F2S
sigma_lw = sp.simplify(F2S * (2 * S**2 - a * DS) / (DS * (a * F2S + 2 * Om2S) - 2 * F2S * S**2))
identity = sp.simplify(sigma_lw + cs2)
res['D_longwave_identity'] = dict(sigma_plus_cs2_simplified=str(identity),
                                  numeric_check_S_0p8=float((sigma_lw + cs2).subs(S, sp.Rational(8, 10))),
                                  numeric_check_S_0p3=float((sigma_lw + cs2).subs(S, sp.Rational(3, 10))))

# ---------------------------------------------------------------- (E) spinodal: a D = 2 S^2 along the family (30 digits)
def family(s):
    s = mp.mpf(s)
    F2 = -Uphi_pf(s) / s
    Om = mp.sqrt(s * s + mp.mpf('0.1') * F2)
    Dv = 3 * (3 * s * s - 1) - 6 * mp.mpf('0.2') * (s - 1) + F2
    return F2, Om, Dv
S_sp = mp.findroot(lambda s: mp.mpf('0.1') * family(s)[2] - 2 * s * s, mp.mpf('0.336'))
F2sp, Omsp, Dsp = family(S_sp)
res['E_spinodal'] = dict(S=mp.nstr(S_sp, 15), F=mp.nstr(mp.sqrt(F2sp), 15), Omega=mp.nstr(Omsp, 15),
                         n_over_mX3=mp.nstr(Omsp * F2sp * 100, 15), track=dict(F=1.8566256, Omega=0.6765785, n=233.2206))
# Omega minimum along the family is at the same point (dOmega/dn = 0 <=> dOmega/dS = 0)
S_Ommin = mp.findroot(lambda s: mp.diff(lambda t: family(t)[1], s), mp.mpf('0.336'))
res['E_spinodal']['S_at_Omega_min'] = mp.nstr(S_Ommin, 15)

# ---------------------------------------------------------------- (F) hidden tree bound: exact numbers + brute force on the actual potential
kap = mp.mpf('0.0034586847211042523'); lX = mp.mpf('0.001'); lY = mp.mpf(1); lS = mp.mpf(1)
MY = mp.mpf(30); hh = mp.mpf(40); mu = mp.mpf('7.159456159513794')
lY_eff = lY - 27 * kap**4 / (16 * lX**3)
s_c = MY**2 / (2 * mu); r2 = hh**2 / (2 * s_c**2); margin = (1 - r2)**2 * s_c**2 - mu**2 / lY_eff
# ratios from CANDIDATE_M1.json heavy block
Hg = 134284273491.30815
ratios = dict(MY_over_H=4028528204739.244 / Hg, h_over_H=5371370939652.326 / Hg, mu_over_H=961402368973.181 / Hg)
res['F_hidden_bound'] = dict(lamY_eff=mp.nstr(lY_eff, 17), codex=0.758516376953125, s_c=mp.nstr(s_c, 15), r2=mp.nstr(r2, 10),
                             margin=mp.nstr(margin, 17), codex_margin=2445.0408755072076, ratios_from_CANDIDATE_M1=ratios)
# brute force: V(x, y, s; phases) from rp.py's form. Fix |Sigma| = s, minimise over x>=0, y>=0 and the two relative phases.
def Vfull(p, s, phase_s=0.0):
    x, y, thx, thy = p
    X = x * np.exp(1j * thx); Y = y * np.exp(1j * thy); Sg = s * np.exp(1j * phase_s)
    return (float(lX) * x**4 + float(MY)**2 * y**2 + float(lY) * y**4 + float(lS) * (s * s - float(hh)**2 / 2)**2
            + 2 * np.real(float(kap) * np.conj(Y) * X**3 + float(mu) * np.conj(Sg) * Y * Y))
rng = np.random.default_rng(7)
bf = []
for s in [0.0, 20.0, float(s_c) * 0.999, float(s_c), 70.0, 100.0, 200.0, 500.0]:
    best = np.inf
    for _ in range(40):
        p0 = [rng.uniform(0, 400), rng.uniform(0, 100), rng.uniform(-np.pi, np.pi), rng.uniform(-np.pi, np.pi)]
        r = minimize(Vfull, p0, args=(s,), method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-10, 'maxiter': 40000})
        if r.fun < best:
            best, bx = r.fun, r.x
    bound = float(s * s * margin) if s >= float(s_c) else 0.0
    bf.append(dict(s=s, V_min_bruteforce=float(best), analytic_lower_bound=bound, bound_holds=bool(best >= bound - 1e-6 * max(1, abs(bound))),
                   argmin_x_y=[float(abs(bx[0])), float(abs(bx[1]))]))
res['F_hidden_bound']['brute_force_at_fixed_s'] = bf
# global random search for V < 0 (any phases), including large fields
samp = rng.uniform(0, 1, (200000, 3)) * np.array([2000, 400, 400])
ph = rng.uniform(-np.pi, np.pi, (200000, 3))
X = samp[:, 0] * np.exp(1j * ph[:, 0]); Y = samp[:, 1] * np.exp(1j * ph[:, 1]); Sg = samp[:, 2] * np.exp(1j * ph[:, 2])
Vs = (float(lX) * abs(X)**4 + float(MY)**2 * abs(Y)**2 + float(lY) * abs(Y)**4 + float(lS) * (abs(Sg)**2 - float(hh)**2 / 2)**2
      + 2 * np.real(float(kap) * np.conj(Y) * X**3 + float(mu) * np.conj(Sg) * Y * Y))
res['F_hidden_bound']['global_random_min_V'] = float(Vs.min())
res['F_hidden_bound']['note'] = ('the bound uses the whole light quartic lambda_X|X|^4 to absorb the kappa term; the remaining light terms '
                                 'g^2 phi^2|X|^2 >= 0 and V_phi >= 0 are then nonnegative. No charge-six or other UV operator is included (Codex says so).')

# ---------------------------------------------------------------- (G) walls: exact V(-v), V(0), kink tension, roll time (mpmath quad)
v = mp.mpf(10)   # v / m_X
Vphys = lambda ph: mp.mpf('0.03') / 4 * (ph**2 - v**2)**2 - mp.mpf('0.002') * v * (ph - v)**3
dV = Vphys(-v)
t_roll = mp.quad(lambda ph: 1 / mp.sqrt(2 * (dV - Vphys(ph))), [-v, 0, v])
sigma0 = mp.quad(lambda ph: mp.sqrt(2 * mp.mpf('0.03') / 4) * (v**2 - ph**2), [-v, v])
res['G_walls'] = dict(V_minus_v_mX4=mp.nstr(dV, 12), V_zero_mX4=mp.nstr(Vphys(0), 12), sigma0_symmetric_kink_mX3=mp.nstr(sigma0, 15),
                      sigma0_closed_form=mp.nstr(4 * mp.sqrt(mp.mpf('0.03') / 2) * v**3 / 3, 15),
                      roll_time_mX_inverse=mp.nstr(t_roll, 12), track=dict(roll=1.905866578012243, sigma0=163.29931618554522),
                      m_phi_over_mX=mp.nstr(mp.sqrt(mp.mpf('0.06')) * v, 12),
                      m_phi_MeV=mp.nstr(mp.sqrt(mp.mpf('0.06')) * mp.mpf('0.924089637375362') * 1000, 10),
                      m_X_MeV=mp.nstr(mp.mpf('0.0924089637375362') * 1000, 10))

json.dump(res, open(OUT, 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
