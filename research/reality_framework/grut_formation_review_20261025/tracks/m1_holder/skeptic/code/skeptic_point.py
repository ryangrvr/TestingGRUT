#!/usr/bin/env python3
"""m1_holder_skeptic (Claude, adversarial reviewer of track m1_holder; own code, independent of holder_point.py).

Recomputes parts 1-2 of track m1_holder with a DIFFERENT method:
  * (xi, alpha) algebra in 40-digit mpmath from the GeV pins of d4747cb (no float shortcuts), unit check in GeV;
  * stationary point from the three real amplitude equations at delta = -pi/2 with mpmath.findroot (Newton, 40 digits),
    independent of scipy BFGS/LM;
  * Hessian by mpmath numerical differentiation (mp.diff, 30 digits) of V written in the six real canonical fields,
    eigenvalues by mpmath.eigsy; independent of sympy and of Codex's A/B block formula;
  * rephasing-invariance check of alpha (numerical: rephase Y, Sigma, recompute arg(lock/mu) after making kappa, mu real);
  * the fixed-R lower-bound chain: each inequality re-derived symbolically (sympy, exact rationals where possible) and
    the final bound evaluated with mpmath;
  * NEW: size of the Wilson coefficient in terms of the declared cutoff Lambda = 1e15 GeV, and the dimension-6
    operator that any tree-level scalar-singlet completion of the dimension-5 holder induces, evaluated at the
    prepared point (estimate, [E]).
Units: fields and H_i = 1, R in H_i^2, V in H_i^4 unless a GeV check is stated.
"""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 40
HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
RES.mkdir(parents=True, exist_ok=True)

# pins [Q d4747cb CANDIDATE_M1.json]
H_GEV = mp.mpf('134284273491.30815')
MY_GEV = mp.mpf('4028528204739.244')
H_VEV_GEV = mp.mpf('5371370939652.326')
MU_GEV = mp.mpf('961402368973.181')
KAP = mp.mpf('0.0034586847211042523')
LX, LY, LS = mp.mpf('0.001'), mp.mpf(1), mp.mpf(1)
XIX = -mp.mpf(5) / 6
MPL_GEV = mp.mpf('2.435e18')
LAMBDA_GEV = mp.mpf('1e15')          # declared EFT cutoff [Q d4747cb]
M = MY_GEV / H_GEV
HV = H_VEV_GEV / H_GEV
MU = MU_GEV / H_GEV
ETA2 = (H_GEV / MPL_GEV) ** 2

CODEX = dict(xi=30.746728738599167, alpha=2.192144074179184,
             amps=[87.19556213287578, 5.75424796328559, 28.387181944361622], V0=-32480.585311406205,
             hess=[13.091408829453263, 14.589249978666821, 538.1466434275541, 1546.6597219731073, 3239.885251158698],
             bound=-285369.6842985777, B=0.5283522987365723, planck=1.0000000000480667)
TRACK = dict(xi=30.746728738599167, alpha=2.192144074179184,
             amps=[87.19556213287582, 5.754247963285605, 28.38718194436162], V0=-32480.58531140619,
             hess=[13.091408829452917, 14.589249978667548, 538.1466434275553, 1546.659721973108, 3239.8852511586965],
             loop_xi=0.001511259525661174, required_over_loop=20345.1)

out = {}

# ----------------------------------------------------------------------------------------------------------------
# 1. (xi, alpha) in GeV units, 40 digits
# ----------------------------------------------------------------------------------------------------------------
# V ⊃ (xi R/M_Y)(e^{i alpha} Sigma* Y^2 + h.c.) = 2 Re[(xi e^{i alpha}/M_Y) R Sigma* Y^2]; with mu Sigma* Y^2 + h.c. the
# coefficient of Sigma* Y^2 is u(R) = mu + (xi e^{i alpha}/M_Y) R.  Requirement [Q f8186fa]: u(12 H^2) = 10 i H.
lock_gev = (10j * H_GEV - MU_GEV) / (12 * H_GEV ** 2)           # GeV^-1
xi = abs(lock_gev) * MY_GEV
alpha = mp.arg(lock_gev)
lock_H = lock_gev * H_GEV                                        # dimensionless in H_i units
out['part1'] = dict(
    mu_over_H=float(MU), M_Y_over_H=float(M), h_over_H=float(HV),
    lock_GeV_inverse=[float(lock_gev.real), float(lock_gev.imag)], lock_over_H=[float(lock_H.real), float(lock_H.imag)],
    xi_own=float(xi), xi_codex=CODEX['xi'], xi_track=TRACK['xi'], xi_rel_diff_vs_codex=float(xi / CODEX['xi'] - 1),
    alpha_own=float(alpha), alpha_codex=CODEX['alpha'], alpha_diff_vs_codex=float(alpha - CODEX['alpha']),
    xi_40_digits=mp.nstr(xi, 30), alpha_40_digits=mp.nstr(alpha, 30),
    check_u_at_R12_over_H=[float((MU + lock_H * 12).real), float((MU + lock_H * 12).imag)],
    # brief's "mu = 7.2 H_i" vs pin
    mu_brief_7p2_vs_pin_rel=float(mp.mpf('7.2') / MU - 1),
    # if someone used mu = 7.2 exactly, what (xi, alpha) would follow
    xi_if_mu_7p2=float(abs((10j - mp.mpf('7.2')) / 12) * M), alpha_if_mu_7p2=float(mp.arg((10j - mp.mpf('7.2')) / 12)),
    # Wilson coefficient relative to the declared cutoff Lambda = 1e15 GeV: (xi/M_Y) = c_Lambda/Lambda
    c_Lambda=float(xi * LAMBDA_GEV / MY_GEV),
    holder_over_mu_at_R12=float(abs(lock_H) * 12 / MU),
    dim5_term_at_R3_over_mu=float(abs(lock_H) * 3 / MU),
    R12_over_MY2=float(12 / M ** 2), R12_over_Lambda2=float(12 * H_GEV ** 2 / LAMBDA_GEV ** 2),
    d_delta_eq_dlnR_at_R12=float(12 * mp.im(lock_H / (MU + 12 * lock_H))),
    delta_eq_at_R=[[float(R), float(mp.arg(MU + lock_H * R) + mp.pi - 2 * mp.pi * (mp.arg(MU + lock_H * R) + mp.pi > mp.pi))] for R in (12, 9, 6, 3, 1, 0)],
)

# one-loop size of xi (the track's estimate): c lambda_Y mu/(16 pi^2 M_Y).  Independent re-derivation of the parametric
# form: the dim-3 vertex mu Sigma* Y^2 dressed by a Y loop in curved space gives a term ~ mu lambda_Y R/(16 pi^2 M_Y^2)
# Sigma* Y^2 (heat-kernel a_1 coefficient (1/6 - xi_Y) R with the quartic vertex); so xi_loop ~ lambda_Y mu (1/6)/(16 pi^2 M_Y).
loop = LY * MU / (16 * mp.pi ** 2 * M)
out['part1']['xi_one_loop_estimate_c_eq_1'] = float(loop)
out['part1']['xi_one_loop_estimate_with_1_over_6'] = float(loop / 6)
out['part1']['required_over_loop_c1'] = float(xi / loop)
out['part1']['required_over_loop_with_1_over_6'] = float(6 * xi / loop)

# ----------------------------------------------------------------------------------------------------------------
# 2. stationary point at R = 12 by mpmath Newton on the three real amplitude equations (X real > 0, Y real < 0,
#    Sigma = -i s): own derivation of the equations from V(x, y, s) along delta = -pi/2, arg Y - 3 arg X = pi.
#    V(x,y,s) = lx x^4 + ly y^4 + M^2 y^2 + ls (s^2 - h^2/2)^2 - 2 kap y x^3 - 2 |u| s y^2 + xiX R x^2   (y>0 here = |Y|)
# ----------------------------------------------------------------------------------------------------------------
R12 = mp.mpf(12)
U12 = abs(MU + lock_H * R12)


def eqs(x, y, s):
    return [4 * LX * x ** 3 - 6 * KAP * y * x ** 2 + 2 * XIX * R12 * x,
            4 * LY * y ** 3 + 2 * M ** 2 * y - 2 * KAP * x ** 3 - 4 * U12 * s * y,
            4 * LS * (s ** 2 - HV ** 2 / 2) * s - 2 * U12 * y ** 2]


# own starting guess, different from the track's: x from the curvature/quartic balance ignoring kappa (x^2 = 5000)
x0 = mp.sqrt(-XIX * R12 / (2 * LX)); y0 = KAP * x0 ** 3 / M ** 2; s0 = HV / mp.sqrt(2)
sol = mp.findroot(lambda x, y, s: eqs(x, y, s), (x0, y0, s0), tol=mp.mpf(10) ** -35, maxsteps=100)
x, y, s = sol
amps_own = [x, y, s]


def V_complex(X, Y, S, R=R12, lock=lock_H):
    u = MU + lock * R
    return (LX * abs(X) ** 4 + LY * abs(Y) ** 4 + M ** 2 * abs(Y) ** 2 + LS * (abs(S) ** 2 - HV ** 2 / 2) ** 2
            + 2 * mp.re(KAP * mp.conj(Y) * X ** 3 + u * mp.conj(S) * Y ** 2) + XIX * R * abs(X) ** 2)


def V_real(r, R=R12, lock=lock_H):
    X = (r[0] + 1j * r[1]) / mp.sqrt(2); Y = (r[2] + 1j * r[3]) / mp.sqrt(2); S = (r[4] + 1j * r[5]) / mp.sqrt(2)
    return V_complex(X, Y, S, R, lock)


Z = [x, -y, -1j * s]
V0 = V_complex(*Z)
r0 = [mp.sqrt(2) * x, mp.mpf(0), -mp.sqrt(2) * y, mp.mpf(0), mp.mpf(0), -mp.sqrt(2) * s]
# gradient by mp.diff (check that the complex-field stationary point is stationary in all six real directions)
grad = [mp.diff(lambda t, i=i: V_real([r0[j] + (t if j == i else 0) for j in range(6)]), 0) for i in range(6)]
# Hessian by mp.diff, second derivatives of V_real (30 digits)
Hm = mp.matrix(6, 6)
for i in range(6):
    for j in range(i, 6):
        if i == j:
            d = mp.diff(lambda t, i=i: V_real([r0[k] + (t if k == i else 0) for k in range(6)]), 0, 2)
        else:
            d = mp.diff(lambda a, b, i=i, j=j: V_real([r0[k] + (a if k == i else 0) + (b if k == j else 0) for k in range(6)]), (0, 0), (1, 1))
        Hm[i, j] = Hm[j, i] = d
ev, Q = mp.eigsy(Hm)
ev = sorted([ev[i] for i in range(6)])
phys_own = [float(e) for e in ev[1:]]
T = mp.matrix(6, 6)
for i, q in enumerate([1, 3, 6]):
    T[2 * i, 2 * i + 1] = -q; T[2 * i + 1, 2 * i] = q
gvec = T * mp.matrix(r0); gvec = gvec / mp.norm(gvec)
gauge_res = float(mp.norm(Hm * gvec))
VR = XIX * x ** 2 + 2 * mp.re(lock_H * mp.conj(Z[2]) * Z[1] ** 2)
out['part2_point'] = dict(
    amplitudes_own=[float(v) for v in amps_own], amplitudes_codex=CODEX['amps'], amplitudes_track=TRACK['amps'],
    amplitudes_30_digits=[mp.nstr(v, 25) for v in amps_own],
    rel_diff_vs_codex=[float(a / b - 1) for a, b in zip(amps_own, CODEX['amps'])],
    equation_residuals=[float(v) for v in eqs(x, y, s)], six_real_gradient_mp=[float(g) for g in grad],
    V_own=float(V0), V_codex=CODEX['V0'], V_track=TRACK['V0'], V_rel_diff_vs_codex=float(V0 / CODEX['V0'] - 1),
    hessian_eigenvalues_mpdiff=[float(e) for e in ev], physical_own=phys_own, physical_codex=CODEX['hess'], physical_track=TRACK['hess'],
    physical_rel_diff_vs_codex=[a / b - 1 for a, b in zip(phys_own, CODEX['hess'])],
    gauge_null_residual=gauge_res, all_physical_positive=all(e > 0 for e in phys_own),
    dV_dR=float(VR), effective_Planck_ratio=float(1 - 2 * ETA2 * VR), effective_Planck_codex=CODEX['planck'],
    X_over_MY=float(x / M), Y_over_MY=float(y / M), S_over_MY=float(s / M),
    X_over_Lambda=float(x * H_GEV / LAMBDA_GEV),
    kappa2_X6_over_MY2_vs_lambdaX_X4=float(KAP ** 2 * x ** 2 / (M ** 2 * LX)),
    holder_term_at_point_over_H4=float(2 * mp.re(lock_H * R12 * mp.conj(Z[2]) * Z[1] ** 2)),
    mu_term_at_point_over_H4=float(2 * mp.re(MU * mp.conj(Z[2]) * Z[1] ** 2)),
    kappa_term_at_point_over_H4=float(2 * mp.re(KAP * mp.conj(Z[1]) * Z[0] ** 3)),
    MY2_Y2_over_H4=float(M ** 2 * y ** 2), lambdaY_Y4_over_H4=float(LY * y ** 4), lambdaX_X4_over_H4=float(LX * x ** 4), xiR_X2_over_H4=float(XIX * 12 * x ** 2),
)

# rephasing invariance of alpha: random U(1)^3 rephasing of (X, Y, Sigma); make kappa and mu real positive by field
# rephasing; the remaining lock phase equals alpha (mod 2 pi) -> alpha is rephasing-invariant
rng = np.random.default_rng(7)
diffs = []
for _ in range(5):
    a, b, c = rng.uniform(-np.pi, np.pi, 3)
    kap_p = KAP * mp.exp(1j * (3 * a - b)); mu_p = MU * mp.exp(1j * (2 * b - c)); lock_p = lock_H * mp.exp(1j * (2 * b - c))
    # undo: choose b' , c' to make kap_p, mu_p real positive (a' = 0): b' = arg(kap_p), c' = 2 b' - arg(mu_p)
    bp = mp.arg(kap_p); cp = 2 * bp - mp.arg(mu_p)
    lock_pp = lock_p * mp.exp(1j * (2 * bp - cp))   # note: after rephasing with (0, bp, cp): factor e^{i(2 bp - cp)} ... sign
    # general rule: under Y -> e^{i b} Y, S -> e^{i c} S, the coefficient of S* Y^2 picks e^{i(2b - c)}; to make mu real we need
    # 2b - c = -arg(mu_p); kappa (Y* X^3) picks e^{-i b} (with a = 0): b = arg(kap_p).
    b2 = mp.arg(kap_p); c2 = 2 * b2 + mp.arg(mu_p)
    mu_fixed = mu_p * mp.exp(1j * (2 * b2 - c2)); kap_fixed = kap_p * mp.exp(-1j * b2); lock_fixed = lock_p * mp.exp(1j * (2 * b2 - c2))
    assert abs(mp.im(mu_fixed)) < 1e-30 and abs(mp.im(kap_fixed)) < 1e-30
    diffs.append(float((mp.arg(lock_fixed) - alpha + mp.pi) % (2 * mp.pi) - mp.pi))
out['part2_point']['alpha_rephasing_invariance_max_abs_diff'] = max(abs(d) for d in diffs)

# ----------------------------------------------------------------------------------------------------------------
# 3. the fixed-R lower bound: symbolic re-derivation of each step
# ----------------------------------------------------------------------------------------------------------------
xs, ys, ss, c, k_, u_, A_, ax_, lam, hh = sp.symbols('x y s c kappa u A a_x lambda h', positive=True)
# step 1: min_x [c x^4 - 2 kappa y x^3]
f1 = c * xs ** 4 - 2 * k_ * ys * xs ** 3
crit = sp.solve(sp.diff(f1, xs), xs)
vals = [sp.simplify(f1.subs(xs, cr)) for cr in crit if cr != 0]
step1 = sp.simplify(vals[0])                               # expect -27 kappa^4 y^4/(16 c^3)
# step 2: min over y^2 of B y^4 - 2 u s y^2  -> -(u s)^2/B
Bsym = sp.symbols('B', positive=True); w = sp.symbols('w', nonnegative=True)
f2 = Bsym * w ** 2 - 2 * u_ * ss * w
step2 = sp.simplify(f2.subs(w, sp.solve(sp.diff(f2, w), w)[0]))
# step 3: min over x^2 of A x^4 - a_x x^2 -> -a_x^2/(4A)
f3 = A_ * w ** 2 - ax_ * w
step3 = sp.simplify(f3.subs(w, sp.solve(sp.diff(f3, w), w)[0]))
# step 4: min over s^2 of lambda (s^2 - h^2/2)^2 - (u^2/B) s^2
cc = sp.symbols('C', positive=True)
f4 = lam * (w - hh ** 2 / 2) ** 2 - cc * w
step4 = sp.simplify(f4.subs(w, sp.solve(sp.diff(f4, w), w)[0]))
codex_form = sp.simplify(lam * hh ** 4 / 4 - (lam * hh ** 2 + cc) ** 2 / (4 * lam))
A = mp.mpf('0.2') * LX; cnum = mp.mpf('0.8') * LX
B = LY - 27 * KAP ** 4 / (16 * cnum ** 3)
ax = max(-XIX * R12, 0)
bound = LS * HV ** 4 / 4 - ax ** 2 / (4 * A) - (LS * HV ** 2 + U12 ** 2 / B) ** 2 / (4 * LS)
# optimise the reserve fraction (Codex chose 0.2): the bound as a function of the reserve
def bound_of_reserve(res):
    res = mp.mpf(res); A = res * LX; cn = (1 - res) * LX; Bv = LY - 27 * KAP ** 4 / (16 * cn ** 3)
    if Bv <= 0:
        return None
    return float(LS * HV ** 4 / 4 - ax ** 2 / (4 * A) - (LS * HV ** 2 + U12 ** 2 / Bv) ** 2 / (4 * LS))
scan = {str(r): bound_of_reserve(r) for r in ('0.05', '0.1', '0.2', '0.3', '0.4', '0.5')}
# minimal B > 0 condition: 27 kappa^4/(16 c^3) < lambda_Y  ->  c > (27 kappa^4/(16 lambda_Y))^(1/3)
c_min = (27 * KAP ** 4 / (16 * LY)) ** (mp.mpf(1) / 3)
out['part2_bound'] = dict(
    step1_min_x=str(step1), step1_matches=bool(sp.simplify(step1 + 27 * k_ ** 4 * ys ** 4 / (16 * c ** 3)) == 0),
    step2_min_y2=str(step2), step3_min_x2=str(step3), step4_min_s2=str(step4), step4_equals_codex_form=bool(sp.simplify(step4 - codex_form) == 0),
    B_own=float(B), B_codex=CODEX['B'], bound_own=float(bound), bound_codex=CODEX['bound'], bound_rel_diff=float(bound / CODEX['bound'] - 1),
    bound_vs_reserve=scan, c_min_for_B_positive=float(c_min), lambda_X_margin_over_c_min=float(LX / c_min),
    comment='each inequality re-derived symbolically; Codex chain exact; of the tested reserves 0.2 gives the tightest bound (0.1: -375217, 0.3: -382175; '
            '0.4 and 0.5 make B < 0); the proof needs (1-reserve) lambda_X > c_min = 6.2e-4, i.e. a margin of only 1.6 in lambda_X (1.42 in kappa)',
)

# ----------------------------------------------------------------------------------------------------------------
# 4. NEW [E]: what a tree-level scalar-singlet completion of the dim-5 holder induces.
#    L ⊃ -1/2 M_S^2 S^2 + S J,  J = g1 M_S R + (g2 Sigma* Y^2 + h.c.)  ->  V_eff = -J^2/(2 M_S^2):
#    cross term  -(g1 g2/M_S) R (Sigma* Y^2 + h.c.)  => xi e^{i alpha}/M_Y = -g1 g2/M_S  => |g1 g2| = xi M_S/M_Y
#    square term -(g2 Sigma* Y^2 + h.c.)^2/(2 M_S^2) = -2 |g2|^2 |Sigma|^2 |Y|^4 cos^2(.)/M_S^2  (dimension 6)
#    also -g1^2 R^2/2 (an R^2 term, coefficient g1^2/2 in units where the operator is M_S^2 S^2 ... dimensionless)
# ----------------------------------------------------------------------------------------------------------------
def completion(MS_over_MY, g1):
    MS = MS_over_MY * M
    g2 = xi * MS / (M * g1)
    # worst-case (cos^2 = 1) induced |Y|^4 coefficient at the prepared point relative to lambda_Y
    d_lamY = 2 * g2 ** 2 * s ** 2 / MS ** 2
    return dict(M_S_over_MY=float(MS_over_MY), g1=float(g1), g2=float(g2), induced_minus_delta_lambdaY_over_lambdaY=float(d_lamY / LY),
                B_after=float(B - d_lamY), induced_term_at_point_over_H4=float(-2 * g2 ** 2 * s ** 2 * y ** 4 / MS ** 2),
                dS_tadpole_over_H=float(g1 * 12 / MS))
comp = [completion(ms, g1) for ms in (1.0, 5.0, 30.0 / 30.0 * (LAMBDA_GEV / MY_GEV)) for g1 in (1.0, 4 * np.pi, 60.0)]
# minimal g1 so that B stays positive (perturbative-g2 requirement): 2 g2^2 s^2/M_S^2 < B with g2 = xi M_S/(M g1) -> g1^2 > 2 xi^2 s^2/(M^2 B)
g1_min = mp.sqrt(2 * xi ** 2 * s ** 2 / (M ** 2 * B))
out['part2_completion_estimate'] = dict(
    rows=comp, g1_min_for_B_positive_any_MS=float(g1_min),
    g1g2_needed_MS_eq_MY=float(xi), g1g2_needed_MS_eq_Lambda=float(xi * LAMBDA_GEV / MY_GEV),
    comment='[E] tree-level singlet exchange: the induced -(2 Re g2 Sigma* Y^2)^2/(2 M_S^2) lowers the Y quartic at the prepared point by 2 g2^2 s^2/M_S^2 '
            '(independent of M_S once g1 g2 = xi M_S/M_Y is imposed, for fixed g1): B > 0 needs g1 > g1_min ~ 60; with g1 <= 4 pi the fixed-R bound proof fails.',
)

# ----------------------------------------------------------------------------------------------------------------
# 5. a few extra consistency items
# ----------------------------------------------------------------------------------------------------------------
# (a) m_A at the prepared point for e = 1e-6 (track quotes 2.7e-4 H_i): m_A^2 = 2 e^2 sum q^2 |z|^2
e = mp.mpf('1e-6')
mA = mp.sqrt(2 * e ** 2 * (1 * x ** 2 + 9 * y ** 2 + 36 * s ** 2))
# (b) the H_i refit law
ratio_track = mp.mpf('1.0180645761996148'); ratio_codex = mp.mpf('1.0180645288398746')
out['extras'] = dict(
    m_A_over_H=float(mA),
    H_ratio_from_codex_ratio=float(ratio_codex ** (-mp.mpf(2) / 3)), H_new_GeV_from_codex_ratio=float(H_GEV * ratio_codex ** (-mp.mpf(2) / 3)),
    H_ratio_from_track_ratio=float(ratio_track ** (-mp.mpf(2) / 3)),
    H_ratio_fully_cooled_codex=float((ratio_codex * (1 + mp.mpf('0.008409936530420273'))) ** (-mp.mpf(2) / 3)),
    # sensitivity: a window change of the ratio by 1.5e-3 moves H_i by
    dH_over_H_per_1e3_ratio=float(-mp.mpf(2) / 3 * mp.mpf('1e-3')),
)
json.dump(out, open(RES / 'skeptic_point.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
