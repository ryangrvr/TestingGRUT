#!/usr/bin/env python3
"""Checks for track grut_native (packet grut_vacuum_principle_20261021).

Run:  python3 -I code/grut_native_checks.py results/grut_native_checks.json
Deterministic (seed 20261021). numpy, scipy, sympy only. About 10 s.

Sections (labels match NOTES.md):
  A  HT minisuperspace: the compensated vacuum shift is a boundary term; generator form.     [R1]
  B  Page-Wootters / access toy: matter shift = clock phase; which state rules are invariant. [R1, K6]
  C  Shift-invariant state rules f(H)/Tr f(H): only exponentials (Cauchy).                  [R1a, K5]
  D  Gibbs/KMS state of the whole w.r.t. the 4-volume clock; self-consistent KMS scaling.    [R1a, K2]
  E  K1: exponent p in rho/M^4 = (V M^4)^(-p) across conventions (retrospective).          [R1b]
  F  One global variable: nonlinear sigma does not decouple rho_vac; linear sigma = fixed V0. [R2a]
  G  Fixed finite 4-volume vs constant Lambda_eff > 0 and present acceleration (K3);
     closed-dust V4(lambda) map (conjugate-pair relocation).                                  [R2a]
  H  Relational subtraction: KP (1/4<T>) vs Tseytlin (S/V) residuals.                        [R2b]
  I  Energy-budget closure (K4): constraint identically satisfied; matter-only version sign.  [R3]
"""
import json
import sys
import math
import platform

import numpy as np
import scipy
import sympy as sp
from scipy import integrate, optimize

SEED = 20261021
rng = np.random.default_rng(SEED)
OUT = {"seed": SEED, "versions": {"python": platform.python_version(), "numpy": np.__version__,
                                   "scipy": scipy.__version__, "sympy": sp.__version__},
       "checks": {}, "numbers": {}}


def check(name, ok, detail=None):
    OUT["checks"][name] = {"pass": bool(ok), "detail": detail}
    print(("PASS " if ok else "FAIL ") + name + ("" if detail is None else f"  {detail}"))


def num(name, value):
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    OUT["numbers"][name] = value
    return value


# ---------------------------------------------------------------------------------------------
# A. Henneaux-Teitelboim FRW minisuperspace (unit comoving volume; reduced Planck mass M)
# ---------------------------------------------------------------------------------------------
t = sp.symbols("t", real=True)
M, k, c, rv = sp.symbols("M k c rho_vac", real=True)
a = sp.Function("a")(t)
N = sp.Function("N")(t)
T = sp.Function("T")(t)
Lam = sp.Function("Lambda")(t)
rm = sp.Function("rho_m")


def L_HT(rho_vac, Lmb):
    return (-3 * M**2 * a * sp.diff(a, t)**2 / N + 3 * M**2 * k * N * a
            - N * a**3 * (rm(a) + rho_vac + Lmb) + Lmb * sp.diff(T, t))


L0 = L_HT(rv, Lam)
Lc = L_HT(rv + c, Lam - c)
diffL = sp.simplify(Lc - L0)
check("A1_compensated_shift_is_minus_c_dT", sp.simplify(diffL + c * sp.diff(T, t)) == 0, str(diffL))
# Euler-Lagrange equations: dLambda/dt = 0 (vary T), dT/dt = N a^3 (vary Lambda)
EL_T = sp.diff(sp.diff(L0, sp.diff(T, t)), t) - sp.diff(L0, T)
EL_L = sp.diff(L0, Lam)
check("A2_EL_T_gives_Lambda_const", sp.simplify(EL_T - sp.diff(Lam, t)) == 0, str(EL_T))
check("A3_EL_Lambda_gives_dT_eq_Na3", sp.simplify(EL_L - (sp.diff(T, t) - N * a**3)) == 0, str(EL_L))
# Lapse equation (Friedmann constraint) depends on rho_vac and Lambda only through their sum
EL_N = sp.diff(L0, N)
EL_Nc = sp.diff(Lc, N)
check("A4_constraint_depends_on_sum_only", sp.simplify(EL_N - EL_Nc) == 0, str(sp.simplify(EL_N)))

# ---------------------------------------------------------------------------------------------
# B. Page-Wootters / access toy for R1 (K6)
#   Unimodular-time Schroedinger equation i dPsi/dT = (h0 + rho_vac) Psi; Lambda_eff = -eig(h0).
# ---------------------------------------------------------------------------------------------
d = 8
X = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
h0 = (X + X.conj().T) / 2
E, V = np.linalg.eigh(h0)
Lambda_eff = -E                                   # eigenvector label, rho_vac independent
phi = rng.normal(size=d) + 1j * rng.normal(size=d)
phi /= np.linalg.norm(phi)
Ts = np.linspace(0, 30, 301)
Yo = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
O = (Yo + Yo.conj().T) / 2                       # a "rest" observable


def psi_T(Tv, rho_vac):
    return V @ (phi * np.exp(-1j * (E + rho_vac) * Tv))


cshift = 7.31e3                                   # a large constant (stands in for a heavy threshold)
dev_prob = 0.0
dev_phase = 0.0
for Tv in Ts:
    p0, p1 = psi_T(Tv, 0.0), psi_T(Tv, cshift)
    dev_prob = max(dev_prob, abs(np.vdot(p0, O @ p0) - np.vdot(p1, O @ p1)))
    dev_phase = max(dev_phase, np.linalg.norm(p1 - np.exp(-1j * cshift * Tv) * p0))
check("B1_conditional_expectations_invariant", dev_prob < 1e-9, f"max dev {dev_prob:.2e}")
check("B2_shift_is_pure_clock_phase", dev_phase < 1e-9, f"max dev {dev_phase:.2e}")
# access structure: an observer without access to T sees the T-averaged (Lambda-decohered) state
rho_access = V @ np.diag(abs(phi)**2) @ V.conj().T
Tlong = np.linspace(0, 4000, 40001)
acc = np.zeros((d, d), complex)
for Tv in Tlong:
    p = psi_T(Tv, cshift)
    acc += np.outer(p, p.conj())
acc /= len(Tlong)
check("B3_T_average_is_Lambda_mixture_and_shift_invariant",
      np.linalg.norm(acc - rho_access) < 5e-3, f"||avg - mixture|| = {np.linalg.norm(acc - rho_access):.2e}")
# rules: (i) bare-label rule "select the sector whose multiplier Lambda = Lambda_eff - rho_vac is
# closest to zero"; (ii) anchored rule "select Lambda_eff closest to zero" (uses h0's zero).
def bare_rule(rho_vac):
    lam_bare = Lambda_eff - rho_vac
    return Lambda_eff[np.argmin(abs(lam_bare))]


def anchored_rule(rho_vac):
    return Lambda_eff[np.argmin(abs(Lambda_eff))]


sel = {f"rho_vac={rvv}": (float(bare_rule(rvv)), float(anchored_rule(rvv))) for rvv in (0.0, 1.0, -2.0, cshift)}
num("B4_rule_selections_(bare,anchored)", sel)
check("B4_bare_label_rule_not_invariant", len({v[0] for v in sel.values()}) > 1, sel)
check("B5_anchored_rule_invariant", len({v[1] for v in sel.values()}) == 1, sel)

# ---------------------------------------------------------------------------------------------
# C. Which equivariant state rules rho = f(H)/Tr f(H) are invariant under H -> H + cI ? (K5)
# ---------------------------------------------------------------------------------------------
d2 = 6
Z = rng.normal(size=(d2, d2)) + 1j * rng.normal(size=(d2, d2))
H = (Z + Z.conj().T) / 2


def rule_state(f, Hm):
    w, U = np.linalg.eigh(Hm)
    fw = f(w)
    return (U * fw) @ U.conj().T / fw.sum()


rules = {
    "gibbs_beta0.7": lambda w: np.exp(-0.7 * w),
    "square": lambda w: w**2 + 1e-12,
    "logistic": lambda w: 1 / (1 + np.exp(w)),
    "lorentzian": lambda w: 1 / (1 + w**2),
    "power_law_on_shifted_positive": lambda w: (w - w.min() + 1.0)**-2,   # anchored at E0
    "ground_projector": lambda w: (w == w.min()).astype(float),           # anchored at E0
}
res = {}
for name, f in rules.items():
    devs = [np.linalg.norm(rule_state(f, H) - rule_state(f, H + cc * np.eye(d2))) for cc in (0.3, -1.7, 5.0)]
    res[name] = max(devs)
num("C1_shift_deviation_by_rule", res)
check("C1_gibbs_invariant", res["gibbs_beta0.7"] < 1e-12, res["gibbs_beta0.7"])
check("C2_fixed_nonexponential_rules_not_invariant",
      min(res["square"], res["logistic"], res["lorentzian"]) > 1e-3,
      {kk: res[kk] for kk in ("square", "logistic", "lorentzian")})
check("C3_anchored_rules_invariant",
      max(res["power_law_on_shifted_positive"], res["ground_projector"]) < 1e-12,
      {kk: res[kk] for kk in ("power_law_on_shifted_positive", "ground_projector")})
# Symbolic: f(x+c) = g(c) f(x) for all x, c  ->  f'(x) = g'(0) f(x)  ->  f = A exp(-beta x)
x = sp.symbols("x", real=True)
beta = sp.symbols("beta", real=True)
fF = sp.Function("f")
sol = sp.dsolve(sp.Eq(fF(x).diff(x), -beta * fF(x)))
check("C4_cauchy_ode_solution_exponential", sp.simplify(sol.rhs.diff(x) + beta * sol.rhs) == 0, str(sol))

# ---------------------------------------------------------------------------------------------
# D. Gibbs/KMS state of the whole with respect to the 4-volume clock (K2)
#   unimodular-time generator eigenvalue = rho_vac - Lambda_eff, so exp(-beta h) ~ exp(+beta Lambda_eff).
#   Density of states: (a) Gibbons-Hawking exp(S), S = 24 pi^2 / lam (lam = Lambda_eff/M^4);
#                      (b) flat on lam > 0.
# ---------------------------------------------------------------------------------------------
A_GH = 24 * math.pi**2
num("D0_GH_exponent_at_rho_obs", A_GH / 7.16e-121)


def lnZ_GH(beta_, lam_min, lam_max=1.0):
    # ln int_{lam_min}^{lam_max} exp(A/lam + beta lam) dlam, evaluated stably in u = 1/lam
    g = lambda u: A_GH * u + beta_ / u - 2 * np.log(u)          # dlam = du/u^2
    umin, umax = 1 / lam_max, 1 / lam_min
    gmax = max(g(umax), g(umin))
    W = umax - umin
    f_s = lambda s: np.exp(g(umax - s) - gmax)                    # peak at s = 0, width ~ 1/A
    edges = [0.0] + [e for e in (0.05, 0.2, 1.0) if e < W] + [W]
    val = sum(integrate.quad(f_s, lo, hi, limit=400)[0] for lo, hi in zip(edges[:-1], edges[1:]))
    return gmax + math.log(val)


tab = []
for lam_min in (1e-1, 1e-2, 1e-3, 1e-4):
    lz = lnZ_GH(+1.0, lam_min)
    lzb = lnZ_GH(-1.0, lam_min)
    tab.append({"lam_min": lam_min, "lnZ(beta=+1)": lz, "lnZ(beta=-1)": lzb, "A/lam_min": A_GH / lam_min})
num("D1_lnZ_GH_vs_cutoff", tab)
check("D1_GH_weighted_Gibbs_diverges_as_cutoff_to_0",
      all(abs(r[key] / r["A/lam_min"] - 1) < 0.01 for r in tab for key in ("lnZ(beta=+1)", "lnZ(beta=-1)")),
      "lnZ = A/lam_min (1 + O(1e-2)) for both signs of beta: no normalizable state as lam_min -> 0")
# fraction of probability within a factor 2 of the cutoff
lam_min = 1e-3
p_near = math.exp(lnZ_GH(-1.0, lam_min, 2 * lam_min) - lnZ_GH(-1.0, lam_min))
num("D2_GH_Gibbs_prob_within_factor2_of_cutoff", p_near)
check("D2_GH_Gibbs_concentrates_at_cutoff", p_near > 0.999, p_near)
# flat density: normalizable only for beta < 0; mean = 1/|beta| (value relocated to beta)
lamS, bS = sp.symbols("lambda b", positive=True)
Zflat = sp.integrate(sp.exp(-bS * lamS), (lamS, 0, sp.oo))
mean_flat = sp.simplify(sp.integrate(lamS * sp.exp(-bS * lamS), (lamS, 0, sp.oo)) / Zflat)
check("D3_flat_density_mean_equals_1_over_absbeta", sp.simplify(mean_flat - 1 / bS) == 0, str(mean_flat))
# self-consistent KMS: beta_T = kappa M^-a H^(a-4); require <Lambda_eff> = 1/beta_T = 3 M^2 H^2
Hs, Ms, kap, aa = sp.symbols("H M kappa a_exp", positive=True)
results_D4 = {}
for a_exp in (0, 1, 2, 3, 4):
    eq = sp.Eq(Ms**a_exp * Hs**(4 - a_exp) / kap, 3 * Ms**2 * Hs**2)
    s = sp.solve(eq, Hs) if a_exp != 2 else "identity iff kappa = 1/3" if sp.simplify(
        eq.lhs.subs(kap, sp.Rational(1, 3)) - eq.rhs) == 0 else "no solution"
    results_D4[f"a={a_exp}"] = str(s)
num("D4_selfconsistent_KMS_solutions", results_D4)
check("D4_selfconsistent_KMS_is_planckian_or_identity",
      all(("M" in v and "kappa" in v) or "identity" in v or v == "[]" for v in results_D4.values()), results_D4)
num("D4_kappa_needed_for_observed_H_a0", (math.sqrt(3.48e-121))**2 / 3)

# ---------------------------------------------------------------------------------------------
# E. K1: the exponent p in rho/M^4 = (V M^4)^(-p) that reproduces rho_obs (retrospective)
# ---------------------------------------------------------------------------------------------
GeV_per_100kmsMpc = 2.1332e-42
h = 0.6736
Om, Or = 0.3153, 4.18e-5 / h**2
OL = 1 - Om - Or
Mred = 2.435e18
H0 = h * GeV_per_100kmsMpc
rho_obs = 3 * OL * (H0 / Mred)**2
num("E0_rho_obs_over_M4", rho_obs)
num("E0_inputs", {"h": h, "Omega_m": Om, "Omega_r": Or, "Omega_L": OL, "M_red_GeV": Mred})
Ez = lambda aa_: math.sqrt(Om / aa_**3 + Or / aa_**4 + OL)
# conformal time tau(a) = int_0^a da'/(a'^2 E) in units 1/H0; past light cone 4-volume
tau_of = lambda aa_: integrate.quad(lambda u: 1 / (u**2 * Ez(u)), 0, aa_, limit=200)[0]
tau0 = tau_of(1.0)
# V = (4 pi/3) int_0^tau0 a^4 (tau0 - tau)^3 dtau = (4 pi/3) int_0^1 a^2/E (tau0 - tau(a))^3 da
agrid = np.concatenate([np.logspace(-8, -2, 400), np.linspace(1e-2, 1, 4000)[1:]])
taus = integrate.cumulative_trapezoid(1 / (agrid**2 * np.array([Ez(u) for u in agrid])), agrid, initial=0)
taus += tau_of(agrid[0])
integrand = agrid**2 / np.array([Ez(u) for u in agrid]) * (tau0 - taus)**3
Vcone_H0 = 4 * math.pi / 3 * integrate.trapezoid(integrand, agrid)
t0_H0 = integrate.quad(lambda u: 1 / (u * Ez(u)), 0, 1, limit=200)[0]
num("E1_Vcone_H0^4", Vcone_H0)
num("E1_t0_H0", t0_H0)
VM4 = {"past_light_cone": Vcone_H0 * (Mred / H0)**4,
       "H0^-4": (Mred / H0)**4,
       "t0^4": (t0_H0)**4 * (Mred / H0)**4}
pt = {}
for name, vm4 in VM4.items():
    p_red = math.log(rho_obs) / (-math.log(vm4))
    p_nonred = math.log(rho_obs / (8 * math.pi)**2) / (-math.log(vm4 * (8 * math.pi)**2))
    alpha_red = rho_obs * math.sqrt(vm4)
    pt[name] = {"V_M4": vm4, "p_reduced": p_red, "p_nonreduced": p_nonred, "alpha_red_for_1sigma": alpha_red}
num("E2_exponent_table", pt)
allp = [v for r in pt.values() for v in (r["p_reduced"], r["p_nonreduced"])]
check("E2_K1_p_in_[0.45,0.55]_all_conventions", all(0.45 <= p <= 0.55 for p in allp),
      f"p range [{min(allp):.4f}, {max(allp):.4f}]")
# how loose is the window? amplitude factor corresponding to delta p = 0.05 at ln(V M^4) ~ 553
num("E3_amplitude_factor_spanned_by_dp_0.05", math.exp(0.05 * math.log(VM4["past_light_cone"])))
check("E3_K1_window_is_weak", math.exp(0.05 * math.log(VM4["past_light_cone"])) > 1e10,
      "a window of +/-0.05 in p admits amplitude factors > 1e10: K1 only confirms rho ~ M^2 H^2")
# exponent if Delta V ~ V (sharp-clock/Fourier minimum): rho ~ 1/V
num("E4_rho_if_DeltaV_eq_V_over_M4", 1 / VM4["past_light_cone"])
# finite clock resolution needed for the observed spread: d_T ~ V/DeltaV ~ sqrt(V M^4)
num("E5_clock_readings_needed_sqrt_VM4", math.sqrt(VM4["past_light_cone"]))

# ---------------------------------------------------------------------------------------------
# F. One global variable (R2a): S = int sqrt(-g)(M^2R/2 - Lambda + L_m) + sigma(Lambda/mu^4).
#    Lambda-equation: sigma'((Leff - rho_vac)/mu^4) = mu^4 V4(Leff).  dLeff/drho_vac = 0 iff sigma''=0.
# ---------------------------------------------------------------------------------------------
Leff, rho, mu = sp.symbols("L_eff rho_vac mu", real=True)
sig1 = sp.Function("sigma1")          # sigma'
V4 = sp.Function("V4")
F = sig1((Leff - rho) / mu**4) - mu**4 * V4(Leff)
dL_drho = sp.simplify(-sp.diff(F, rho) / sp.diff(F, Leff))
num("F1_dLeff_drho_vac", str(dL_drho))
# direct evaluation with explicit sigma: linear sigma' = v (constant) vs quadratic sigma' = v + q x
v_, q_ = sp.symbols("v q", real=True)
for label, s1expr in (("linear_sigma", lambda z: v_ + 0 * z), ("nonlinear_sigma", lambda z: v_ + q_ * z)):
    Fx = s1expr((Leff - rho) / mu**4) - mu**4 * V4(Leff)
    dx = sp.simplify(-sp.diff(Fx, rho) / sp.diff(Fx, Leff))
    OUT["numbers"][f"F2_dLeff_drho_{label}"] = str(dx)
check("F2_linear_sigma_decouples_rho_vac", OUT["numbers"]["F2_dLeff_drho_linear_sigma"] == "0",
      "sigma linear: V4(Leff) = v/mu^4, so Leff does not depend on rho_vac")
Fq = (v_ + q_ * (Leff - rho) / mu**4) - mu**4 * V4(Leff)
dq = sp.simplify(-sp.diff(Fq, rho) / sp.diff(Fq, Leff))
check("F3_nonlinear_sigma_leaks_rho_vac", sp.simplify(dq) != 0, str(dq))

# ---------------------------------------------------------------------------------------------
# G. Fixed finite 4-volume vs constant Lambda_eff > 0 (K3), and the closed-dust V4(lambda) map.
# ---------------------------------------------------------------------------------------------
# (i) E^2 a^3 = Om + Ok a + OL a^3 has no root for a > 1 when 3 OL > -Ok (then V4 = infinity)
bad = 0
grid = 0
for Ok in np.linspace(-0.05, 0.05, 41):
    for OLg in np.linspace(0.60, 0.75, 31):
        Omg = 1 - Ok - OLg
        aa_ = np.logspace(0, 6, 2000)
        fval = Omg + Ok * aa_ + OLg * aa_**3
        grid += 1
        if fval.min() <= 0:
            bad += 1
check("G1_no_recollapse_for_OL_in_[0.60,0.75]_Ok_in_[-0.05,0.05]", bad == 0, f"{bad}/{grid} grid points recollapse")
num("G1_note", "constant Lambda_eff>0 with present acceleration => V4 infinite => no finite V0 solution")
# (ii) closed dust universe, units 3M^2 = 1 and A = rho_m a^3/(3M^2) = 1:  adot^2 = A/a + lam a^2 - 1
A_d = 1.0
lam_E = 4 / (27 * A_d**2)
num("G2_lambda_Einstein_static", lam_E)


def V4_closed(lam):
    f = lambda aa_: A_d / aa_ + lam * aa_**2 - 1
    # turnaround a_max: smallest positive root of f
    amax = optimize.brentq(f, 1e-9 + 0.0, 2.0 * A_d + 10.0) if lam <= 0 else optimize.brentq(f, 1e-9, 1.5 * A_d)
    # V4 = 2 * 2 pi^2 int_0^amax a^3 / sqrt(f) da ; substitute a = amax (1 - s^2) to remove the root singularity
    def integrand(s):
        aa_ = amax * (1 - s**2)
        val = f(aa_)
        return aa_**3 / math.sqrt(max(val, 1e-300)) * 2 * s * amax
    val, err = integrate.quad(integrand, 0, 1, limit=400)
    return 4 * math.pi**2 * val, amax


Gtab = []
for frac in (-10.0, -1.0, -0.1, 0.0, 0.5, 0.9, 0.99, 0.999):
    lam = frac * lam_E
    V4v, amax = V4_closed(lam)
    Gtab.append({"lambda/lambda_E": frac, "V4": V4v, "a_max": amax})
num("G2_closed_dust_V4_vs_lambda", Gtab)
V4s = [r["V4"] for r in Gtab]
check("G2_V4_monotone_in_lambda_and_diverges_at_lambda_E",
      all(np.diff(V4s) > 0) and V4s[-1] > 3 * V4s[3], [round(v, 3) for v in V4s])
# heavy threshold: the global multiplier absorbs c; V4 depends on Lambda_eff only
num("G3_note", "V4(lambda; A) monotone -> one-to-one map V0 <-> Lambda_eff at fixed A; rho_vac absent (F2) [D]")

# ---------------------------------------------------------------------------------------------
# H. Relational subtraction (R2b): KP  M^2 G = T - <T>/4 g   vs  Tseytlin S/V:
#    M^2 G = T - <M^2 R/2 + L_m> g.  Residual Lambda_eff (w.r.t. non-vacuum tau, L_tau):
#    KP: <tau>/4 (requires <R> = 0);  Tseytlin: <tau/2 - L_tau>.
# ---------------------------------------------------------------------------------------------
tauA, LtA, RA, cA, LeffS = sp.symbols("tau_avg L_avg R_avg c Leff", real=True)
# Tseytlin: Leff = <M^2 R/2 + L_tau>, trace: M^2 <R> = 4 Leff - <tau>
solT = sp.solve(sp.Eq(LeffS, (4 * LeffS - tauA) / 2 + LtA), LeffS)[0]
check("H1_Tseytlin_residual", sp.simplify(solT - (tauA / 2 - LtA)) == 0, str(solT))
# vacuum invariance: full T = tau - 4c, L = L_tau - c: subtraction term must absorb c exactly
kp_sub = (tauA - 4 * cA) / 4          # (1/4)<T>
ts_sub_full = (4 * (tauA / 2 - LtA) - tauA) / 2 + (LtA - cA)   # <M^2R/2> + <L> with <L> = L_tau - c
check("H2_both_subtractions_absorb_c",
      sp.simplify(kp_sub - (tauA / 4 - cA)) == 0 and sp.simplify(ts_sub_full - ((tauA / 2 - LtA) - cA)) == 0,
      "subtraction = residual - c in both cases")
rho_s, w_s = sp.symbols("rho w", real=True)
tau_f = -rho_s + 3 * w_s * rho_s
Htab = {}
for fluid, wv in (("dust", 0), ("radiation", sp.Rational(1, 3))):
    for Lname, Lval in (("L=p", w_s * rho_s), ("L=-rho", -rho_s)):
        Htab[f"{fluid},{Lname},Tseytlin"] = str(sp.simplify((tau_f / 2 - Lval).subs(w_s, wv)))
    Htab[f"{fluid},KP"] = str(sp.simplify((tau_f / 4).subs(w_s, wv)))
num("H3_residual_per_unit_avg_rho", Htab)
check("H3_KP_and_Tseytlin_differ_on_matter", Htab["dust,KP"] != Htab["dust,L=-rho,Tseytlin"], Htab)

# ---------------------------------------------------------------------------------------------
# I. Energy-budget closure (R3, K4)
# ---------------------------------------------------------------------------------------------
# (i) the gravitational Hamiltonian constraint vanishes on the flat dust + Lambda solution for every Lambda
tt, Ls, Ms2, Cm = sp.symbols("t Lambda M C_m", positive=True)
w0 = sp.sqrt(Ls / (3 * Ms2**2)) * sp.Rational(3, 2)
aa_s = (Cm / Ls)**sp.Rational(1, 3) * sp.sinh(w0 * tt)**sp.Rational(2, 3)
rho_m_s = Cm / aa_s**3
Hsym = sp.diff(aa_s, tt) / aa_s
constraint = sp.simplify(3 * Ms2**2 * Hsym**2 - rho_m_s - Ls)
check("I1_hamiltonian_constraint_zero_for_every_Lambda", sp.simplify(constraint.rewrite(sp.exp)) == 0, str(constraint))
# (ii) k = +1: rho_m + rho_L = 3M^2 (H^2 + 1/a^2) > 0 on every slice -> the slice version cannot hold
num("I2_note", "k=+1: rho_m + rho_L = 3M^2(H^2+1/a^2) > 0 on every slice [D]")
# (iii) history-average version, flat dust + negative Lambda: <rho_m>_V4 / |rho_L| = 2 for every |rho_L|
#   a^3 = (C/|L|) sin^2(3 w t/2),  w = sqrt(|L|/(3M^2));  a^3 rho_m = C constant.
def ratio_dust(absL, C=1.0, M2=1.0):
    w = math.sqrt(absL / (3 * M2))
    Tc = 2 * math.pi / (3 * w)
    a3 = lambda s: (C / absL) * math.sin(1.5 * w * s)**2
    num_ = C * Tc
    den = integrate.quad(a3, 0, Tc)[0]
    return (num_ / den) / absL


rt = {f"|rho_L|={x_}": ratio_dust(x_) for x_ in (1e-3, 1.0, 1e4)}
num("I3_avg_rho_m_over_absrhoL_dust", rt)
check("I3_ratio_is_2_at_every_scale", all(abs(v - 2) < 1e-8 for v in rt.values()), rt)


# (iv) dust + radiation, flat, negative Lambda: ratio R(r) = <rho_m + rho_r>/|rho_L| as r = rho_r/rho_m at turnaround varies
def ratio_mixed(r):
    # units |rho_L| = 1, 3M^2 = 1. At turnaround a_t = 1: rho_m(1) + rho_r(1) = 1, rho_r = r rho_m.
    rm1 = 1 / (1 + r)
    rr1 = r / (1 + r)
    H2 = lambda aa_: rm1 / aa_**3 + rr1 / aa_**4 - 1.0
    # V4 ~ 2 int_0^1 a^3 da / (a sqrt(H2)) ;  matter+rad energy weight: 2 int a^3 (rho) da/(a sqrt(H2))
    def I(fn):
        g = lambda s: fn(1 - s**2) / ((1 - s**2) * math.sqrt(max(H2(1 - s**2), 1e-300))) * 2 * s
        return integrate.quad(g, 0, 1, limit=400)[0]
    numer = I(lambda aa_: aa_**3 * (rm1 / aa_**3 + rr1 / aa_**4))
    denom = I(lambda aa_: aa_**3)
    return numer / denom


Rt = {f"r={r}": ratio_mixed(r) for r in (0.0, 0.01, 0.1, 1.0, 10.0)}
num("I4_ratio_mixed_vs_r", Rt)
check("I4_ratio_scale_free_and_above_1", all(v > 1.0 for v in Rt.values()) and abs(Rt["r=0.0"] - 2) < 1e-3, Rt)
num("I5_note", "rho_L = -<rho_matter> <= 0 for WEC matter [D]; observed rho_L > 0 (retrospective)")

# ---------------------------------------------------------------------------------------------
n_pass = sum(v["pass"] for v in OUT["checks"].values())
OUT["summary"] = f"{n_pass}/{len(OUT['checks'])} checks pass"
print(OUT["summary"])
with open(sys.argv[1] if len(sys.argv) > 1 else "grut_native_checks.json", "w") as fh:
    json.dump(OUT, fh, indent=1, default=str)
