#!/usr/bin/env python3
"""Checks for track `unimodular` (packet grut_vacuum_principle_20261021, issue #13, Codex 02b4515).

Standard physics only.  Nothing here is a new theory, an empirical validation or a confirmatory
test.  Every number built from H0, Omega_m, Omega_Lambda, BBN or CMB levels is RETROSPECTIVE
(already-exposed data).  Conventions: signature (-,+,+,+), hbar = c = 1,
reduced Planck mass M = (8 pi G)^(-1/2) = 2.435e18 GeV, l_r = 1/M; non-reduced M_P = G^(-1/2).

Sections
  A  trace-free (unimodular) equations: a constant vacuum term drops out; Lambda is an
     integration constant fixed by slice data, not by the theory
  B  Henneaux-Teitelboim minisuperspace: the vacuum shift is a boundary phase e^{-i c V4};
     the Lambda-translation lemma (kinematic closure structure cannot fix <Lambda>);
     what breaks the translation (Lambda-dependent gravitational weights)
  C  Ng-van Dam / Smolin "on-shell action" condition: dependence on the on-shell matter Lagrangian
  D  everpresent Lambda magnitude: past-light-cone 4-volume for radiation, matter and Planck-2018
     LCDM; V^(-1/2) against rho_Lambda / M^4 in two unit conventions
  E  Gaussian (linear) statistics of the everpresent ansatz: sigma(Omega_Lambda) per epoch,
     lag correlation sqrt(V1/V2)
  F  small vectorised Monte Carlo of the published "Model 1" recursion (Ahmed et al. 2004;
     Das et al. 2023) to check E with Friedmann feedback and the H^2 < 0 terminations
  G  Barrow-type CMB estimate for horizon-scale independent Lambda patches (order of magnitude)
  H  heavy-threshold bookkeeping: what is stable and what is a value statement

Run:  python3 unimodular_checks.py      (numpy, scipy, sympy; about one minute; seed 20261021)
"""
import json
import math
import os

import numpy as np
import sympy as sp
from scipy import integrate, special

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "unimodular_checks.json")
SEED = 20261021
rng = np.random.default_rng(SEED)
R = {}
checks = []


def record(name, value, tol, note=""):
    ok = bool(abs(value) <= tol)
    checks.append({"name": name, "max_abs_deviation": float(value), "tol": tol, "pass": ok,
                   "note": note})
    return ok


def record_bool(name, ok, note=""):
    checks.append({"name": name, "pass": bool(ok), "note": note})
    return ok


# ------------------------------------------------------------------------------------------
# A. Trace-free equations
# ------------------------------------------------------------------------------------------
# A1: TF(X) = X - (1/4) g tr(g^-1 X) is blind to X -> X - C g, for a random Lorentzian metric.
C = sp.symbols("C", real=True)


def rand_sym_rational(n, lo=-5, hi=5):
    m = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            v = sp.Rational(int(rng.integers(lo, hi + 1)), int(rng.integers(1, 4)))
            m[i, j] = v
            m[j, i] = v
    return m


while True:
    g = rand_sym_rational(4)
    g = g + sp.diag(-12, 9, 9, 9)          # make it Lorentzian and well conditioned
    ev = np.linalg.eigvalsh(np.array(g, dtype=float))
    if (ev < 0).sum() == 1 and abs(np.linalg.det(np.array(g, dtype=float))) > 1:
        break
T = rand_sym_rational(4)
ginv = g.inv()


def TF(X):
    return X - sp.Rational(1, 4) * g * (ginv * X).trace()


diffA1 = sp.simplify(TF(T - C * g) - TF(T))
record_bool("A1_tracefree_part_blind_to_C_g", diffA1 == sp.zeros(4, 4),
            "TF(T - C g) - TF(T) = 0 symbolically in C, random rational Lorentzian g")

# A2: FRW. TF 00-equation and conservation imply d/dt [3 M^2 (H^2 + k/a^2) - rho] = 0,
#     and rho may contain any rho_vac (p_vac = -rho_vac) without changing the TF equation.
t, M, k, rv = sp.symbols("t M k rho_vac", real=True)
a = sp.Function("a")(t)
rho = sp.Function("rho")(t)      # non-vacuum density
p = sp.Function("p")(t)          # non-vacuum pressure
Hs = sp.diff(a, t) / a
add = sp.diff(a, t, 2)
# components (mostly plus): R_00 = -3 a''/a ; R = 6(a''/a + H^2 + k/a^2)
Rsc = 6 * (add / a + Hs ** 2 + k / a ** 2)
TF00_R = -3 * add / a + Rsc / 4
rho_tot, p_tot = rho + rv, p - rv
Ttrace = -rho_tot + 3 * p_tot
TF00_T = rho_tot + Ttrace / 4
tf_eq = sp.simplify(M ** 2 * TF00_R - TF00_T)                  # = 0 on solutions
record_bool("A2a_TF_equation_independent_of_rho_vac", sp.diff(tf_eq, rv) == 0,
            "M^2 TF_00(R) - TF_00(T) has no rho_vac dependence")
add_sol = sp.solve(sp.Eq(tf_eq, 0), add)[0]
Lam_hat = 3 * M ** 2 * (Hs ** 2 + k / a ** 2) - rho_tot       # candidate integration constant
dLam = sp.diff(Lam_hat, t)
dLam = dLam.subs(sp.diff(a, t, 2), add_sol)
dLam = dLam.subs(sp.diff(rho, t), -3 * Hs * (rho + p))       # non-vacuum conservation
record_bool("A2b_Lambda_hat_conserved", sp.simplify(dLam) == 0,
            "d/dt[3M^2(H^2+k/a^2) - rho - rho_vac] = 0 using the TF equation and conservation")
rho_eff = sp.simplify(Lam_hat + rv)                           # what gravitates: rho_vac + Lambda_hat
record_bool("A2c_rho_eff_equals_slice_data", sp.simplify(rho_eff - (3 * M ** 2 * (Hs ** 2 + k / a ** 2) - rho)) == 0,
            "rho_vac + Lambda_hat = 3M^2(H^2+k/a^2) - rho_nonvac: set by slice data, not by rho_vac")
# general-tensor version of the same integration constant (Weinberg 7.6-7.7):
# M^2 R + T = 4 Lambda_hat_tot and M^2 G = T - Lambda_hat_tot g;  with T = tau - rho_vac g,
# rho_eff = Lambda_hat_tot + rho_vac = (M^2 R + tau^mu_mu)/4.  Check the FRW reduction:
tau_trace = -rho + 3 * p
expr = sp.simplify((M ** 2 * Rsc + tau_trace) / 4 - rho_eff)
expr = sp.simplify(expr.subs(sp.diff(a, t, 2), add_sol))
record_bool("A2d_rho_eff_equals_(M2R+tau)/4", expr == 0, "FRW reduction of rho_eff = (M^2 R + tau)/4")

# A3: numerical integration of the TF ODE (second order in a) with dust + radiation + a huge
#     rho_vac: the expansion history is identical for rho_vac = 0 and rho_vac = 1e6 (code units),
#     and the gravitating constant equals 3H^2 - rho at the initial slice.
def tf_rhs(tt, y, rv_val, om, orad):
    av, adv = y
    rho_n = 3 * (om * av ** -3 + orad * av ** -4)
    p_n = 3 * orad * av ** -4 / 3
    # a'' = a (H^2 - (rho+p)/2) for k = 0, M = 1 ; rho_vac never enters (rho_vac + p_vac = 0)
    H = adv / av
    return [adv, av * (H ** 2 - (rho_n + p_n) / 2.0 + 0.0 * rv_val)]


om_, orad_ = 0.3, 1e-4
a0, H0c = 0.05, None
rho0 = 3 * (om_ * a0 ** -3 + orad_ * a0 ** -4)
Lam_target = 2.1            # chosen initial-slice value of the gravitating constant (code units)
H0c = math.sqrt((rho0 + Lam_target) / 3.0)
sols = {}
for rv_val in (0.0, 1.0e6):
    s = integrate.solve_ivp(tf_rhs, (0, 1.2), [a0, a0 * H0c], args=(rv_val, om_, orad_),
                            rtol=1e-11, atol=1e-13, dense_output=True)
    sols[rv_val] = s
tt = np.linspace(0, 1.2, 50)
dev = np.max(np.abs(sols[0.0].sol(tt)[0] - sols[1.0e6].sol(tt)[0]))
record("A3a_TF_history_independent_of_rho_vac", dev, 1e-12, "max |a(t; rho_vac=0) - a(t; 1e6)|")
yy = sols[0.0].sol(tt)
Hh = yy[1] / yy[0]
Lam_along = 3 * Hh ** 2 - 3 * (om_ * yy[0] ** -3 + orad_ * yy[0] ** -4)
record("A3b_integration_constant_constant_and_equal_initial_slice", np.max(np.abs(Lam_along - Lam_target)), 1e-7,
       "3H^2 - rho_nonvac along the TF solution equals its initial-slice value 2.1")
R["A3_note"] = ("In the trace-free theory the gravitating constant is whatever the initial slice says: "
                "any value is a solution. Fixing it from (H0, rho_m0) is a measurement, not a prediction.")

# ------------------------------------------------------------------------------------------
# B. Henneaux-Teitelboim minisuperspace and the Lambda-translation lemma
# ------------------------------------------------------------------------------------------
N_, Lam, cc = sp.symbols("N Lambda c", real=True)
aT = sp.Function("a")(t)
TT = sp.Function("T")(t)
LamF = sp.Function("Lam")(t)
Nf = sp.Function("N")(t)
rhom0, rhor0 = sp.symbols("rho_m0 rho_r0", positive=True)


def L_HT(rvac, lamf):
    rho_a = rhom0 * aT ** -3 + rhor0 * aT ** -4
    return (-3 * M ** 2 * aT * sp.diff(aT, t) ** 2 / Nf + 3 * M ** 2 * k * Nf * aT
            - Nf * aT ** 3 * (rho_a + rvac + lamf) + lamf * sp.diff(TT, t))


L0 = L_HT(rv, LamF)
eqs = sp.euler_equations(L0, [aT, TT, LamF, Nf], t)
eq_T = sp.simplify(eqs[1].lhs)        # from T:  -d/dt Lam = 0
eq_L = sp.simplify(eqs[2].lhs)        # from Lam: T' - N a^3 = 0
eq_N = sp.simplify(eqs[3].lhs)        # constraint
R["B1_HT_equations"] = {"dL/dT": str(eq_T), "dL/dLambda": str(eq_L), "dL/dN": str(eq_N)}
record_bool("B1a_Lambda_constant_and_dT_equals_4volume_rate",
            sp.simplify(eq_T + sp.diff(LamF, t)) == 0 and sp.simplify(eq_L - (sp.diff(TT, t) - Nf * aT ** 3)) == 0,
            "dLambda/dt = 0 and dT/dt = N a^3 (T counts 4-volume)")
# shift (rho_vac, Lambda) -> (rho_vac + c, Lambda - c): Lagrangian changes by -c dT/dt only
Lsh = L_HT(rv + cc, LamF - cc)
dL = sp.simplify(sp.expand(Lsh - L0))
record_bool("B1b_vacuum_shift_is_boundary_term", sp.simplify(dL + cc * sp.diff(TT, t)) == 0,
            "L(rho_vac+c, Lambda-c) - L(rho_vac, Lambda) = -c dT/dt: action changes by -c * (4-volume)")
# the constraint depends only on Lambda_eff = Lambda + rho_vac
eq_N_sh = sp.simplify(sp.euler_equations(Lsh, [aT, TT, LamF, Nf], t)[3].lhs)
record_bool("B1c_constraint_depends_only_on_Lambda_plus_rho_vac", sp.simplify(eq_N_sh - eq_N) == 0,
            "Friedmann constraint unchanged under the compensated shift")

# B2: canonical form. With p_T = Lambda the constraint reads  a^3 p_T + K(a, p_a) + a^3 rho_vac = 0,
# so the unimodular-time Schroedinger generator is h = -(K/a^3) - rho_vac: h_c = h - c I exactly.
# Hence psi_c(T) = e^{+icT} psi(T) (sign convention-dependent): Codex's fixed-background identity
# with the external t replaced by unimodular time T.  Numerical illustration with a random
# Hermitian stand-in for -K/a^3 (finite truncation; the identity does not depend on K).
d = 12
X = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
hK = (X + X.conj().T) / 2
psi0 = rng.normal(size=d) + 1j * rng.normal(size=d)
psi0 /= np.linalg.norm(psi0)
Tval, cval = 3.7, 11.3
w, U = np.linalg.eigh(hK)
ev_T = lambda hshift: U @ (np.exp(-1j * (w + hshift) * Tval) * (U.conj().T @ psi0))
psiA, psiB = ev_T(0.0), ev_T(-cval)
record("B2a_probabilities_in_unimodular_time_blind_to_c", np.max(np.abs(np.abs(psiA) ** 2 - np.abs(psiB) ** 2)), 1e-12)
record("B2b_relative_phase_is_e^{icT}", abs(np.vdot(psiA, psiB) - np.exp(1j * cval * Tval)), 1e-12)

# B3: Lambda-translation lemma.  L^2(dLambda), [T, Lambda] = i, W_c = e^{icT} translates Lambda -> Lambda + c.
# Every translation-covariant functional (norm, variance, sigma_T, entropy of |A|^2,
# uncertainty product, |psi(T)|^2) is unchanged; only <Lambda> moves.
Lg = np.linspace(-40, 40, 2 ** 14)
dLg = Lg[1] - Lg[0]


def packet(mean, sig, chirp=0.3):
    A = np.exp(-(Lg - mean) ** 2 / (4 * sig ** 2) + 1j * chirp * (Lg - mean) ** 2)
    return A / math.sqrt(np.sum(np.abs(A) ** 2) * dLg)


def stats(A):
    P = np.abs(A) ** 2 * dLg
    mL = np.sum(Lg * P)
    sL = math.sqrt(np.sum((Lg - mL) ** 2 * P))
    # T representation: psi(T) = int dLambda A(Lambda) e^{-i Lambda T}
    Tg = np.fft.fftshift(np.fft.fftfreq(Lg.size, d=dLg)) * 2 * np.pi
    psiT = np.fft.fftshift(np.fft.fft(A * np.exp(-1j * Lg[0] * 0)))
    PT = np.abs(psiT) ** 2
    PT /= PT.sum()
    mT = np.sum(Tg * PT)
    sT = math.sqrt(np.sum((Tg - mT) ** 2 * PT))
    Pn = P[P > 1e-300]
    ent = float(-np.sum(Pn * np.log(Pn / dLg)))
    return dict(norm=float(P.sum()), mean=float(mL), sigL=float(sL), sigT=float(sT), prod=float(sL * sT),
                ent=ent, PT=PT)


sA = stats(packet(0.0, 1.3))
sB = stats(packet(7.5, 1.3))
record("B3a_norm", abs(sA["norm"] - sB["norm"]), 1e-12)
record("B3b_sigma_Lambda", abs(sA["sigL"] - sB["sigL"]), 1e-10)
record("B3c_sigma_T", abs(sA["sigT"] - sB["sigT"]), 1e-8)
record("B3d_entropy_of_|A|^2", abs(sA["ent"] - sB["ent"]), 1e-9)
record("B3e_|psi(T)|^2_profile", float(np.max(np.abs(sA["PT"] - sB["PT"]))), 1e-10)
record("B3f_mean_moves_by_c", abs((sB["mean"] - sA["mean"]) - 7.5), 1e-10)
s_min = stats(packet(0.0, 1.3, chirp=0.0))
record("B3g_unchirped_gaussian_saturates_sigL_sigT_half", abs(s_min["prod"] - 0.5), 1e-6)
R["B3_lemma"] = {"sigma_L": sA["sigL"], "sigma_T": sA["sigT"], "product": sA["prod"],
                 "product_min_uncertainty": s_min["prod"], "shift_c": 7.5,
                 "statement": ("A principle built only from the kinematic closure structure of quantum "
                               "unimodular gravity (Hilbert space L^2(dLambda), the HT inner product, unitarity in "
                               "unimodular time, normalizability, uncertainty relations, entropies) is invariant under "
                               "Lambda -> Lambda + c and cannot select <Lambda>. It can select translation-invariant "
                               "features (spread, shape). Fixing the mean needs a Lambda-dependent weight from the "
                               "gravitational dynamics.")}

# B4: weights that break the translation (supplied by gravitational dynamics, not by closure).
# Euclidean de Sitter action S_E = -24 pi^2 M^4 / rho_L  (8 pi G = 1/M^2) -> Baum-Hawking weight exp(+24 pi^2/x),
# x = rho_L / M^4.  Not translation covariant, peaks at x -> 0+, not normalizable at 0+.
x_obs = None  # filled after D
# B5: everpresent spread from a saturated Gaussian (Das et al. 2023 eq. 3.5): sigma_rho * sigma_V = 1/2,
# sigma_V = sqrt(N) l^4 = sqrt(V) l^2  ->  sigma_rho = (1/2) / (l^2 sqrt V).  In reduced units
# (rho/M^4 versus V M^4) the coefficient is alpha_red = (1/2) (l_r / l)^2.
def alpha_red_of_l(l_over_lr):
    return 0.5 / l_over_lr ** 2


R["B5_alpha_red_of_discreteness_length"] = {
    "l = l_r (reduced Planck length)": alpha_red_of_l(1.0),
    "l = l_P = l_r/sqrt(8 pi)": alpha_red_of_l(1 / math.sqrt(8 * math.pi)),
    "Sorkin 2007 ansatz Lambda = (1/l^2)/sqrt(V), l^2 ~ 3 kappa (his text): alpha_red": 1.0 / 3.0,
}

# ------------------------------------------------------------------------------------------
# C. Ng-van Dam / Smolin condition: on-shell action per unit 4-volume = 0
# ------------------------------------------------------------------------------------------
rhoL, rr, LmS = sp.symbols("rho_L rho L_m", real=True)
# S = int sqrt(-g)[M^2 R/2 - rho_L + L_m]; trace: M^2 R = 4 rho_L - T; on-shell density = rho_L + L_m - T/2
onshell = sp.simplify((4 * rhoL - sp.Symbol("Ttr")) / 2 - rhoL + LmS)
record_bool("C1_onshell_density", sp.simplify(onshell - (rhoL + LmS - sp.Symbol("Ttr") / 2)) == 0,
            "M^2R/2 - rho_L + L_m = rho_L + L_m - T/2 on shell (trace equation M^2 R = 4 rho_L - T)")
# Using M^2 R/2 = 2 rho_L - T/2:   density = rho_L + L_m - T/2
dens = lambda Lm, Ttr: rhoL + Lm - Ttr / 2
cases = {
    "dust, L_m = p = 0": dens(0, -rr),
    "dust, L_m = -rho": dens(-rr, -rr),
    "radiation, L_m = p = rho/3": dens(rr / 3, 0),
    "radiation, L_m = -rho": dens(-rr, 0),
}
solC = {kk: str(sp.solve(sp.Eq(vv, 0), rhoL)[0]) for kk, vv in cases.items()}
R["C2_S_onshell_zero_gives_rho_L"] = solC
record_bool("C2_sign_depends_on_onshell_matter_Lagrangian",
            sp.solve(sp.Eq(cases["dust, L_m = p = 0"], 0), rhoL)[0] == -rr / 2 and
            sp.solve(sp.Eq(cases["dust, L_m = -rho"], 0), rhoL)[0] == rr / 2,
            "rho_L = -<rho>/2 (L=p) versus +<rho>/2 (L=-rho) for dust: the 'prediction' is not well defined")
R["C3_explicit_dS_dLambda"] = "dS/d(rho_L) at fixed history = -V4 != 0: no stationary point in rho_L at finite V4; " \
                              "the rho_L integral with flat measure instead imposes delta(V4 - V_fixed)"

# ------------------------------------------------------------------------------------------
# D. Everpresent Lambda magnitude
# ------------------------------------------------------------------------------------------
GeV_per_s = 6.582119569e-25            # hbar in GeV s
Mpc_km = 3.0856775814913673e19
M_red = 2.435e18                         # GeV
M_P = M_red * math.sqrt(8 * math.pi)     # 1.2209e19 GeV
Gyr_s = 3.15576e16

# Planck 2018 TT,TE,EE+lowE+lensing (retrospective inputs)
h = 0.6736
Om = 0.3153
Tcmb = 2.7255
Neff = 3.046
og = 2.4728e-5 * (Tcmb / 2.7255) ** 4 / h ** 2
Or = og * (1 + 0.2271 * Neff)
OL = 1 - Om - Or
H0_GeV = 100 * h / Mpc_km * GeV_per_s
R["D0_inputs"] = {"h": h, "Omega_m": Om, "Omega_r": Or, "Omega_L": OL, "H0_GeV": H0_GeV,
                  "status": "retrospective: Planck 2018 base LCDM; massless-neutrino radiation approximation"}

# D1: power-law eras, analytic H^4 V via Beta functions, and quadrature check
def H4V_powerlaw(pw):
    # a = t^p:  V = (4pi/3) int_0^t a^3 (int_t'^t dt''/a)^3 dt' ;  H = p/t
    f = lambda tp: tp ** (3 * pw) * ((1 - tp ** (1 - pw)) / (1 - pw)) ** 3
    val = integrate.quad(f, 0, 1, epsabs=1e-14, epsrel=1e-12)[0]
    return (4 * math.pi / 3) * val * pw ** 4


H4V_rad_exact = (64 * math.pi / 3) * special.beta(5, 4) / 16
H4V_mat_exact = 108 * math.pi * special.beta(9, 4) * 16 / 81
record("D1a_H4V_radiation", H4V_powerlaw(0.5) - H4V_rad_exact, 1e-12, "4 pi / 840 = 0.0149600")
record("D1b_H4V_matter", H4V_powerlaw(2 / 3) - H4V_mat_exact, 1e-12, "matter era 0.0338484")
record("D1c_radiation_equals_4pi_over_840", H4V_rad_exact - 4 * math.pi / 840, 1e-15)
R["D1_H4V"] = {"radiation": H4V_rad_exact, "matter": H4V_mat_exact,
               "sigma_Omega_per_alpha_red_radiation": 1 / (3 * math.sqrt(H4V_rad_exact)),
               "sigma_Omega_per_alpha_red_matter": 1 / (3 * math.sqrt(H4V_mat_exact))}

# D2: Planck-2018 LCDM: conformal time, past-light-cone 4-volume V(a) (units H0 = 1)
def E(aa):
    return np.sqrt(Or * aa ** -4 + Om * aa ** -3 + OL)


lna = np.linspace(math.log(1e-14), math.log(1.0), 400001)
aa = np.exp(lna)
dtau_dlna = 1.0 / (aa * E(aa))                       # dtau = da/(a^2 H) = dlna/(a H)
tau = integrate.cumulative_trapezoid(dtau_dlna, lna, initial=0.0) + aa[0] / math.sqrt(Or)
dt_dlna = 1.0 / E(aa)
tcos = integrate.cumulative_trapezoid(dt_dlna, lna, initial=0.0) + aa[0] ** 2 / (2 * math.sqrt(Or))
# moments I_j = int a^4 tau^j dtau
w4 = aa ** 4 * dtau_dlna
I = []
for j in range(4):
    Ij = integrate.cumulative_trapezoid(w4 * tau ** j, lna, initial=0.0)
    Ij += Or ** 2 * tau[0] ** (5 + j) / (5 + j)                  # radiation-era head start
    I.append(Ij)
Vcone = (4 * math.pi / 3) * (tau ** 3 * I[0] - 3 * tau ** 2 * I[1] + 3 * tau * I[2] - I[3])
H4V = E(aa) ** 4 * Vcone
# convergence of the cancellation: compare with direct quadrature at a few epochs
def V_direct(a_obs):
    tau_obs = np.interp(math.log(a_obs), lna, tau)
    f = lambda la: math.exp(la) ** 4 * (tau_obs - np.interp(la, lna, tau)) ** 3 / (math.exp(la) * E(math.exp(la)))
    return (4 * math.pi / 3) * integrate.quad(f, math.log(1e-14), math.log(a_obs), limit=400, epsrel=1e-10)[0]


for a_obs in (1e-9, 1 / 1090.8, 0.5, 1.0):
    vd = V_direct(a_obs)
    vm = np.interp(math.log(a_obs), lna, Vcone)
    record(f"D2_V_moments_vs_direct_a={a_obs:.4g}", abs(vm / vd - 1), 2e-4)
i8 = int(np.argmin(np.abs(aa - 1e-8)))
record("D2_radiation_limit_H4V", abs(H4V[i8] / H4V_rad_exact - 1), 2e-3, "a = 1e-8, deep radiation era")

epochs = {"BBN (a = 1e-9, T ~ 0.2 MeV)": 1e-9, "matter-radiation equality": Or / Om,
          "recombination z = 1089.8": 1 / 1090.8, "z = 2": 1 / 3, "z = 1": 0.5, "z = 0.5": 1 / 1.5,
          "today z = 0": 1.0}
ep = {}
for name, ae in epochs.items():
    i = int(np.argmin(np.abs(aa - ae)))
    ep[name] = {"a": float(aa[i]), "H4V": float(H4V[i]),
                "sigma_Omega_Lambda_per_alpha_red": float(1 / (3 * math.sqrt(H4V[i]))),
                "H0^4 V": float(Vcone[i]), "t*H0": float(tcos[i])}
R["D2_epochs_LCDM"] = ep
H4V0 = float(H4V[-1])
V0_H0units = float(Vcone[-1])
t0H0 = float(tcos[-1])
R["D2_today"] = {"H0^4 V_cone": V0_H0units, "t0 H0": t0H0, "t0_Gyr": t0H0 / (H0_GeV / GeV_per_s) / Gyr_s,
                 "tau0 H0 (comoving particle horizon, units c/H0)": float(tau[-1])}

# magnitudes in two conventions
x_obs = 3 * OL * (H0_GeV / M_red) ** 2                   # rho_L / M^4 (reduced)
VM4 = V0_H0units * (M_red / H0_GeV) ** 4                 # V in reduced-Planck units
VP4 = V0_H0units * (M_P / H0_GeV) ** 4                   # V in G = 1 Planck units
rho_over_MP4 = x_obs * (M_red / M_P) ** 4
R["D3_magnitudes"] = {
    "rho_L_obs / M^4 (reduced, 8piG=1)": x_obs,
    "rho_L_obs / M_P^4 (G=1)": rho_over_MP4,
    "Lambda l_P^2 (G=1 curvature)": 8 * math.pi * rho_over_MP4,
    "(H0/M)^2": (H0_GeV / M_red) ** 2,
    "(t_P/t0)^2 (Barrow's estimate, t_P = sqrt(G))": (1 / M_P / (t0H0 / H0_GeV)) ** 2,
    "V_cone in reduced units, V M^4": VM4,
    "V_cone in G=1 units, V/l_P^4": VP4,
    "(V M^4)^(-1/2)": VM4 ** -0.5,
    "(V/l_P^4)^(-1/2)": VP4 ** -0.5,
    "alpha_red* = (rho_L/M^4) (V M^4)^(1/2)  [1-sigma coefficient, reduced convention]": x_obs * math.sqrt(VM4),
    "alpha_G1* = (rho_L/M_P^4) (V/l_P^4)^(1/2)  [G=1 convention, Das et al. alpha]": rho_over_MP4 * math.sqrt(VP4),
    "N = V/l_P^4 (Sorkin's 'number of elements' if l = l_P)": VP4,
}
record("D3a_alpha_conventions_differ_by_8pi",
       abs(x_obs * math.sqrt(VM4) / (rho_over_MP4 * math.sqrt(VP4)) - 8 * math.pi), 1e-9)
record("D3b_alpha_red*_closed_form_3OmegaL_sqrt(H4V)", abs(x_obs * math.sqrt(VM4) - 3 * OL * math.sqrt(H4V0)), 1e-12)
# alternative 4-volume proxies (all in H0 units)
proxies = {"past light cone (computed)": V0_H0units, "H0^-4": 1.0, "t0^4": t0H0 ** 4,
           "(4pi/3) H0^-3 t0 (Hubble volume x age)": 4 * math.pi / 3 * t0H0,
           "(4pi/3) tau0^3 t0 (particle-horizon volume x age, comoving=physical today)": 4 * math.pi / 3 * float(tau[-1]) ** 3 * t0H0}
R["D4_proxies"] = {kk: {"H0^4 V": vv, "(V M^4)^(-1/2)": (vv * (M_red / H0_GeV) ** 4) ** -0.5,
                        "ratio rho_obs/(that)": x_obs / (vv * (M_red / H0_GeV) ** 4) ** -0.5}
                   for kk, vv in proxies.items()}
# discreteness length implied by alpha_red* in the saturated-Gaussian normalisation (B5)
ar_star = x_obs * math.sqrt(VM4)
R["D5_discreteness_length_for_1sigma"] = {"l/l_r": math.sqrt(0.5 / ar_star),
                                          "l/l_P": math.sqrt(0.5 / ar_star) * math.sqrt(8 * math.pi)}
# B4 numbers now that x_obs is known
R["B4_Baum_Hawking"] = {"24 pi^2 M^4 / rho_L_obs (exponent = dS entropy)": 24 * math.pi ** 2 / x_obs,
                        "ln[w(rho/2)/w(rho)]": 24 * math.pi ** 2 / x_obs,
                        "note": "weight exp(+24 pi^2 M^4/rho_L) on rho_L>0 is maximal as rho_L -> 0+, not normalisable; "
                                "it fixes the sign (+) and drives the scale to zero"}

# ------------------------------------------------------------------------------------------
# E. Gaussian statistics: sigma(Omega) per epoch and lag correlations (fixed LCDM background)
# ------------------------------------------------------------------------------------------
def corr(a1, a2):
    return math.sqrt(np.interp(math.log(a1), lna, Vcone) / np.interp(math.log(a2), lna, Vcone))


R["E1_lag_correlation_sqrt(V1/V2)"] = {
    "z=0.5 vs 0": corr(1 / 1.5, 1.0), "z=1 vs 0": corr(0.5, 1.0), "z=2 vs 0": corr(1 / 3, 1.0),
    "rec vs 0": corr(1 / 1090.8, 1.0), "BBN vs 0": corr(1e-9, 1.0),
    "per e-fold in a, matter era (analytic a^-3)": math.exp(-3.0),
    "per e-fold in a, radiation era (analytic a^-4)": math.exp(-4.0)}
# analytic: random walk S(V) with Var = alpha^2 V; rho = S/V; Cov(rho1, rho2) = alpha^2/V2 (V1 < V2)
record("E1a_matter_era_corr_power", abs(math.log(corr(0.3, 0.6)) / math.log(0.5) - 3) , 0.35,
       "late LCDM is not pure matter; loose check of the a^3 law")
# E2: linear-Gaussian joint probability (fixed LCDM background, independent epochs since corr ~ 1e-9,
#     no H^2 <= 0 termination, no Omega nonlinearity): an optimistic estimate maximised over alpha.
from scipy.stats import norm
s0 = ep["today z = 0"]["sigma_Omega_Lambda_per_alpha_red"]
srec = ep["recombination z = 1089.8"]["sigma_Omega_Lambda_per_alpha_red"]
sbbn = ep["BBN (a = 1e-9, T ~ 0.2 MeV)"]["sigma_Omega_Lambda_per_alpha_red"]
best = None
for al in np.logspace(-3, 0.5, 4000):
    p0 = norm.cdf(0.74, 0, s0 * al) - norm.cdf(0.64, 0, s0 * al)
    pr = 2 * norm.cdf(0.05 / (srec * al)) - 1
    pb = 2 * norm.cdf(0.07 / (sbbn * al)) - 1
    if best is None or p0 * pr * pb > best[0]:
        best = (p0 * pr * pb, al, p0, pr, pb)
R["E2_linear_gaussian_joint_max"] = {"P_max": best[0], "alpha_red_at_max": best[1], "P(today in [0.64,0.74])": best[2],
                                     "P(|Omega_rec|<0.05)": best[3], "P(|Omega_BBN|<0.07)": best[4],
                                     "BBN-compatible alpha (2 sigma: 0.07/(2*sigma_BBN))": 0.07 / (2 * sbbn),
                                     "Omega_0 = 0.685 in units of sigma_Omega(today) at that alpha":
                                         0.685 / (s0 * 0.07 / (2 * sbbn)),
                                     "status": "optimistic: ignores termination, which removes all realizations at "
                                               "alpha_red >= 0.3 when started before BBN (section F)"}

# ------------------------------------------------------------------------------------------
# F. Small Monte Carlo of Model 1 (Ahmed-Dodelson-Greene-Sorkin 2004; Das-Nasiri-Yazdi 2023)
# ------------------------------------------------------------------------------------------
om_h2 = Om * h ** 2
or_h2 = Or * h ** 2
fm = om_h2 / (om_h2 + or_h2)
fr = or_h2 / (om_h2 + or_h2)
H_m0_kms = 100 * math.sqrt(om_h2 + or_h2)               # Hubble rate today if rho_L = 0


def model1(alpha_red, nreal, a_start=1e-10, per_efold=80, seed=1, agrid=None):
    """Model 1 in units rho / (M^2 H_m0^2), V H_m0^4. rho_L = S/V, dS = alpha sqrt(dV) xi.
    Stationary start: S_i ~ N(0, alpha^2 V_i) (as if the walk ran since a = 0).
    Grid: log-uniform in a (per_efold steps per e-fold) unless agrid (increasing a values, last = 1) is given.
    Termination: H^2 <= 0 at a grid point (Model 1 has no dynamics beyond that point).
    Returns arrays sampled at a_BBN = 1e-9, a_rec, a = 0.5 and a = 1."""
    g = np.random.default_rng(seed)
    if agrid is None:
        x = np.linspace(math.log(a_start), 0.0, int(round(-math.log(a_start) * per_efold)) + 1)
    else:
        x = np.log(np.asarray(agrid, float))
    av = math.exp(x[0])
    tau_ = np.full(nreal, av / math.sqrt(fr))
    # nested (Cauchy) form: V = 8 pi J3, J0 = int a^4 dtau, J_{n+1} = int J_n dtau  (no cancellation)
    J = [np.full(nreal, fr ** 2 * tau_[0] ** (5 + j) / den) for j, den in enumerate((5, 30, 210, 1680))]
    Vv = 8 * math.pi * J[3]
    S = alpha_red * np.sqrt(Vv) * g.standard_normal(nreal)
    rhoL_ = S / Vv
    alive = np.ones(nreal, bool)
    samples = {"bbn": 1e-9, "rec": 1 / 1090.8, "z1": 0.5, "today": 1.0}
    out = {kk: np.full(nreal, np.nan) for kk in samples}
    Vout = {kk: np.full(nreal, np.nan) for kk in samples}
    targets = {kk: (int(np.argmin(np.abs(x - math.log(v)))) if v >= av else -1) for kk, v in samples.items()}

    def rho_bg(aa_):
        return 3 * (fm * aa_ ** -3 + fr * aa_ ** -4)

    for n in range(len(x) - 1):
        a1, a2 = math.exp(x[n]), math.exp(x[n + 1])
        dx = x[n + 1] - x[n]
        tot1 = rho_bg(a1) + rhoL_
        tot2 = rho_bg(a2) + rhoL_
        bad = (tot1 <= 0) | (tot2 <= 0)
        alive &= ~bad
        H1 = np.sqrt(np.where(alive, tot1, 1.0) / 3)
        H2 = np.sqrt(np.where(alive, tot2, 1.0) / 3)
        f1 = np.where(alive, 1 / (a1 * H1), 0.0)       # dtau/dx; frozen once terminated
        f2 = np.where(alive, 1 / (a2 * H2), 0.0)
        tau_ = tau_ + 0.5 * dx * (f1 + f2)
        Jold = [jj.copy() for jj in J]
        J[0] = Jold[0] + 0.5 * dx * (a1 ** 4 * f1 + a2 ** 4 * f2)
        for j in range(1, 4):
            J[j] = Jold[j] + 0.5 * dx * (Jold[j - 1] * f1 + J[j - 1] * f2)
        Vn = 8 * math.pi * J[3]
        dV = np.clip(Vn - Vv, 0, None)
        S = S + alpha_red * np.sqrt(dV) * g.standard_normal(nreal)
        Vv = Vn
        rhoL_ = np.where(alive, S / Vv, rhoL_)
        for kk, idx in targets.items():
            if n + 1 == idx:
                tot = rho_bg(a2) + rhoL_
                out[kk] = np.where(alive & (tot > 0), rhoL_ / np.where(tot > 0, tot, 1), np.nan)
                Vout[kk] = np.where(alive, Vv, np.nan)
    tot_today = rho_bg(1.0) + rhoL_
    alive &= tot_today > 0
    H0_kms = np.where(alive, H_m0_kms * np.sqrt(np.where(tot_today > 0, tot_today, 1) / 3), np.nan)
    return dict(alive=alive, Om=out, V=Vout, H0=H0_kms, rhoL_today=np.where(alive, rhoL_, np.nan))


# F0: volume integrator against direct quadrature on the Lambda-free background (alpha -> 0)
lna_m = np.linspace(math.log(1e-12), 0.0, 200001)
am = np.exp(lna_m)
Em = np.sqrt(fm * am ** -3 + fr * am ** -4)
taum = integrate.cumulative_trapezoid(1 / (am * Em), lna_m, initial=0.0) + am[0] / math.sqrt(fr)
fV = lambda la: math.exp(la) ** 3 * (taum[-1] - np.interp(la, lna_m, taum)) ** 3 / np.interp(la, lna_m, Em)
V_ref = (4 * math.pi / 3) * integrate.quad(fV, math.log(1e-12), 0.0, limit=500, epsrel=1e-10)[0]
conv = {}
for pe in (20, 40, 80, 160):
    r_chk = model1(1e-12, 2, a_start=1e-12, per_efold=pe, seed=5)
    conv[pe] = float(np.nanmean(r_chk["V"]["today"]) / V_ref - 1)
R["F0_volume_convergence (V/V_ref - 1 by steps per e-fold)"] = conv
R["F0_Hm0^4 V_today_Lambda_free"] = V_ref
record("F0_model1_volume_integrator_vs_quadrature_80/efold", abs(conv[80]), 2e-3,
       "nested trapezoid V(a=1) against direct quadrature; error O(dx^2): 2.3%, 0.58%, 0.15%, 0.04%")


def summarize(r, label):
    alive = r["alive"]
    Om0 = r["Om"]["today"][alive]
    win = (r["Om"]["today"] > sel["today_window"][0]) & (r["Om"]["today"] < sel["today_window"][1])
    ok_rec = np.abs(r["Om"]["rec"]) < sel["rec_abs_lt"]
    ok_bbn = np.abs(r["Om"]["bbn"]) < sel["bbn_abs_lt"]
    joint = alive & win & ok_rec & ok_bbn
    q = (lambda arr: [float(v) for v in np.percentile(arr, [16, 50, 84])]) if alive.sum() > 20 else (lambda arr: None)
    d = {"label": label,
         "fraction_terminated_H2<=0_before_today": float(1 - alive.mean()),
         "survivors": int(alive.sum()),
         "Omega_L_today_16/50/84pct(survivors)": q(Om0),
         "P(Omega_L_today>0 | survive)": float(np.mean(Om0 > 0)) if alive.sum() else None,
         "P(Omega_L_today in [0.64,0.74]) (all realizations)": float(np.mean(alive & win)),
         "n_today_window": int((alive & win).sum()),
         "P(|Omega_L(rec)|<0.05) (all realizations)": float(np.mean(ok_rec)),
         "P(|Omega_L(BBN)|<0.07) (all realizations)": float(np.mean(ok_bbn)),
         "n_joint(survive & today window & rec & BBN)": int(joint.sum()),
         "H0_km/s/Mpc 16/50/84pct(survivors)": q(r["H0"][alive])}
    m = alive & np.isfinite(r["Om"]["z1"])
    if m.sum() > 200:
        d["corr(Omega_z1, Omega_0) empirical (survivors)"] = float(np.corrcoef(r["Om"]["z1"][m], r["Om"]["today"][m])[0, 1])
        d["sqrt(V_z1/V_0) (linear prediction)"] = float(np.nanmean(np.sqrt(r["V"]["z1"][m] / r["V"]["today"][m])))
    return d


sel = {"rec_abs_lt": 0.05, "today_window": (0.64, 0.74), "bbn_abs_lt": 0.07}
alphas = [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.72]
NREAL = 20000
mc = {}
for al in alphas:
    r = model1(al, NREAL, a_start=1e-10, per_efold=80, seed=SEED + int(al * 1000))
    d = summarize(r, f"alpha_red={al}")
    d["alpha_G1 (= alpha_red/8pi, Das et al. convention)"] = al / (8 * math.pi)
    mc[al] = d
    if al == 0.02:
        small_std = float(np.std(r["Om"]["today"][r["alive"]]))
R["F_model1_MC"] = {"settings": {"a_start": 1e-10, "steps_per_efold": 80, "n_realizations": NREAL,
                                 "start": "stationary S_i ~ N(0, alpha^2 V_i)", "Omega_m h^2": om_h2,
                                 "Omega_r h^2": or_h2, "H_m0 (rho_L = 0) km/s/Mpc": H_m0_kms,
                                 "selection (illustrative; cf. ZAS: Omega_de ~ 0.7 at z~0, ~0 at z~1000 at 5% level)": sel,
                                 "status": "computation of a published phenomenological model; all selections use "
                                           "already-exposed data (retrospective)"},
                    "by_alpha_red": {str(kk): vv for kk, vv in mc.items()}}
# small-alpha linear check on the Lambda-free background: std(Omega_0) = alpha / (3 sqrt(H_m0^4 V_ref))
pred_small = 0.02 / (3 * math.sqrt(V_ref))
R["F_small_alpha_check"] = {"MC std Omega today (alpha_red 0.02)": small_std, "linear prediction": pred_small}
record("F1_small_alpha_std_matches_linear", abs(small_std / pred_small - 1), 0.05,
       "20000 realizations; statistical error ~0.5%; nonlinearity of Omega small at alpha 0.02")
cz = mc[0.05]
record("F2_lag_correlation_matches_sqrt(V1/V2)",
       abs(cz["corr(Omega_z1, Omega_0) empirical (survivors)"] - cz["sqrt(V_z1/V_0) (linear prediction)"]), 0.03,
       "random-walk prediction Corr(rho(a1), rho(a2)) = sqrt(V1/V2); 20000 samples, se ~0.007")
record_bool("F3_survivor_bias_to_positive_Lambda",
            mc[0.2]["P(Omega_L_today>0 | survive)"] > mc[0.02]["P(Omega_L_today>0 | survive)"] and
            mc[0.2]["P(Omega_L_today>0 | survive)"] > 0.5,
            "termination removes negative excursions; cf. Ahmed-Sorkin 2012: 50-60% positive")
# step-size convergence of the termination fraction
r40 = model1(0.15, 8000, per_efold=40, seed=77)
r80 = model1(0.15, 8000, per_efold=80, seed=77)
r160 = model1(0.15, 8000, per_efold=160, seed=77)
R["F_convergence_alpha0.15"] = {"terminated_40/efold": float(1 - r40["alive"].mean()),
                                "terminated_80/efold": float(1 - r80["alive"].mean()),
                                "terminated_160/efold": float(1 - r160["alive"].mean())}
record("F4_conv_termination_fraction_80_vs_160", abs(r80["alive"].mean() - r160["alive"].mean()), 0.03,
       "termination fraction stable to 3 points (different noise paths; se ~0.6 points)")
# start-epoch scan (5000 realizations each)
scan = {}
for a_st in (1e-10, 1e-6, 1e-4, 1e-3):
    for al in (0.15, 0.2, 0.25, 0.3, 0.5, 0.72):
        r = model1(al, 5000, a_start=a_st, per_efold=80, seed=int(1e6 * a_st) + int(al * 100) + 3)
        alive = r["alive"]
        win = alive & (r["Om"]["today"] > 0.64) & (r["Om"]["today"] < 0.74)
        scan[f"a_start={a_st:g}, alpha_red={al}"] = {"survive": float(alive.mean()),
                                                      "P(today in [0.64,0.74])": float(win.mean())}
R["F_start_epoch_scan"] = scan
# coarse equal-steps-in-a recipe (Das et al. 2023 sec. 4.2 uses a_i = i/n): the termination fraction depends
# on the grid because H^2 is only checked at grid points.  Same stationary start at a = 1e-4.
coarse = {}
for al in (0.2, 0.25, 0.3):
    agr = np.arange(1, 10001) / 10000.0
    rc = model1(al, 5000, agrid=agr, seed=900 + int(al * 100))
    rf = model1(al, 5000, a_start=1e-4, per_efold=80, seed=900 + int(al * 100))
    coarse[f"alpha_red={al}"] = {"survive_equal_steps_in_a(n=10^4, from a=1e-4)": float(rc["alive"].mean()),
                                 "survive_log_grid_80/efold(from a=1e-4)": float(rf["alive"].mean())}
R["F_grid_dependence"] = coarse

# ------------------------------------------------------------------------------------------
# G. Barrow-type CMB estimate (order of magnitude)
# ------------------------------------------------------------------------------------------
# uniform sphere of radius R = 1/H with density contrast delta_rho: Phi_centre = -2 pi G delta_rho R^2
# = -delta_rho R^2/(4 M^2) = -(3/4) dOmega ; Sachs-Wolfe dT/T = Phi/3 = -dOmega/4.
dT_T_bound = 1e-5
dOmega_max = 4 * dT_T_bound
sig_rec = ep["recombination z = 1089.8"]["sigma_Omega_Lambda_per_alpha_red"]
R["G_Barrow_type"] = {"dT/T per unit dOmega (horizon patch, SW)": 0.25, "dT/T bound used": dT_T_bound,
                      "dOmega_rec max": dOmega_max, "sigma_Omega_per_alpha at rec": sig_rec,
                      "alpha_red max (independent horizon patches)": dOmega_max / sig_rec,
                      "alpha_red needed for 1 sigma today": ar_star,
                      "ratio needed/allowed": ar_star / (dOmega_max / sig_rec),
                      "status": "estimate [E]; conditional on spatially independent Lambda in disjoint past cones, "
                                "which Zwane-Afshordi-Sorkin dispute; no covariant inhomogeneous model exists"}

# ------------------------------------------------------------------------------------------
# H. Heavy thresholds
# ------------------------------------------------------------------------------------------
Lt, lbar, dl, Dm = sp.symbols("Lambda_tilde lambda_bar delta_lambda Delta_m", real=True)
Leff = Lt + lbar + dl
# prescription (i): <lambda> = 0 for the multiplier -> <Leff> = Lt, shifts by Delta under a threshold
pi_shift = sp.simplify((Leff.subs(Lt, Lt + Dm) - Leff).subs(dl, 0))
# prescription (ii): <Leff> = 0 imposed on the total -> lbar = -Lt is re-solved after each threshold
lbar_ii = sp.solve(sp.Eq(Lt + lbar, 0), lbar)[0]
pii_shift = sp.simplify((Lt + Dm + sp.solve(sp.Eq(Lt + Dm + lbar, 0), lbar)[0]) - (Lt + lbar_ii))
record_bool("H1_prescription_i_shifts_by_threshold", pi_shift == Dm)
record_bool("H2_prescription_ii_stable_because_total_is_stated", pii_shift == 0,
            "stable only because the total mean is imposed: a value statement (obstruction scorecard R1)")
me = 0.51099895e-3
rho_L_GeV4 = x_obs * M_red ** 4
R["H_threshold_example"] = {"electron m^4/(64 pi^2) / rho_L_obs": me ** 4 / (64 * math.pi ** 2) / rho_L_GeV4,
                            "spread delta_lambda": "independent of Lambda_tilde: at fixed V the factor e^{-i Lambda_tilde V} "
                                                   "is a global phase (Das et al. 2023 eq. 2.14)"}

# ------------------------------------------------------------------------------------------
R["checks"] = checks
R["n_checks"] = len(checks)
R["n_pass"] = sum(1 for c in checks if c["pass"])
with open(OUT, "w") as fh:
    json.dump(R, fh, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating,)) else str(o))
print(json.dumps({k: R[k] for k in ["D1_H4V", "D2_today", "D3_magnitudes", "D4_proxies", "D5_discreteness_length_for_1sigma",
                                    "E1_lag_correlation_sqrt(V1/V2)", "G_Barrow_type", "F_convergence_alpha0.15",
                                    "F_small_alpha_check", "F_grid_dependence", "C2_S_onshell_zero_gives_rho_L", "B4_Baum_Hawking"]},
                 indent=1, default=str))
for kk, vv in R["F_model1_MC"]["by_alpha_red"].items():
    print(kk, json.dumps(vv, default=str))
print(json.dumps(R["D2_epochs_LCDM"], indent=1))
print(f"checks: {R['n_pass']}/{R['n_checks']} pass")
for c in checks:
    if not c["pass"]:
        print("FAIL", c)
