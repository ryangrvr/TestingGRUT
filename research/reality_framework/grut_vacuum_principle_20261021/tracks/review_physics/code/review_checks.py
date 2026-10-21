"""Adversarial physics review: independent checks (review_physics track).

Sections
  S1 conventions and scale arithmetic (rho_Lambda/M^4, thresholds, N-bound, graviton-loop cutoff)
  S2 symbolic: Delta T for -C int sqrt(-g); trace; sequestering residual; trace-free algebra;
     Tseytlin / Smolin on-shell residuals; Gibbs uniqueness is not re-done (elementary)
  S3 Wheeler-DeWitt minisuperspace with Lambda > 0: limit-circle test at x -> infinity
     (are E-solutions normalizable?  does a box impose a self-adjoint extension?)
  S4 unimodular / trace-free FRW with a vacuum-energy phase transition (time-dependent rho_vac)
  S5 Lombriser mirror-collapse toy (sign and magnitude)
  S6 past-light-cone 4-volume, alpha_red, sigma(Omega) per epoch
  S7 others-track arithmetic (ADDG, DHK, CKN, Weinberg 500/729, BP ratio, Barrow-Shaw tension)
  S8 grut_native checks (closed-dust V4, lambda_E, flat dust Lambda<0 ratio)
Units hbar = c = 1; M = reduced Planck mass.
"""
import json
import math
import sys

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

R = {}
checks = []


def check(name, ok, **vals):
    checks.append({"name": name, "pass": bool(ok), **{k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in vals.items()}})


# ------------------------------------------------------------------ S1 scales
M = 2.435e18                      # GeV
MP = math.sqrt(8 * math.pi) * M
h = 0.674
H0 = 2.1332e-42 * h               # GeV  (100 km/s/Mpc = 2.1332e-42 GeV)
OL = 0.685
rhoL = 3 * OL * H0**2 * M**2
R["S1"] = {
    "H0_GeV": H0, "(H0/M)^2": (H0 / M) ** 2, "rhoL_GeV4": rhoL, "rhoL^(1/4)_meV": rhoL**0.25 * 1e12,
    "rhoL/M^4": rhoL / M**4, "rhoL/MP^4": rhoL / MP**4,
    "ln d = 24 pi^2 M^4/rhoL": 24 * math.pi**2 * M**4 / rhoL,
}
for nm, m in (("electron", 0.511e-3), ("higgs", 125.0), ("top", 173.0)):
    R["S1"]["m^4/64pi^2 / rhoL (" + nm + ")"] = m**4 / (64 * math.pi**2) / rhoL
Mc = (rhoL * M**2) ** (1 / 6)
Mc_loop = (rhoL * M**2 * (16 * math.pi**2) ** 2) ** (1 / 6)
R["S1"]["graviton-loop cutoff M_c (M_c^6/M^2 = rhoL) GeV"] = Mc
R["S1"]["same with (16pi^2)^2 GeV"] = Mc_loop
R["S1"]["(1 TeV)^6/M^2 / rhoL"] = 1e18 / M**2 / rhoL
check("S1 rhoL/M^4 = 7.16e-121 (obstruction F7)", abs(rhoL / M**4 / 7.16e-121 - 1) < 0.01, val=rhoL / M**4)
check("S1 graviton-loop cutoff 23 MeV / 0.12 GeV (sequestering)", abs(Mc - 0.023) < 0.001 and abs(Mc_loop - 0.12) < 0.01,
      Mc=Mc, Mc_loop=Mc_loop)

# ------------------------------------------------------------------ S2 symbolic
g = sp.diag(-1, 1, 1, 1) * sp.Symbol("A", positive=True) ** 2  # conformally flat point value is enough
C = sp.Symbol("C")
# T_mn = -(2/sqrt(-g)) dS/dg^mn ; d sqrt(-g)/d g^mn = -(1/2) sqrt(-g) g_mn
# S = -C int sqrt(-g):  dS/dg^mn = -C * (-(1/2) sqrt(-g) g_mn) = (C/2) sqrt(-g) g_mn
T_shift = -2 * (C / 2) * g
R["S2"] = {"Delta_T for S_m -> S_m - C int sqrt(-g)": "T_mn -> T_mn - C g_mn (rho = +C)",
           "check": str(sp.simplify(T_shift + C * g) == sp.zeros(4, 4))}
check("S2 Delta T = -C g for -C int sqrt(-g)", sp.simplify(T_shift + C * g) == sp.zeros(4, 4))
rho, p = sp.symbols("rho p")
ginv = sp.diag(-1, 1, 1, 1)
Tud = sp.diag(-rho, p, p, p)       # T^mu_nu
check("S2 trace T = -rho + 3p in (-,+,+,+)", sp.simplify(Tud.trace() - (-rho + 3 * p)) == 0)
# sequestering residual: M^2 G^mu_nu = T^mu_nu - (1/4) delta <T>;  T = -V delta + tau
V, tauavg = sp.symbols("V tauavg")
rhs = -V + sp.Symbol("tau00") - sp.Rational(1, 4) * (-4 * V + tauavg)
check("S2 sequestering: V_vac cancels, Lambda_eff = <tau>/4",
      sp.simplify(rhs - (sp.Symbol("tau00") - tauavg / 4)) == 0)
# Tseytlin S/V: M^2 G = T - <M^2 R/2 + L_m> g; trace M^2 R = 4 Leff - tau ; Leff = <M^2R/2 + L_tau>
Leff, tav, Lt = sp.symbols("Leff tav Lt")
sol = sp.solve(sp.Eq(Leff, (4 * Leff - tav) / 2 + Lt), Leff)[0]
R["S2"]["Tseytlin Lambda_eff"] = str(sol)
check("S2 Tseytlin Lambda_eff = <tau>/2 - <L_tau>", sp.simplify(sol - (tav / 2 - Lt)) == 0)

# ------------------------------------------------------------------ S3 WDW Lambda > 0
# -1/2 psi'' + (x^2/2 - lam x^4/2) psi = E psi,  psi(0) = 0
# Weyl limit-point/limit-circle at infinity: for -d^2/dx^2 + q with q ~ -x^4 (power > 2) the
# endpoint is limit circle: EVERY solution is square integrable near infinity (WKB |psi|^2 ~ x^-2).
lam = 0.05


def wdw_norm(E, X):
    f = lambda x, y: [y[1], 2 * (0.5 * x**2 - 0.5 * lam * x**4 - E) * y[0], y[0] ** 2]
    s = solve_ivp(f, (0, X), [0.0, 1.0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
    return s


s3 = {}
for E in (-20.0, 0.0, 3.7, 25.0):
    sol = wdw_norm(E, 24.0)
    Ns = [float(sol.sol(X)[2]) for X in (8.0, 12.0, 16.0, 20.0, 24.0)]
    # WKB envelope: psi ~ A x^-1 -> tail of int psi^2 from X to inf ~ A^2 / X
    s3[str(E)] = {"N(X) for X=8,12,16,20,24": Ns,
                  "N(24)-N(16)": Ns[4] - Ns[2], "N(16)-N(8)": Ns[2] - Ns[0]}
R["S3"] = {"lam": lam, "norm_growth": s3}
# a norm-convergent tail ~ 1/X gives ratio (1/16-1/24)/(1/8-1/16) = 1/3; a non-normalizable
# (continuum-like) solution would give a ratio >= 1
for v in s3.values():
    v["increment_ratio"] = v["N(24)-N(16)"] / v["N(16)-N(8)"]
conv = all(v["increment_ratio"] < 0.6 for v in s3.values())
check("S3 Lambda>0: solutions for every E have converging norm (limit circle; not 'no normalizable states')", conv)


# box levels: Dirichlet at L is one self-adjoint-extension family; spectrum discrete, moves with L
def box_levels(L, n=4000, nlev=3):
    from scipy.linalg import eigh_tridiagonal
    xs = np.linspace(0, L, n + 2)[1:-1]
    hh = xs[1] - xs[0]
    Vx = 0.5 * xs**2 - 0.5 * lam * xs**4
    w = eigh_tridiagonal(1.0 / hh**2 + Vx, -0.5 / hh**2 * np.ones(n - 1), select="i", select_range=(0, nlev - 1))[0]
    return [float(v) for v in w]


R["S3"]["box_levels"] = {str(L): box_levels(L) for L in (12.0, 16.0)}
R["S3"]["V(L)"] = {str(L): 0.5 * L**2 - 0.5 * lam * L**4 for L in (12.0, 16.0)}

# ------------------------------------------------------------------ S4 trace-free FRW with phase transition
# trace-free FRW: Hdot - k/a^2 = -(rho + p)/(2 M^2) (vacuum drops out); total energy conserved.
# rho_vac(t) = V0 + Delta*(1 - S(t)) with S a smooth step; latent heat -> radiation.
# Diagnostic: X(t) = 3 M^2 (H^2 + k/a^2) - rho_nonvac  (the 'gravitating constant' read on a slice)
M2 = 1.0
Delta, V0, tstar, w = 5.0, 1000.0, 1.0, 0.02
Sstep = lambda t: 0.5 * (1 + math.tanh((t - tstar) / w))
dS = lambda t: 0.5 / w / math.cosh((t - tstar) / w) ** 2


def rhs_tf(t, y):
    a, H, rr = y
    p = rr / 3                     # vacuum has rho + p = 0
    dH = -(rr + p) / (2 * M2)
    drr = -4 * H * rr + Delta * dS(t)   # latent heat from vacuum into radiation
    return [a * H, dH, drr]


rr0, a0 = 30.0, 1.0
Lhat = 2.0                          # integration constant chosen by initial data
# set initial H so that 3H^2 = rr0 + rho_eff0 with rho_eff0 = Lhat + rho_vac(0) = Lhat + V0 + Delta
H0tf = math.sqrt((rr0 + Lhat + V0 + Delta) / 3)
sol = solve_ivp(rhs_tf, (0, 2.0), [a0, H0tf, rr0], method="DOP853", rtol=1e-12, atol=1e-12, dense_output=True)
X = lambda t: 3 * M2 * sol.sol(t)[1] ** 2 - sol.sol(t)[2]
R["S4"] = {"rho_eff before (t=0.5)": X(0.5), "rho_eff after (t=1.5)": X(1.5),
           "change": X(1.5) - X(0.5), "Delta_vac": -Delta,
           "note": "trace-free equations: rho_eff = Lambda_hat + rho_vac(t); a vacuum change gravitates as in GR"}
check("S4 unimodular: phase transition changes rho_eff by -Delta (no protection for changes)",
      abs((X(1.5) - X(0.5)) + Delta) < 1e-6, change=X(1.5) - X(0.5))

# ------------------------------------------------------------------ S5 Lombriser mirror toy
def lomb(om, orad, ac):
    ol = 1 - om - orad
    E = lambda a: math.sqrt(orad / a**4 + om / a**3 + ol)
    t = quad(lambda a: 1 / (a * E(a)), 0, ac, limit=200)[0]
    v = quad(lambda a: a**2 / E(a), 0, ac, limit=200)[0]
    return -0.25 * om * t / v / ol


R["S5"] = {"ratio a_c=1 (Om=0.315)": lomb(0.315, 9.2e-5, 1.0),
           "ratio a_c=1 (Om=0.308)": lomb(0.308, 9.2e-5, 1.0),
           "a_c for |ratio|=1 (0.315)": brentq(lambda a: abs(lomb(0.315, 9.2e-5, a)) - 1, 0.3, 1.0),
           "a_c for |ratio|=1 (0.308)": brentq(lambda a: abs(lomb(0.308, 9.2e-5, a)) - 1, 0.3, 1.0)}
check("S5 Lombriser toy -0.408 and a_c = 0.721", abs(R["S5"]["ratio a_c=1 (Om=0.315)"] + 0.408) < 0.003
      and abs(R["S5"]["a_c for |ratio|=1 (0.315)"] - 0.721) < 0.003)

# ------------------------------------------------------------------ S6 past light cone 4-volume
hP, omP = 0.6736, 0.3153
OrP = 9.2e-5
OLP = 1 - omP - OrP
EP = lambda a: math.sqrt(OrP / a**4 + omP / a**3 + OLP)


def conformal_time(a):
    return quad(lambda x: 1 / (x**2 * EP(x)), 0, a, limit=400, epsabs=0, epsrel=1e-12)[0]


def Vcone(a_obs):
    tau0 = conformal_time(a_obs)
    # V = (4 pi/3) int_0^tau0 a^4 (tau0 - tau)^3 dtau = (4pi/3) int_0^a_obs a^4 (tau0 - tau(a))^3 da/(a^2 E)
    f = lambda a: a**2 * (tau0 - conformal_time(a)) ** 3 / EP(a)
    # tabulate tau(a) for speed
    ag = np.logspace(-9, math.log10(a_obs), 4001)
    tg = np.array([conformal_time(x) for x in ag[::50]])
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(np.log(ag[::50]), tg)
    f = lambda la: math.exp(3 * la) * (tau0 - cs(la)) ** 3 / EP(math.exp(la))
    val = quad(f, math.log(1e-9), math.log(a_obs), limit=400, epsrel=1e-10)[0]
    return (4 * math.pi / 3) * val


V0c = Vcone(1.0)
alpha_red = 3 * OLP * math.sqrt(V0c)
R["S6"] = {"H0^4 V_cone(today)": V0c, "alpha_red_1sigma": alpha_red, "alpha_G1": alpha_red / (8 * math.pi),
           "sigma(Omega)/alpha today": 1 / (3 * math.sqrt(V0c)),
           "radiation H^4V exact 4pi/840": 4 * math.pi / 840, "matter H^4V exact": 16 * (4 * math.pi / 3) * math.factorial(8) * 6 / math.factorial(12)}
arec = 1 / 1090.8
Vr = Vcone(arec)
R["S6"]["H^4 V at recombination"] = Vr * EP(arec) ** 4
R["S6"]["sigma(Omega)/alpha rec"] = 1 / (3 * math.sqrt(Vr * EP(arec) ** 4))
check("S6 H0^4 V_cone = 0.1229, alpha_red = 0.720", abs(V0c - 0.1229) < 0.0005 and abs(alpha_red - 0.720) < 0.003,
      V=V0c, alpha=alpha_red)

# ------------------------------------------------------------------ S7 others arithmetic
E1 = 1e12
R["S7"] = {
    "ADDG Mbar (E=TeV^4) GeV": math.sqrt(E1 * M**2 / rhoL),
    "ADDG Mbar (E=M^4) GeV": math.sqrt(M**6 / rhoL),
    "ADDG 10 TeV factor": 1 + (1e4) ** 4 / (64 * math.pi**2) / E1,
    "DHK L*H0": 1 / math.sqrt(3 * OL),
    "CKN Lambda_UV reduced meV": math.sqrt(M * H0) * 1e12,
    "Weinberg rhoV_max/rho_obs (z_c=4)": (500 / 729) * (1.686 * 5) ** 3 * 0.315 / 0.685,
    "BP rhoL/rho_rad(1 MeV)": rhoL / (math.pi**2 / 30 * 10.75 * 1e-12),
    "Barrow-Shaw tension sigma (-0.00619 vs 0.0023+-0.0011)": (0.0023 + 0.00619) / 0.0011,
}
# Weinberg coefficient: Lambda/3 <= K^3/(27 (GM)^2), K = (5/3) d H^2 r^2, GM = H^2 r^3/2
d, Hs, r = sp.symbols("d H r", positive=True)
coef = sp.simplify(((sp.Rational(5, 3) * d * Hs**2 * r**2) ** 3 / (27 * (Hs**2 * r**3 / 2) ** 2)) / (d**3 * Hs**2))
R["S7"]["Weinberg coefficient"] = str(coef)
check("S7 Weinberg bound coefficient 500/729", coef == sp.Rational(500, 729))

# ------------------------------------------------------------------ S8 grut_native
# closed dust (3M^2=1, A=1): a(eta) = (1 - cos eta)/2, V4 = 2 pi^2 int a^4 deta
V4 = 2 * math.pi**2 * quad(lambda e: ((1 - math.cos(e)) / 2) ** 4, 0, 2 * math.pi)[0]
R["S8"] = {"closed dust V4": V4, "35 pi^3/32": 35 * math.pi**3 / 32}
check("S8 closed-dust V4 = 35 pi^3/32", abs(V4 - 35 * math.pi**3 / 32) < 1e-9)
# flat dust + Lambda<0: <rho_m>/|rho_L| = 2 (analytic: a^3 ∝ sin^2)
lamb = 0.7
T = 2 * math.pi / (3 * math.sqrt(lamb))
num = quad(lambda t: 1.0, 0, T)[0]                          # a^3 rho_m = A = const (A = 1)
den = quad(lambda t: math.sin(1.5 * math.sqrt(lamb) * t) ** 2 / lamb, 0, T)[0]
R["S8"]["flat dust L<0 ratio"] = num / den / lamb
check("S8 flat dust Lambda<0 <rho>/|rhoL| = 2", abs(num / den / lamb - 2) < 1e-9)

R["checks"] = checks
R["n_pass"] = sum(c["pass"] for c in checks)
R["n_checks"] = len(checks)
path = sys.argv[1] if len(sys.argv) > 1 else "review_checks.json"
with open(path, "w") as f:
    json.dump(R, f, indent=1, default=str)
print(json.dumps(R, indent=1, default=str))
