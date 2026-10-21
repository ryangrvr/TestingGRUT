"""Referee checks for the vacuum-principle pass (track review_scope).

Deterministic; numpy/scipy only; no Monte Carlo; about 10 s.
Run: python3 -I code/review_checks.py results/review_checks.json

A. WDW minisuperspace of the `obstruction` track (k=+1, radiation E, Lambda):
   -(1/2) psi'' + (x^2/2 - lam x^4/2) psi = E psi,  psi(0)=0 or psi'(0)=0.
   Question: for lam > 0, are solutions square integrable at large x (Weyl limit-circle)?
B. Everpresent Lambda, Model 1 (rho = S/V, S Brownian in V), linear-Gaussian estimate of the
   statistic pre-registered by the `grut_native` track (bin-to-bin structure of rho_DE(z)),
   evaluated against already-exposed DESI DR2 w0wa central values. Retrospective [E].
C. Arithmetic re-checks of headline numbers quoted by the tracks.
"""
import json
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import erf

R = {}
checks = []


def record(name, value, passed):
    checks.append({"name": name, "value": value, "pass": bool(passed)})


# --------------------------------------------------------------------------------------------
# A. WDW endpoint at large x for lam > 0
# --------------------------------------------------------------------------------------------
def norm_growth(lam, E, psi0, dpsi0, xs_out, max_step=0.002):
    def rhs(x, y):
        V = 0.5 * x * x - 0.5 * lam * x ** 4
        return [y[1], 2.0 * (V - E) * y[0], y[0] ** 2]
    sol = solve_ivp(rhs, (0.0, xs_out[-1]), [psi0, dpsi0, 0.0], method="DOP853",
                    rtol=1e-10, atol=1e-12, max_step=max_step, dense_output=True)
    out = []
    for X in xs_out:
        y = sol.sol(X)
        out.append(float(y[2]))
    return out


XS = [10.0, 15.0, 20.0, 25.0, 30.0]
A = {}
for E in [1.0, 3.0, -20.0]:
    for lab, (p0, dp0) in {"psi(0)=0": (0.0, 1.0), "psi'(0)=0": (1.0, 0.0)}.items():
        N = norm_growth(0.05, E, p0, dp0, XS)
        # tail behaviour: for |psi|^2 ~ C/x^2 the increments N(X2)-N(X1) ~ C(1/X1 - 1/X2)
        inc = [N[i + 1] - N[i] for i in range(len(N) - 1)]
        pred_ratio = (1 / 25 - 1 / 30) / (1 / 20 - 1 / 25)  # expected inc[3]/inc[2]
        Ctail = inc[3] / (1 / 25 - 1 / 30)          # |psi|^2 ~ Ctail / x^2 on average
        N_inf = N[4] + Ctail / 30.0                  # extrapolated total norm
        A[f"lam=0.05,E={E},{lab}"] = {"N(X)": dict(zip([str(x) for x in XS], N)),
                                       "increments": inc,
                                       "inc_ratio_last_two": inc[3] / inc[2],
                                       "inc_ratio_if_x^-2_tail": pred_ratio,
                                       "N_extrapolated_infinity": N_inf,
                                       "N_inf_over_N30": N_inf / N[4]}
        record(f"A1 lam=0.05 E={E} {lab}: norm increments follow the x^-2 tail (ratio 2/3), so N(inf) is finite",
               inc[3] / inc[2], abs(inc[3] / inc[2] - pred_ratio) < 0.02)
# contrast: Lambda < 0 and Lambda = 0, non-eigenvalue E -> exponential growth (limit point)
for lam, E in [(-0.1, 2.0), (0.0, 2.0)]:
    N = norm_growth(lam, E, 0.0, 1.0, [4.0, 6.0, 8.0], max_step=0.002)
    A[f"lam={lam},E={E},psi(0)=0"] = {"N(X)": dict(zip(["4", "6", "8"], N))}
    record(f"A2 lam={lam} E={E}: N(8)/N(4) huge (not square integrable)", N[2] / N[0], N[2] / N[0] > 1e6)
R["A_wdw_large_x"] = A
R["A_note"] = ("For lam>0 both independent solutions are square integrable at large x for every E "
               "(Weyl limit-circle): no normalizability quantization exists until a boundary condition "
               "at large x (a self-adjoint extension parameter) is supplied; with one, the spectrum is "
               "discrete and unbounded below. The Dirichlet box of the obstruction track is one such "
               "supplied condition, which is why its levels move with the box.")

# --------------------------------------------------------------------------------------------
# B. Everpresent Lambda: pre-registered bin statistic vs exposed data (retrospective estimate)
# --------------------------------------------------------------------------------------------
h, Om, Orad, OL = 0.6736, 0.3153, 9.2e-5, 0.6846  # Planck 2018 base LCDM as in the unimodular track
Ok = 1.0 - Om - Orad - OL


def Hof(a):
    return math.sqrt(Orad / a ** 4 + Om / a ** 3 + Ok / a ** 2 + OL)


def tau_of(a):  # conformal time, H0 = 1
    return quad(lambda x: 1.0 / (x * x * Hof(x)), 0.0, a, limit=400, epsabs=0, epsrel=1e-11)[0]


def V_cone(a0):  # past-light-cone 4-volume, H0^-4 units: (4 pi/3) int a^4 (tau0 - tau)^3 dtau
    t0 = tau_of(a0)
    f = lambda a: a ** 2 / Hof(a) * (t0 - tau_of(a)) ** 3
    # split for accuracy
    pts = [0.0, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.3, a0]
    pts = [p for p in pts if p <= a0]
    tot = sum(quad(f, pts[i], pts[i + 1], limit=200, epsrel=1e-9)[0] for i in range(len(pts) - 1))
    return 4.0 * math.pi / 3.0 * tot


zs = [0.0, 0.5, 1.0, 2.0]
Bv = {}
for z in zs:
    a = 1.0 / (1.0 + z)
    V = V_cone(a)
    E2 = Hof(a) ** 2
    H4V = E2 ** 2 * V
    Bv[str(z)] = {"V_H0^4": V, "H^4V": H4V, "E^2": E2,
                  "sigmaOmega_over_alpha": 1.0 / (3.0 * math.sqrt(H4V))}
R["B_past_cone"] = Bv
record("B0 H0^4 V_cone(today) reproduces unimodular 0.1229", Bv["0.0"]["H^4V"],
       abs(Bv["0.0"]["H^4V"] - 0.1229) < 5e-4)
corr_z1 = math.sqrt(Bv["1.0"]["V_H0^4"] / Bv["0.0"]["V_H0^4"])
R["B_corr_z1_z0"] = corr_z1
record("B0b corr(z=1, z=0) = sqrt(V1/V0) reproduces 0.19", corr_z1, abs(corr_z1 - 0.19) < 0.01)


def sigma_rho_DE0(alpha, z):  # std of rho_Lambda(z) in units of rho_DE0 = OL rho_crit0
    b = Bv[str(z)]
    return alpha * b["sigmaOmega_over_alpha"] * b["E^2"] / OL


# exposed central values (DESI DR2 Table 5, w0wa, DESI+CMB+SN; quoted via the others track)
cpl = {"LCDM": (-1.0, 0.0), "DESI+CMB+Pantheon+": (-0.838, -0.62),
       "DESI+CMB+Union3": (-0.667, -1.09), "DESI+CMB+DESY5": (-0.752, -0.86)}
expo = {}
for k, (w0, wa) in cpl.items():
    expo[k] = {str(z): (1 + z) ** (3 * (1 + w0 + wa)) * math.exp(-3 * wa * z / (1 + z)) for z in zs[1:]}
R["B_exposed_rhoDE_over_rhoDE0_central"] = expo


def p_bins(alpha, delta, bins=(0.5, 1.0, 2.0), delta_last=None):
    """P that every backward increment |rho(z_k) - rho(z_{k-1})| <= delta (rho_DE0 units).
    Model 1: given rho at larger V, rho at smaller V has mean equal to it and variance
    sigma_k^2 - sigma_{k-1}^2 (Brownian-bridge Markov property)."""
    P = 1.0
    prev = 0.0
    for z in bins:
        s2 = sigma_rho_DE0(alpha, z) ** 2 - sigma_rho_DE0(alpha, prev) ** 2
        s = math.sqrt(s2)
        d = delta_last if (delta_last is not None and z == bins[-1]) else delta
        P *= erf(d / (s * math.sqrt(2.0)))
        prev = z
    return P


def p_today(alpha, lo=0.64, hi=0.74):
    s = alpha * Bv["0.0"]["sigmaOmega_over_alpha"]
    return 0.5 * (erf(hi / (s * math.sqrt(2))) - erf(lo / (s * math.sqrt(2))))


alpha_star = 3 * OL * math.sqrt(Bv["0.0"]["H^4V"])
R["B_alpha_star_1sigma_today"] = alpha_star
record("B1 alpha_red* = 3 Omega_L sqrt(H0^4 V) reproduces 0.72", alpha_star, abs(alpha_star - 0.72) < 0.005)
Bres = {}
for delta in [0.2, 0.3, 0.5]:
    pb_star = p_bins(alpha_star, delta)
    pb_star_2 = p_bins(alpha_star, delta, bins=(0.5, 1.0))
    grid = np.exp(np.linspace(math.log(0.005), math.log(5.0), 2000))
    joint = np.array([p_today(a) * p_bins(a, delta) for a in grid])
    joint2 = np.array([p_today(a) * p_bins(a, delta, bins=(0.5, 1.0)) for a in grid])
    joint3 = np.array([p_today(a) * p_bins(a, delta, delta_last=1.0) for a in grid])
    i = int(np.argmax(joint)); j = int(np.argmax(joint2)); k3 = int(np.argmax(joint3))
    Bres[str(delta)] = {"P_bins(alpha*) z=0.5,1,2": pb_star, "P_bins(alpha*) z=0.5,1": pb_star_2,
                        "max_alpha P_today*P_bins (z=0.5,1,2)": float(joint[i]), "at_alpha": float(grid[i]),
                        "max_alpha P_today*P_bins (z=0.5,1)": float(joint2[j]), "at_alpha_2bins": float(grid[j]),
                        "P_bins(alpha*) z=0.5,1 + z=2 at delta=1.0": p_bins(alpha_star, delta, delta_last=1.0),
                        "max_alpha P_today*P_bins (z=0.5,1 + z=2 at delta=1.0)": float(joint3[k3])}
R["B_prereg_statistic_vs_exposed"] = Bres
record("B2 at alpha* and delta=0.3 rho_DE0, P(bins z<=2 near-constant) < 1% (pre-registered failure level)",
       Bres["0.3"]["P_bins(alpha*) z=0.5,1,2"], Bres["0.3"]["P_bins(alpha*) z=0.5,1,2"] < 0.01)
record("B3 even maximized over alpha, P(today window AND bins z<=2) < 1% for delta = 0.3",
       Bres["0.3"]["max_alpha P_today*P_bins (z=0.5,1,2)"], Bres["0.3"]["max_alpha P_today*P_bins (z=0.5,1,2)"] < 0.01)
record("B4 robust form: z<=1 bins only, delta = 0.5 rho_DE0, maximized over alpha, joint with today window < 1e-3",
       Bres["0.5"]["max_alpha P_today*P_bins (z=0.5,1)"], Bres["0.5"]["max_alpha P_today*P_bins (z=0.5,1)"] < 1e-3)
R["B_note"] = ("Linear-Gaussian Model 1 on the Planck-2018 LCDM background; ignores termination, BBN and "
               "CMB (each would only lower the probabilities) and back-reaction of rho_Lambda on V. delta is "
               "a deliberately generous tolerance on bin-to-bin changes of rho_DE in units of rho_DE0; the "
               "exposed DESI DR2 w0wa central values sit within a few to ~20 percent of constancy for z<=2.")

# --------------------------------------------------------------------------------------------
# C. Arithmetic re-checks
# --------------------------------------------------------------------------------------------
Mred = 2.435e18          # GeV
MP = 1.22089e19          # GeV, non-reduced
H0 = 2.1332e-42 * 0.674  # GeV
rhoL = 3 * 0.685 * (H0 / Mred) ** 2  # in M^4 units
C = {}
C["rhoL_over_M4"] = rhoL
C["rhoL_GeV4"] = rhoL * Mred ** 4
C["rhoL^(1/4)_meV"] = (rhoL * Mred ** 4) ** 0.25 * 1e12
VM4 = 1.01e240
C["K1_p_reduced_pastcone"] = -math.log(7.15e-121) / math.log(VM4)
me = 0.51099895e-3
C["electron_one_dof_m4_over_64pi2_over_rhoL"] = me ** 4 / (64 * math.pi ** 2) / (rhoL * Mred ** 4)
# Barrow-Shaw Omega_k0 vs exposed Omega_K summaries (sigma distances)
bs = {"printed": -0.0055, "planck_like": -0.00637, "desi_like": -0.00619}
ok = {"DESI DR2+CMB LCDM+OmegaK (Table 5)": (0.0023, 0.0011),
      "DESI DR2+CMB+Pantheon+ w0wa+OmegaK (Table 5)": (0.0011, 0.0013),
      "Planck 2018 TT,TE,EE+lowE+lensing+BAO (recall)": (0.0007, 0.0019)}
C["BarrowShaw_tension_sigma"] = {f"{kb} vs {ko}": (mo - vb) / so
                                for kb, vb in bs.items() for ko, (mo, so) in ok.items()}
# Wang-Zhu-Unruh: exponent needed for H = Lambda_UV exp(-x) with Lambda_UV = M_P
C["WZU_exponent_needed_ln(MP/H0)"] = math.log(MP / H0)
# Arkani-Hamed-Hall-Kolda-Murayama scale (1 TeV)^2/M, compare rho_L^(1/4)
C["TeV^2/M_red_meV"] = (1e3) ** 2 / Mred * 1e12
C["TeV^2/M_P_meV"] = (1e3) ** 2 / MP * 1e12
R["C_arithmetic"] = C
record("C1 rho_L^(1/4) = 2.24 meV", C["rhoL^(1/4)_meV"], abs(C["rhoL^(1/4)_meV"] - 2.24) < 0.01)
record("C2 K1 p = 0.5006", C["K1_p_reduced_pastcone"], abs(C["K1_p_reduced_pastcone"] - 0.5006) < 2e-4)
record("C3 electron one-dof threshold 4.3e30 rho_L", C["electron_one_dof_m4_over_64pi2_over_rhoL"],
       abs(C["electron_one_dof_m4_over_64pi2_over_rhoL"] / 4.3e30 - 1) < 0.02)

R["checks"] = checks
R["n_checks"] = len(checks)
R["n_pass"] = sum(c["pass"] for c in checks)
R["versions"] = {"numpy": np.__version__, "scipy": __import__("scipy").__version__,
                 "python": sys.version.split()[0]}
out = sys.argv[1] if len(sys.argv) > 1 else "review_checks.json"
with open(out, "w", encoding="utf-8") as fh:
    json.dump(R, fh, indent=1)
print(f"{R['n_pass']}/{R['n_checks']} checks pass")
for c in checks:
    print(("PASS " if c["pass"] else "FAIL ") + c["name"], c["value"] if not isinstance(c["value"], dict) else "")
