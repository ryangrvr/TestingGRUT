"""v_critic: small cross-module consistency checks and reference-history bookkeeping.

Own code (numpy/scipy only). It reads module RESULTS JSON only as quoted data (never module code):
  modules/charge_phase/results/thermal_phase.json  (n_local(T_h) column of the leading-order (T, n) phase diagram)
Everything else is closed-form arithmetic from BRIEF.md inputs. Labels in the output: [D] derived here,
[C] computed here, [Q] quoted, [E] order-of-magnitude estimate.

Run:  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 critic_checks.py   (< 5 s)
"""
import json
import math
import os

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.linalg import eigh_tridiagonal

HERE = os.path.dirname(os.path.abspath(__file__))
PKT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = {}

# ---------------------------------------------------------------- inputs (BRIEF.md) ----------------------------
v = 1000.0            # GeV
gX, lphi, lX, hD, lD = 0.1, 0.01, 0.001, 0.01, 0.001
MP = 2.435323e18      # reduced Planck mass, GeV
hbar_s = 6.582119569e-25   # GeV s
gstar = 106.75
sigma = 2 * math.sqrt(2 * lphi) / 3 * v ** 3      # GeV^3
mX = gX * v
mstar = v * min(gX, (lphi * lX) ** 0.25)
Y_B = 8.674e-11       # [Q] charge_phase / relic_yield
s_over_T3 = 2 * math.pi ** 2 * gstar / 45         # visible entropy density / T^3
OUT["inputs"] = dict(sigma_GeV3=sigma, mX=mX, mstar=mstar, Y_B=Y_B)


# ---------------------------------------------------------------- 1. T_dom conventions ---------------------------
def H_of_T(T, g):
    return math.sqrt(math.pi ** 2 * g / 90) * T ** 2 / MP


def T_of_t_RD(t_s, g):
    """radiation era t = 1/(2H)."""
    t = t_s / hbar_s
    H = 1 / (2 * t)
    return math.sqrt(H * MP / math.sqrt(math.pi ** 2 * g / 90))


g_keV = 3.383   # [Q] Saikawa-Shirai g_rho0 used by cosmo_observable
A = 0.8
t_w = 3 * MP ** 2 / (4 * A * sigma) * hbar_s
OUT["T_dom"] = {
    "t_w_noBR_s [D] (network_sim formula)": t_w,
    "T at t_w (no backreaction) keV [C]": T_of_t_RD(t_w, g_keV) * 1e6,
    "T for sigma/(M_P^2 H)=1 (inflation_spectator criterion) keV [C]":
        math.sqrt(sigma / MP ** 2 * MP / math.sqrt(math.pi ** 2 * g_keV / 90)) * 1e6,
    "T for 2 A sigma H = 3 M_P^2 H^2 (A=0.8, no BR) keV [C]":
        math.sqrt(2 * A * sigma / (3 * MP ** 2) * MP / math.sqrt(math.pi ** 2 * g_keV / 90)) * 1e6,
    "note": "cosmo_observable production value with backreaction is 6.367 keV at t = 2.71e4 s [Q]; t_w = 3.9e4 s "
            "is the no-backreaction time (5.82 keV); 8.0 keV is a one-wall-per-Hubble-patch criterion.",
}

# ---------------------------------------------------------------- 2. charge exit bookkeeping --------------------
tp = json.load(open(os.path.join(PKT, "modules", "charge_phase", "results", "thermal_phase.json")))
rows = [(r["T"], r["n_local"]) for r in tp["rows"] if r.get("n_local") is not None]
Tg = np.array([r[0] for r in rows])
nl = np.array([r[1] for r in rows])
OUT["quoted_n_local_table [Q charge_phase]"] = {"T_h/v": Tg.tolist(), "n_local/v^3": nl.tolist()}

Tc_res = (1.5571397, 1.5641468)   # [Q] thermal_potential adopted range
Tc0_lo = 1.5575                    # [Q] charge_phase leading order


def nQ_over_Th3(r, xi):
    """n_Q/T_h^3 with T_h = xi T_SM, both entropies separately conserved [D]."""
    return r * Y_B * s_over_T3 / xi ** 3


def exit_Th(r, xi):
    """T_h (units v) where n_Q/T_h^3 = n_local(T_h)/T_h^3, using the quoted table (T_h <= 1.545 grid; interpolation
    in T).  Returns None if the crossing lies above the last tabulated point (exit is then within ~0.01 v of T_c0,
    i.e. effectively thermal)."""
    c = nQ_over_Th3(r, xi)
    Tpos = Tg[Tg > 0]
    f = nl[Tg > 0] / Tpos ** 3 - c
    # f decreases with T (n_local/T^3 falls); find first sign change
    for i in range(len(Tpos) - 1):
        if f[i] >= 0 and f[i + 1] < 0:
            # linear in T
            T0, T1 = Tpos[i], Tpos[i + 1]
            return float(T0 + (T1 - T0) * f[i] / (f[i] - f[i + 1]))
    if f[-1] >= 0:
        return None   # beyond the table: exit just below T_c0
    return float(Tpos[0])


def T_n_cold(r):
    """relic_yield / charge_phase cold-sector T_n = [n_c/(r Y_B s/T^3)]^{1/3} (T in GeV) [D]."""
    nc = math.sqrt(8 * lX) * (lphi / (2 * gX ** 2)) ** 1.5 * v ** 3
    return (nc / (r * Y_B * s_over_T3)) ** (1 / 3)


ex = {}
for r in (1e-4, 1e-3, 3e-3, 1e-2, 0.05, 0.0897, 0.64, 1.0):
    for xi in (5e-4, 7.58e-4):
        Th = exit_Th(r, xi)
        key = f"r={r:g}, xi={xi:g}"
        rBEC = (xi / 7.09e-3) ** 3
        ex[key] = {
            "nQ/Th^3 [D]": nQ_over_Th3(r, xi),
            "T_h,exit/v [C] (None = within 0.01 v of T_c0: thermal exit)": Th,
            "T_SM,exit TeV [C]": (Th * v / xi / 1e3) if Th else None,
            "relic/charge_phase cold T_n TeV [D]": T_n_cold(r) / 1e3,
            "condensed at T_c? (r > r_BEC=(xi/7.09e-3)^3) [D]": bool(r > rBEC),
            "r_BEC": rBEC,
        }
OUT["charged_exit_bookkeeping"] = ex

# r threshold at xi_max for a distinctly charged exit (T_h,exit <= 1.0 v and <= 0.75 v)
for Tcut in (1.0, 0.75):
    j = np.argmin(abs(Tg - Tcut))
    need = nl[j] / Tg[j] ** 3
    OUT[f"r needed at xi=7.58e-4 for exit at T_h <= {Tcut} v [C]"] = need / nQ_over_Th3(1.0, 7.58e-4)
# xi boundary for 'charge holds phi=0 continuously from T_c' using the maximum of n_local/T^3 just below T_c0
# (relic_yield used the cold n_c/T_c^3):
cold_ratio = 0.0316228 / 1.5476 ** 3
for Tref in (1.5, 1.3):
    j = np.argmin(abs(Tg - Tref))
    OUT[f"n_local/T^3 at T_h={Tref} [Q-derived]"] = float(nl[j] / Tg[j] ** 3)
OUT["n_c,cold/T_c^3 used by relic_yield [D]"] = cold_ratio

# ---------------------------------------------------------------- 3. thermal nucleation in the metastable window
# thin-wall O(3) bounce for a +-v vacuum bubble inside the homogeneous phi=0 charged phase (canonical units
# ell = 1/(sqrt(lphi) v), rho in lphi v^4, n in sqrt(lphi) v^3).  P_hom = mu n - rho on branch 2 [D].
b = lX / lphi
sigQ = 0.680013          # [Q] charge_phase / v_charge (Q-matter/vacuum tension, sqrt(lphi) v^3 units)
nuc = {}
for nbar in (0.3162, 0.5, 1.0, 1.5):
    rho = 0.25 + 0.75 * b ** (1 / 3) * nbar ** (4 / 3)
    mu = b ** (1 / 3) * nbar ** (1 / 3)
    eps = -(mu * nbar - rho)            # vacuum favoured by eps (lphi v^4)
    # S3/T = 16 pi sig^3/(3 eps^2 T); sig^3/eps^2 carries (sqrt(l)v^3)^3/(l v^4)^2 = v/sqrt(l)
    S3_over_T_times_T = 16 * math.pi * sigQ ** 3 / (3 * eps ** 2) * (1 / math.sqrt(lphi))   # in units of v
    S4 = 27 * math.pi ** 2 * sigQ ** 4 / (2 * eps ** 3) / lphi
    nuc[f"n={nbar * 0.1:.4f} v^3"] = {"eps (lphi v^4)": eps, "S4 [E]": S4,
                                      "S3/T at T_h=0.15v [E]": S3_over_T_times_T / 0.15,
                                      "S3/T at T_h=0.4v [E]": S3_over_T_times_T / 0.4}
OUT["thin_wall_nucleation_in_charged_phase [E]"] = nuc

# ---------------------------------------------------------------- 4. X bound state on the planar wall -----------
# -psi'' + mX^2 tanh^2(z/delta) psi = w^2 psi ; delta = sqrt(2/lphi)/v ; Poschl-Teller l(l+1) = mX^2 delta^2 = 2 gX^2/lphi
delta = math.sqrt(2 / lphi) / v
Lz, Nz = 60 * delta, 6000
z = np.linspace(-Lz / 2, Lz / 2, Nz)
dz = z[1] - z[0]
for name, depth2, m2 in (("X", (gX * v) ** 2, (gX * v) ** 2), ("D", (hD * v) ** 2, (gX * v) ** 2 + (hD * v) ** 2)):
    Vz = m2 - depth2 / np.cosh(z / delta) ** 2
    d = 2 / dz ** 2 + Vz
    e = -np.ones(Nz - 1) / dz ** 2
    w2 = eigh_tridiagonal(d, e, select="i", select_range=(0, 2))[0]
    lam_PT = 0.5 * (-1 + math.sqrt(1 + 4 * depth2 * delta ** 2))
    OUT[f"wall bound states of {name} [C]"] = {
        "PT l": lam_PT, "omega_lowest_GeV": float(math.sqrt(w2[0])), "bulk_mass_GeV": math.sqrt(m2),
        "analytic_lowest_GeV (PT)": math.sqrt(m2 - (lam_PT / delta) ** 2),
        "reflectionless (integer l)": abs(lam_PT - round(lam_PT)) < 1e-9}

# ---------------------------------------------------------------- 5. spectator-condensate relic for H_inf >> v [E]
# After inflation the light spectators (X, phi, D) carry <U> = c H^4 (c = 0.042-0.094, hidden_heating [Q]).
# Instant reheating: rho_h/rho_r = (c/3)(H/M_P)^2; quartic oscillations redshift like radiation; quanta of energy
# E ~ k_osc ~ sqrt(lX) |X| ~ 0.08 H (estimate); hidden number frozen (no number-changing; X, D, h stable).
spec = {}
for Hinf in (1e3, 1e6, 1e8, 1e10, 1e11, 1e13):
    c = 0.07
    ratio = c / 3 * (Hinf / MP) ** 2
    TRH = (90 / (math.pi ** 2 * gstar)) ** 0.25 * math.sqrt(Hinf * MP)
    E = 0.08 * Hinf
    Y = ratio * 0.75 * TRH / E          # n/s = (rho_h/E)/s, rho_r/s = 3T/4
    spec[f"H_inf={Hinf:g}"] = {"rho_h/rho_r [E]": ratio, "T_RH_instant": TRH, "Y_N [E]": Y,
                                "Y_N / 4.30e-12 (relic bound, relic_yield) [E]": Y / 4.30e-12}
OUT["spectator_condensate_relic [E]"] = spec

# ---------------------------------------------------------------- 6. reference history numbers -------------------
ref = {}
xi = 5e-4
gh = 4
ref["xi"] = xi
ref["B_h for this xi (thermal form, hidden_heating formula) [D]"] = (gh / gstar) * xi ** 4
ref["T_SM at hidden T_c (T_c=1.560 v) GeV [D]"] = 1.560 * v / xi
ref["xi_kin(T_h) at T_h = T_c: 6.6e-5 [Q]; satisfied"] = xi > 6.6e-5
ref["omega_sym = 0.272 (xi/1e-3)^3 [Q relic_yield scaling near the bound]"] = 0.272 * (xi / 1e-3) ** 3
for r in (1e-4, 0.05):
    om_asym = mstar * r * Y_B * 2.7437e8
    ref[f"r={r:g}: omega_asym >= m* r Y_B ... [D]"] = om_asym
    ref[f"r={r:g}: maser product LHS m_X r + 79.28 omega_sym GeV (< 63.93) [D]"] = mX * r + 79.28 * 0.272 * (xi / 1e-3) ** 3
    Th = exit_Th(r, xi)
    ref[f"r={r:g}: T_h,exit/v [C]"] = Th
    ref[f"r={r:g}: T_SM,exit GeV [C]"] = (Th * v / xi) if Th else 1.560 * v / xi
    Tex = (Th * v / xi) if Th else 1.560 * v / xi
    Hx = H_of_T(Tex, gstar)
    ref[f"r={r:g}: H/m_X at exit [C]"] = Hx / mX
    ref[f"r={r:g}: H^-1/delta at exit [C]"] = 1 / (Hx * delta)
# hidden energy carried by the charged condensate at T_c relative to the thermal bath (branch-2 homogeneous) [D]
for r in (0.05, 1.0):
    for x in (5e-4, 7.58e-4):
        Tc = 1.560
        nbar = nQ_over_Th3(r, x) * Tc ** 3 / 0.1    # in sqrt(lphi) v^3 units
        rho_c = lphi * (0.25 + 0.75 * b ** (1 / 3) * nbar ** (4 / 3))   # v^4
        rho_th = math.pi ** 2 / 30 * gh * Tc ** 4
        ref[f"rho_condensate/rho_thermal at T_c, r={r:g}, xi={x:g} [D]"] = rho_c / rho_th
OUT["reference_history"] = ref

# cold complementary branch
cold = {}
for Hinf in (120.0, 138.0, 184.0):
    TRH = (90 / (math.pi ** 2 * gstar)) ** 0.25 * math.sqrt(Hinf * MP)
    cold[f"H_inf={Hinf}: instant T_RH GeV [D]"] = TRH
    cold[f"H_inf={Hinf}: N_roll ~ (3H^2/(l v^2)) ln(2 pi v/H) [Q formula]"] = 3 * Hinf ** 2 / (lphi * v ** 2) * math.log(2 * math.pi * v / Hinf)
cold["T_n(r=0.05) cold TeV [D]"] = T_n_cold(0.05) / 1e3
cold["T_s1(r=0.05) = T_n (n_s1/n_c)^(1/3) TeV [D]"] = T_n_cold(0.05) / 1e3 * (0.0427667 / 0.0316228) ** (1 / 3)
OUT["cold_branch"] = cold

json.dump(OUT, open(os.path.join(HERE, "critic_results.json"), "w"), indent=1, default=float)
print(json.dumps(OUT, indent=1, default=float)[:12000])
