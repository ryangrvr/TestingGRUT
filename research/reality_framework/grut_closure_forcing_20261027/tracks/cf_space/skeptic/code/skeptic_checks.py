"""Skeptic re-derivation of track cf_space (closure in space and capacity).

Independent of ../../code/*.py: constants are rebuilt from SI values, the e-fold bookkeeping is redone with a
piecewise (ln a, ln H) timeline, and the composite (multi-interval) form of TCC is checked directly.

Run:  python3 -I skeptic_checks.py     (writes ../results/skeptic_checks.json)
Labels used in comments: [C] computed here, [Q] quoted (source in SKEPTIC.md), [R] recalled.
"""
import itertools
import json
import math
import os

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import norm

out = {}
LN10 = math.log(10.0)

# ------------------------------------------------------------------ constants from SI (CODATA 2018) [R]
G_SI = 6.67430e-11
HBAR = 1.054571817e-34
C = 2.99792458e8
EV = 1.602176634e-19
KB = 1.380649e-23
GEV_J = 1e9 * EV
M_P = math.sqrt(HBAR * C / (8 * math.pi * G_SI)) * C ** 2 / GEV_J        # reduced Planck mass, GeV
M_PL = math.sqrt(HBAR * C / G_SI) * C ** 2 / GEV_J                        # non-reduced
MPC_M = 3.0856775814913673e22
T0 = 2.7255 * KB / GEV_J                                                    # GeV [R, Fixsen 2009]
NEFF = 3.044
G_S0 = 2 + 7 / 8 * 2 * NEFF * 4 / 11                                       # 3.938 (N_eff = 3.044)
G106 = 106.75


def H0_gev(h=0.6736):
    return h * 100e3 / MPC_M * HBAR / GEV_J        # s^-1 -> GeV


out["constants"] = {"M_P_reduced_GeV": M_P, "M_P_nonreduced_GeV": M_PL, "T0_GeV": T0, "g_s0_Neff3.044": G_S0,
                    "H0_GeV_h0.6736": H0_gev()}


# ------------------------------------------------------------------ piecewise timeline in (ln a, ln H)
def rho_rad(T, g=G106):
    return math.pi ** 2 / 30 * g * T ** 4


def T_of_rho(rho, g=G106):
    return (30 * rho / (math.pi ** 2 * g)) ** 0.25


def H_of_rho(rho, MP=M_P):
    return math.sqrt(rho / 3) / MP


def ln_a0_over_a_rh(T_rh, g_s_rh=G106):
    """entropy conservation from reheating to today"""
    return math.log(T_rh / T0) + math.log(g_s_rh / G_S0) / 3


def N_min(H_end, T_rh=None, w=None, g=G106, h=0.6736, MP=M_P):
    """ln(a_end H_end / a0 H0): e-folds needed for today's Hubble scale to have been inside the Hubble radius at
    the start of a constant-H stage ending at H_end. Optional constant-w phase between H_end and T_rh."""
    rho_end = 3 * MP ** 2 * H_end ** 2
    if T_rh is None:
        T_rh = T_of_rho(rho_end, g)
    rho_rh = rho_rad(T_rh, g)
    ln_aend_over_arh = 0.0 if w is None else -math.log(rho_end / rho_rh) / (3 * (1 + w))
    return math.log(H_end / H0_gev(h)) + ln_aend_over_arh - ln_a0_over_a_rh(T_rh, g)


# ------------------------------------------------------------------ A. heavy-window kinematics, cross-checked
H_i = 1.3335e11
A_s = math.exp(3.044) * 1e-10          # Planck 2018 ln(1e10 A_s) = 3.044 [R]
V_heavy = 3 * M_P ** 2 * H_i ** 2
A = {}
A["V_quarter"] = V_heavy ** 0.25
A["r"] = 2 * H_i ** 2 / (math.pi ** 2 * M_P ** 2 * A_s)
A["T_inst_106.75"] = T_of_rho(V_heavy)
A["N_min_exact"] = N_min(H_i)
# Liddle-Leach (astro-ph/0305263) normalisation: N = 67 - ln(k/a0H0) + 1/4 ln(V^2/(M^4 rho_end)) + 1/12 ln(rho_reh/rho_end) - 1/12 ln g
A["N_min_LiddleLeach67"] = 67 + 0.25 * math.log(V_heavy / M_P ** 4) - math.log(G106) / 12
# k = 0.05/Mpc pivot
k_over_a0H0 = 0.05 / (0.6736 * 100 / 299792.458)
A["N_star_k0.05"] = A["N_min_exact"] - math.log(k_over_a0H0)
out["A_heavy_window"] = A

# ------------------------------------------------------------------ B. TCC single stage, conventions
B = {}
for tag, MP in (("reduced", M_P), ("nonreduced", M_PL)):
    lnH = brentq(lambda x: math.log(MP / math.exp(x)) - N_min(math.exp(x), MP=M_P), -20, 40)
    Hb = math.exp(lnH)
    B[f"standard_{tag}"] = {"H_max": Hb, "V_quarter_max": (3 * M_P ** 2 * Hb ** 2) ** 0.25,
                            "r_max": 2 * Hb ** 2 / (math.pi ** 2 * M_P ** 2 * A_s),
                            "N_TCC_heavy": math.log(MP / H_i), "deficit_heavy": N_min(H_i) - math.log(MP / H_i)}
# most permissive history: w = -1/3 (aH constant) from H_end down to reheating at T_BBN. N_min then equals the
# N_min of a stage that ends exactly at T_BBN, independent of H_end. [C]
for TB in (4.7e-3, 5e-3):
    H_bbn = H_of_rho(rho_rad(TB, 10.75))
    nmin_relax = math.log(H_bbn / H0_gev()) - ln_a0_over_a_rh(TB, 10.75)
    B[f"relaxed_TBBN_{TB*1e3:.1f}MeV"] = {
        "N_min_floor_any_history": nmin_relax,
        "H_end_max_reduced": M_P * math.exp(-nmin_relax),
        "H_end_max_nonreduced": M_PL * math.exp(-nmin_relax),
        "heavy_over_bound_reduced": H_i / (M_P * math.exp(-nmin_relax)),
        "heavy_over_bound_nonreduced": H_i / (M_PL * math.exp(-nmin_relax)),
        "margin_efolds_reduced": nmin_relax - math.log(M_P / H_i),
        "margin_efolds_nonreduced": nmin_relax - math.log(M_PL / H_i),
    }
# pre-inflationary radiation era from Planck density: modes Planckian at t_P are already stretched by (M_P/H)^{1/2}
B["preinflation_radiation_from_Planck_density_N_max_heavy"] = 0.5 * math.log(M_P / H_i)
out["B_TCC_single"] = B


# ------------------------------------------------------------------ C. composite TCC for multi-stage histories
# Literal TCC on every interval [t_i, t_f]: a_f/a_i < M_P/H_f.  Since a is monotone, the strongest instance takes
# t_i = start of the history, so TCC <=> max_t a(t)H(t) < M_P a_start.  For a history that starts with a stage at
# H_1 whose initial Hubble patch contains our present Hubble volume (a_start H_1 <= a0 H0), and ends with any later
# stage(s) and reheating:  max aH >= a_last,end H_last = e^{N_min(last)} a0 H0, so TCC requires
#     e^{N_min(last)} < M_P / H_1 ,   and N_min(last) >= N_min floor (w = -1/3 to BBN) = 19.38.
def timeline_check(H1, N1, T_on2, H2, N2, T_rh2, w_mid=None):
    """Build a two-stage history explicitly and test (i) horizon, (ii) composite TCC by brute force over a grid."""
    pts = []                         # (ln a, ln H) at breakpoints; aH is monotone between breakpoints
    la = 0.0
    pts.append((la, math.log(H1)))
    la += N1
    pts.append((la, math.log(H1)))                                   # end stage 1, instant reheating
    T_rh1 = T_of_rho(3 * M_P ** 2 * H1 ** 2)
    la += math.log(T_rh1 / T_on2)                                    # radiation (g const) to onset of stage 2
    pts.append((la, math.log(H_of_rho(rho_rad(T_on2)))))
    pts.append((la, math.log(H2)))
    la += N2
    pts.append((la, math.log(H2)))
    rho_end2 = 3 * M_P ** 2 * H2 ** 2
    rho_rh2 = rho_rad(T_rh2)
    w = 0.0 if w_mid is None else w_mid
    la += math.log(rho_end2 / rho_rh2) / (3 * (1 + w))
    pts.append((la, math.log(H_of_rho(rho_rh2))))
    la_now = la + ln_a0_over_a_rh(T_rh2)
    pts = np.array(pts)
    lnaH = pts[:, 0] + pts[:, 1]
    horizon_margin = (la_now + math.log(H0_gev())) - (pts[0, 0] + pts[0, 1])   # >= 0 means our patch inside
    tcc_violation = np.max(lnaH) - (pts[0, 0] + math.log(M_P))                # < 0 means TCC satisfied
    return horizon_margin, tcc_violation


Cc = {}
# scan: stage 1 = heavy window (N1 <= 16.72 as single-stage TCC allows); stage 2 any H2, any N2, any reheating
rng = np.random.default_rng(2026)
n_ok_both = 0
best = -np.inf
for _ in range(20000):
    N1 = rng.uniform(1, math.log(M_P / H_i))
    H2 = 10 ** rng.uniform(-23, 10.9)
    T_on2 = T_of_rho(3 * M_P ** 2 * H2 ** 2)
    Tmax2 = T_on2
    T_rh2 = 10 ** rng.uniform(math.log10(5e-3), math.log10(Tmax2)) if Tmax2 > 5e-3 else None
    if T_rh2 is None:
        continue
    N2 = rng.uniform(0, 80)
    hm, tv = timeline_check(H_i, N1, T_on2, H2, N2, T_rh2)
    if hm >= 0:
        best = max(best, -tv)
        if tv < 0:
            n_ok_both += 1
Cc["random_two_stage_histories_with_horizon_solved_by_stage1_patch_and_TCC_ok"] = n_ok_both
Cc["best_TCC_margin_among_horizon_solving_histories_(negative=violated)"] = best
Cc["analytic_requirement"] = "e^{N_min(last stage)} < M_P/H_1; N_min floor 19.38 > ln(M_P/H_i) = 16.72"
# alternative: stage 2 alone solves the horizon problem (stage-1 patch smaller than our Hubble volume)
# relic yield dilution D = e^{3 N2} (s_rh2/s_on2) >= e^{3 N_min(stage 2)}; least dilution at the lowest stage-2 scale
dil = []
for lT in np.linspace(math.log10(5e-3), 9, 200):
    T2 = 10 ** lT
    H2 = H_of_rho(rho_rad(T2, 10.75 if T2 < 0.2 else G106))
    n2 = N_min(H2, g=10.75 if T2 < 0.2 else G106)
    dil.append((3 * n2 / LN10, T2, n2))
dmin = min(dil)
Cc["stage2_only_min_log10_dilution"] = dmin[0]
Cc["stage2_only_at_T_rh2_GeV"] = dmin[1]
Cc["stage2_only_N2"] = dmin[2]
Cc["track_claimed_log10_dilution"] = 3 * (N_min(H_i) - math.log(M_P / H_i)) / LN10
# Y ~ H^{3/2}: H boost to recover the yield
Cc["stage2_only_required_H_boost_log10"] = dmin[0] / 1.5
Cc["stage2_only_required_H_GeV"] = H_i * 10 ** (dmin[0] / 1.5)
# how many stage-1 Hubble patches (at the end of stage 1, N1 = 16.72) make up our Hubble volume if stage 2 alone
# solves the horizon problem with T_rh2 = 5 MeV: ratio of comoving scales (a0H0)^-1 / (a_i1 H1)^-1
T2 = 5e-3
H2 = H_of_rho(rho_rad(T2, 10.75))
# comoving scale of our Hubble volume relative to the stage-1 initial Hubble patch:
#   with stage 2 minimal (a_on2 H2 = a0 H0) and instant reheating in both stages:
#   a_end1 = a_on2 (T_on2/T_rh1)  =>  (a_i1 H1)/(a0 H0) = e^{-N1} (H1/H2) (T_on2/T_rh1)
T_rh1 = T_of_rho(3 * M_P ** 2 * H_i ** 2)
ln_ratio = -math.log(M_P / H_i) + math.log(H_i / H2) + math.log(T2 / T_rh1)
# consistency: this history obeys the composite TCC (max aH over a_start is below M_P)
Cc["stage2_only_composite_TCC_ln_margin_end_stage1_(0_by_construction_N1=N_TCC)"] = 0.0
Cc["stage2_only_composite_TCC_ln_margin_end_stage2"] = math.log(M_P) - (
    N_min(H2, g=10.75) + math.log(H_i) - ln_ratio)                          # ln M_P - ln(a_end2 H2 / a_i1)
Cc["stage2_only_log10_number_of_stage1_patches_in_our_Hubble_volume"] = 3 * ln_ratio / LN10
out["C_composite_TCC"] = Cc

# ------------------------------------------------------------------ D. Banks-Fischler vs the covariant (diamond) count
h = 0.6736
Om, OL = 0.3153, 0.6847
Or = 2.47e-5 / h ** 2 * (1 + 0.2271 * NEFF)
E = lambda a: math.sqrt(Om / a ** 3 + Or / a ** 4 + OL)
chi_EH = quad(lambda a: 1 / (a ** 2 * E(a)), 1, np.inf)[0]                        # future event horizon, c/(a0H0)
chi_P = quad(lambda a: 1 / (a ** 2 * E(a)), 1e-30, 1, limit=400, points=[1e-10, 1e-6, 1e-3])[0]  # particle horizon
s0 = 2 * math.pi ** 2 / 45 * G_S0 * T0 ** 3
H0 = H0_gev(h)
S_H0 = s0 * 4 * math.pi / 3 / H0 ** 3
HL = H0 * math.sqrt(OL)
S_L = 8 * math.pi ** 2 * M_P ** 2 / HL ** 2
D = {"chi_event_horizon": chi_EH, "chi_particle_horizon": chi_P,
     "S_rad_today_Hubble_sphere": S_H0, "S_Lambda": S_L,
     "S_rad_in_eternal_observer_diamond_max": S_H0 * (chi_EH + chi_P) ** 3,
     "ratio_diamond_to_S_Lambda": S_H0 * (chi_EH + chi_P) ** 3 / S_L,
     "BF_Delta_(1/3)ln(S_L/S_H0)": math.log(S_L / S_H0) / 3,
     "note": "radiation entropy inside the causal diamond of an eternal comoving observer is bounded by the comoving "
             "ball of radius chi_EH + chi_P for ANY N >= N_min + ln(chi_EH+chi_P): the covariant N-bound does not cap N"}
out["D_BF_vs_diamond"] = D

# ------------------------------------------------------------------ E. closed-dust countermodel, numerical diamond
E_ = {}
c_over_H0_Gpc = 299792.458 / (h * 100) / 1000
GPC_GEVINV = 3.0856775814913673e25 / (HBAR * C / GEV_J)
for delta in (1.0, 0.1, 0.01):
    Omega_m = 1 + delta
    Rc = c_over_H0_Gpc / math.sqrt(delta)              # curvature radius today, Gpc (a0 = 1)
    a_max = Omega_m / (Omega_m - 1)                    # turnaround scale factor
    # a(eta) = (a_max/2)(1 - cos eta); diamond of observer from eta=0 to 2pi: chi_D = min(eta, 2pi - eta)
    f = lambda eta: -((a_max / 2) * (1 - math.cos(eta))) ** 2 * (math.sin(min(eta, 2 * math.pi - eta, math.pi / 2))) ** 2
    res = minimize_scalar(f, bounds=(0.01, 2 * math.pi - 0.01), method="bounded")
    rmax = math.sqrt(-res.fun) * Rc * GPC_GEVINV
    E_[f"delta={delta}"] = {"eta_at_max": res.x, "area_radius_Gpc": rmax / GPC_GEVINV,
                            "capacity_A_over_4G": 8 * math.pi ** 2 * M_P ** 2 * rmax ** 2}
out["E_countermodel"] = E_

# ------------------------------------------------------------------ F. curvature data, extended (retrospective)
data = {
    "Planck18_TTTEEE_lowE_eq46b_upper_err": (-0.044, 0.018),
    "DESI_DR2_Table5_CMB_only": (-0.0107, 0.0064),
    "DESI_DR2+CMB_LCDM+OmK_Table5": (0.0023, 0.0011),
    "Chudaykin+_2511.20757_CMB+BAO": (0.0024, 0.0012),
    "Chudaykin+_2511.20757_CMBPR4+FS+BAO": (0.0022, 0.0011),
    "Chudaykin+_2511.20757_CMB+FS+BAO": (0.0027, 0.0012),
    "ExpQuintessence_2504.04226_OmK_free": (0.003, 0.001),
    "DESI_DR2+CMB+DESY5_w0wa+OmK": (0.0008, 0.0013),
}
F = {}
for k, (mu, s) in data.items():
    F[k] = {"P_OmK_lt_0": float(norm.cdf(-mu / s)), "z_zero": mu / s}
out["F_curvature"] = F

# ------------------------------------------------------------------ G. TCC + S3 band, recomputed
def f_tcc(Hx):
    return math.log(M_P / Hx) - N_min(Hx)


Hlim = math.exp(brentq(lambda x: f_tcc(math.exp(x)), -20, 40))
Vlim = (3 * M_P ** 2 * Hlim ** 2) ** 0.25
Gd = {"V_quarter_TCC_limit": Vlim}
for x in (1e-3, 1e-4):
    Hx = math.exp(brentq(lambda y: f_tcc(math.exp(y)) - 0.5 * math.log(1 / x), -40, 40))
    Gd[f"V_quarter_for_OmK_floor_{x:g}"] = (3 * M_P ** 2 * Hx ** 2) ** 0.25
Lmin = 0.985 * 2 * 3.1158
Ht = math.exp(brentq(lambda y: f_tcc(math.exp(y)) - math.log(Lmin), -40, 40))
Gd["V_quarter_max_null_topology_T3_ell1"] = (3 * M_P ** 2 * Ht ** 2) ** 0.25
out["G_TCC_plus_compact"] = Gd

# ------------------------------------------------------------------ H. family floors vs relaxed bounds
c6, c6_ceiling, eps = 1.2504e6, 1.486e13, 0.146
Hf = {"H_min_c6": H_i * math.sqrt(c6 / c6_ceiling), "H_min_yield_eps1": H_i * eps ** (2 / 3)}
for tag in ("reduced", "nonreduced"):
    Hb = B["relaxed_TBBN_5.0MeV"][f"H_end_max_{tag}"]
    Hf[f"relaxed_bound_{tag}"] = Hb
    Hf[f"yield_floor_over_relaxed_bound_{tag}"] = Hf["H_min_yield_eps1"] / Hb
out["H_family_floors"] = Hf

# ------------------------------------------------------------------ I. Gauss law, FFT method, Z6 counting
n = 24
rng = np.random.default_rng(7)
Ef = rng.normal(size=(3, n, n, n))
k = 2 * np.pi * np.fft.fftfreq(n)
div_hat = sum(1j * np.sin(k).reshape([-1 if d == i else 1 for d in range(3)]) * np.fft.fftn(Ef[i])
              for i in range(3))                 # central-difference divergence in Fourier space
Ig = {"FFT_zero_mode_of_divE": abs(div_hat[0, 0, 0]),
      "allowed_residues": sum(1 for a, b in itertools.product(range(6), repeat=2) if (a + 3 * b) % 6 == 0),
      "Z6_flat_bundles_T3": len(list(itertools.product(range(6), repeat=3)))}
out["I_gauss"] = Ig

# ------------------------------------------------------------------ J. late-time numbers
J = {"S_Lambda": S_L, "ln_MP_over_HL": math.log(M_P / HL),
     "t_TCC_yr": math.log(M_P / HL) / HL * HBAR / GEV_J / 3.15576e7,
     "rho_Lambda_GeV4": 3 * M_P ** 2 * HL ** 2,
     "premise_rho_vac_min_GeV4": -2.72e-47,
     "LQT_2506.24011_mean_OmL_-0.69_GeV4": -0.69 * 3 * M_P ** 2 * H0_gev(0.68) ** 2,
     "LQT_2506.24011_bestfit_OmL_-1.61_GeV4": -1.61 * 3 * M_P ** 2 * H0_gev(0.68) ** 2}
out["J_late"] = J

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "skeptic_checks.json")
with open(path, "w") as fh:
    json.dump(out, fh, indent=1, default=float)


def show(d, pre=""):
    for kk, v in d.items():
        if isinstance(v, dict):
            show(v, pre + kk + ".")
        else:
            print(f"{pre}{kk} = {v:.6g}" if isinstance(v, float) else f"{pre}{kk} = {v}")


show(out)
