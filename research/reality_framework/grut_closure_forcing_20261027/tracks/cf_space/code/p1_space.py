"""P1: spatial compactness (S^3 or T^3), its global constraints, and the combined P1 + P2 curvature/topology
consequences.  Run: python3 -I p1_space.py   (writes ../results/p1_space.json)
"""
import itertools
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cfs_common import (G_GRUT, G_SM, H0_PLANCK, HEAVY, M_P, H0_GeV, V_of_H, ln_aH_end_over_a0H0, S_dS, H_lambda)

out = {}
GPC_INV_GEV = 3.0856775814913673e25 / 1.97326980e-16   # 1 Gpc in GeV^-1
c_over_H0_Gpc = 299792.458 / H0_PLANCK / 1000.0

# ============================================================== 1. curvature data (retrospective)
# [Q] Planck 2018 VI (1807.06209) eqs. 46b, 47a, 47b; DESI DR2 (2503.14738) Table 5
data = {
    "Planck18_TTTEEE_lowE_eq46b": (-0.044, 0.0165),          # -0.044 +0.018 -0.015 (symmetrised)
    "Planck18_+lensing_eq47a": (-0.0106, 0.0065),
    "Planck18_+lensing+BAO_eq47b": (0.0007, 0.0019),
    "DESI_DR2+CMB_LCDM+OmK_Table5": (0.0023, 0.0011),
    "DESI_DR2+CMB+PantheonPlus_w0wa+OmK": (0.0011, 0.0013),
    "DESI_DR2+CMB+Union3_w0wa+OmK": (0.0006, 0.0013),
    "DESI_DR2+CMB+DESY5_w0wa+OmK": (0.0008, 0.0013),
}
curv = {}
for k, (mu, sig) in data.items():
    curv[k] = {
        "mean": mu, "sigma": sig,
        "P_OmK_lt_0_(S3_side)": float(norm.cdf(-mu / sig)),
        "z_of_zero": mu / sig,
        "S3_OmK_95lower": mu - 1.645 * sig,   # one-sided 95% lower limit on Omega_K
        "S3_Rc_min_Gpc_if_closed": (c_over_H0_Gpc / math.sqrt(abs(min(mu - 1.645 * sig, -1e-12))))
        if mu - 1.645 * sig < 0 else None,
    }
out["curvature_data"] = curv
out["c_over_H0_Gpc"] = c_over_H0_Gpc
# cosmic-variance floor [Q via search]: Waterhouse-Zibin 0804.1771 ~1.5e-5 (1 sigma); VTS 0901.3354 ~1e-4 practical
out["curvature_floor"] = {"cosmic_variance_1sigma": 1.5e-5, "practical_floor": 1e-4}

# ============================================================== 2. distance to last scattering, topology bounds
h = H0_PLANCK / 100.0
Om, Ogam = 0.3153, 2.472e-5 / h ** 2
Orad = Ogam * (1 + 0.2271 * 3.046)
OL = 1 - Om - Orad
zstar = 1089.9
chi_rec_H0 = quad(lambda z: 1 / math.sqrt(Om * (1 + z) ** 3 + Orad * (1 + z) ** 4 + OL), 0, zstar, limit=500)[0]
topo = {
    "chi_rec_in_c_over_H0": chi_rec_H0,
    "chi_rec_Gpc": chi_rec_H0 * c_over_H0_Gpc,
    "d_LSS_Gpc": 2 * chi_rec_H0 * c_over_H0_Gpc,
    # COMPACT Part I (2211.02603): shortest loop through us > 0.985 d_LSS [Q]
    "L_through_us_min_in_c_over_H0": 0.985 * 2 * chi_rec_H0,
    # Planck 2015 XVIII (1502.01593): cubic torus R_i > 0.97 chi_rec (99% CL) -> L = 2 R_i [Q]
    "L_cubic_T3_min_in_c_over_H0": 2 * 0.97 * chi_rec_H0,
}
out["topology"] = topo

# ============================================================== 3. P1 + P2 combined (compact space + capacity cap)
# With a cap N <= N_min + f on inflation (f from the capacity postulate), and an initial compact size
#   S^3: Omega_k,i = 1/(a_i R_c H_i)^2 = 1/ell^2,   T^3: side L_i = ell / H_i,
# one has exactly |Omega_k,0| = |Omega_k,i| exp(-2 (N - N_min)) and L_0 H_0 = ell exp(N - N_min).


def f_TCC(H, g=G_GRUT):
    """TCC headroom: max (N - N_min) = ln(M_P/H) - N_min(H), instant reheating."""
    return math.log(M_P / H) - ln_aH_end_over_a0H0(H, g_rh=g)


def H_at_f(ftarget):
    return 10 ** brentq(lambda lx: f_TCC(10 ** lx) - ftarget, -12, 18)


comb = {"TCC_headroom_note": "f_TCC(H) = (3/4) ln(V_TCC/V) to within dof factors; zero at the TCC limit"}
Hb = H_at_f(0.0)
comb["V_quarter_TCC_limit"] = V_of_H(Hb) ** 0.25
# S^3, ell = 1 (initial S^3 radius = Hubble radius): curvature floor today
s3 = {}
for x in (1e-2, 3e-3, 1e-3, 1e-4, 1.5e-5):
    fH = 0.5 * math.log(1.0 / x)
    Hx = H_at_f(fH)
    s3[f"|OmK0|={x:g}"] = {"f_needed": fH, "V_quarter_GeV": V_of_H(Hx) ** 0.25,
                           "V_quarter_over_limit": V_of_H(Hx) ** 0.25 / comb["V_quarter_TCC_limit"]}
comb["S3_ell1_floor_reached_at"] = s3
comb["S3_note"] = ("under TCC + S^3 + Omega_k,i = 1, |Omega_K,0| >= x only for V^{1/4} above the listed value "
                   "(and below the TCC limit); below it the floor is unobservable")
# DESI DR2 + CMB tension of a forced closed value Omega_K0 = -x
comb["S3_DESI_sigma_for_x"] = {f"{x:g}": (0.0023 + x) / 0.0011 for x in (1e-4, 1e-3, 3e-3, 1e-2)}
comb["S3_Planck+lens+BAO_sigma_for_x"] = {f"{x:g}": (0.0007 + x) / 0.0019 for x in (1e-4, 1e-3, 3e-3, 1e-2)}
# T^3 with ell = 1: L0 H0 <= e^f; topology null (through-us loop) requires L0 H0 > L_min
Lmin = topo["L_through_us_min_in_c_over_H0"]
fT = math.log(Lmin)
HT = H_at_f(fT)
comb["T3_ell1"] = {
    "L0H0_max_at_TCC_limit": 1.0,
    "f_needed_for_L_above_COMPACT_bound": fT,
    "V_quarter_max_GeV_for_null_topology": V_of_H(HT) ** 0.25,
    "ratio_to_TCC_limit": V_of_H(HT) ** 0.25 / comb["V_quarter_TCC_limit"],
}
# same for Banks-Fischler: cap Delta = 26.11 independent of scale -> floors unobservable
Delta = 26.114
comb["BF"] = {"Delta": Delta, "S3_floor_factor": math.exp(-2 * Delta), "T3_ceiling_L0H0_over_ell": math.exp(Delta),
              "note": "BF bounds N - N_min from above only; the curvature floor e^{-2 Delta} ~ 2e-23 |Omega_k,i| "
                      "is far below the cosmic-variance floor"}
out["P1_plus_P2"] = comb

# ============================================================== 4. Gauss law on compact slices (hidden U(1): 1, 3, 6)
rng = np.random.default_rng(13)
n = 16
E = rng.normal(size=(3, n, n, n))
div_periodic = sum(np.roll(E[i], -1, axis=i) - E[i] for i in range(3))   # forward-difference divergence on T^3
# open box with the same field: the sum of the divergence equals the boundary flux, not zero
div_open = sum(float(np.diff(E[i], axis=i).sum()) for i in range(3))   # = net flux through the box faces
gauss = {
    "T3_lattice_total_charge_sum_div": float(div_periodic.sum()),
    "T3_lattice_max_abs_div": float(np.abs(div_periodic).max()),
    "open_box_nonzero_example": float(div_open),
    "homogeneous_note": "for spatially constant E, div E = 0 identically, so rho = 0 pointwise on R^3 or T^3",
}
# Higgs phase: the condensate carries charge in multiples of 6: allowed (N_X, N_Y) residues
allowed = [(a, b) for a, b in itertools.product(range(6), range(6)) if (a + 3 * b) % 6 == 0]
gauss["Z6_allowed_residues_(NX mod 6, NY mod 6)"] = allowed
gauss["fraction_of_residue_classes_allowed"] = len(allowed) / 36.0
# number of flat Z6 bundles (holonomy sectors) on T^3 = |Hom(Z^3, Z6)| = 6^3
gauss["Z6_holonomy_sectors_T3"] = 6 ** 3
# flux quantisation on a 2-torus with a charge-1 field present: integral F = 2 pi n / e; a Z6 string carries 2 pi/(6e)
gauss["net_Z6_string_crossings_per_2cycle_mod"] = 6
out["gauss_law"] = gauss

# ============================================================== 5. countermodel: finite capacity without Lambda > 0
# k = +1 dust universe, Lambda = 0, Omega_m0 = 1 + delta: turnaround radius r_max = (c/H0)(1+delta)/delta^{3/2};
# its largest causal-diamond 2-sphere is the S^3 equator at turnaround: CEB capacity A/4G = 8 pi^2 M_P^2 r_max^2.
cm = {}
for delta in (0.01, 0.1, 1.0):
    rmax = c_over_H0_Gpc * (1 + delta) / delta ** 1.5 * GPC_INV_GEV
    cm[f"delta={delta}"] = {"r_max_Gpc": rmax / GPC_INV_GEV, "CEB_capacity": 8 * math.pi ** 2 * M_P ** 2 * rmax ** 2}
cm["S_dS_observed_Lambda"] = S_dS(H_lambda())
cm["note"] = "finite for every delta > 0: finite capacity does not force Lambda > 0 unless a crunch is also excluded"
out["countermodel_closed_recollapse"] = cm

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "p1_space.json")
with open(path, "w") as f:
    json.dump(out, f, indent=1, default=str)


def show(d, pre=""):
    for k, v in d.items():
        if isinstance(v, dict):
            show(v, pre + k + ".")
        else:
            print(f"{pre}{k} = {v:.6g}" if isinstance(v, float) else f"{pre}{k} = {v}")


show(out)
