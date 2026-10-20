#!/usr/bin/env python3
"""review.py -- Claude's independent review calculations for GitHub issue #13 (Candidate 1 at 2ee7f82).

Written without running or importing Codex's audit.py. Inputs: the Candidate 1 action (benchmark couplings), Pesce et al. 2020
Table 1 as transcribed by Claude (data/pesce2020_table1_claude_transcription.json), Mossa et al. 2020 omega_b.

  R1  classical charge floor: m* from direct 2D minimization of U/|X|^2, closed form, and a pointwise random test of rho >= m*|n_Q|
  R2  flat-FRW floor omega_h,min and Hfloor (nominal and conservative corners)
  R3  angular-distance envelope versus the six masers (primary and stress corners), with convention sensitivity
  R4  vacuum-energy scope: the benchmark one-loop vacuum energy in omega units; the most negative total Lambda for which each
      mismatch survives
  R5  one-loop Coleman-Weinberg estimate of the homogeneous floor m*_eff (MSbar, mu = mX), cross-checked against the archived
      one-loop vacuum shift phi/v = 1.00012193740
Writes ../results/review_results.json.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq, minimize

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results"
OUT.mkdir(exist_ok=True)
t0 = time.time()

# ---------------- benchmark (GRUT_CANDIDATE_1.md section 2 and line 211; empirical protocol parameters)
g, lphi, lX = 0.1, 0.01, 0.001
mD0_over_v, hD, lD = 0.1, 0.01, 0.001
v = 1000.0            # GeV
mX = g * v            # 100 GeV
r = 1.0               # |YQ0|/YB0
c_kms = 299792.458
res = {"benchmark": {"gX": g, "lambda_phi": lphi, "lambda_X": lX, "mD0_over_v": mD0_over_v, "hD": hD, "lambda_D": lD,
                     "v_GeV": v, "mX_GeV": mX, "r": r}}


# ---------------- R1 classical floor
def U(phi, q1, q2, D):
    s2 = (q1 ** 2 + q2 ** 2) / 2.0
    return (lphi / 4) * (phi ** 2 - 1.0) ** 2 + g ** 2 * phi ** 2 * s2 + lX * s2 ** 2 + 0.5 * (mD0_over_v ** 2 + hD ** 2 * phi ** 2) * D ** 2 + (lD / 4) * D ** 4


mstar_closed = min(g, (lphi * lX) ** 0.25)          # units of v
# direct minimization of U(phi, s, D=0)/s^2 over phi and s > 0 (units v = 1)
best = (1e9, None)
for p0 in np.linspace(0, 1.2, 13):
    for ls0 in np.linspace(-6, 2, 17):
        f = lambda x: U(x[0], math.sqrt(2) * math.exp(x[1]), 0.0, 0.0) / math.exp(2 * x[1])
        o = minimize(f, [p0, ls0], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-16, "maxiter": 20000})
        if o.fun < best[0]:
            best = (o.fun, o.x)
mstar_num = math.sqrt(best[0])
# the infimum on the phi=v branch is approached only as s -> 0, value g^2; record both branches
branch_phi0 = math.sqrt(lphi * lX)  # min over s of U(0,s)/s^2 = 2*sqrt(lphi/4 * lX) = sqrt(lphi*lX) -> m* = (lphi lX)^(1/4)
rng = np.random.default_rng(20261020)
viol, worst = 0, 1e9
N = 200000
for _ in range(4):
    phi = rng.normal(0, 1.5, N); q1 = rng.normal(0, 2, N); q2 = rng.normal(0, 2, N); D = rng.normal(0, 1, N)
    dphi, dq1, dq2, dD = (rng.normal(0, 1, N) for _ in range(4))
    grad2 = rng.exponential(1.0, N)  # any non-negative gradient energy
    nQ = q2 * dq1 - q1 * dq2
    rho = 0.5 * (dphi ** 2 + dq1 ** 2 + dq2 ** 2 + dD ** 2) + grad2 + U(phi, q1, q2, D)
    ratio = rho / np.maximum(np.abs(nQ), 1e-300)
    viol += int(np.sum(rho < mstar_closed * np.abs(nQ) - 1e-12))
    worst = min(worst, float(np.min(ratio[np.abs(nQ) > 1e-6])))
# saturation: homogeneous rotating condensate at the optimum (phi = 0, s^2 = sqrt(lphi/(4 lX)), omega = m*)
s_opt2 = math.sqrt(lphi / (4 * lX))
s_opt = math.sqrt(s_opt2)
omega = mstar_closed
n_sat = omega * 2 * s_opt2                     # n = omega * f^2 with f^2 = q1^2+q2^2 = 2 s^2
rho_sat = omega ** 2 * s_opt2 + U(0.0, math.sqrt(2) * s_opt, 0, 0)
res["R1_classical_floor"] = {
    "mstar_over_v_closed_form": mstar_closed, "mstar_over_v_direct_minimization": mstar_num,
    "mu_mstar_over_mX": mstar_closed / g, "mstar_GeV": mstar_closed * v,
    "branch_values_over_v": {"phi=v branch (s->0 limit)": g, "phi=0 branch": branch_phi0 ** 0.5},
    "pointwise_random_test": {"samples": 4 * N, "violations": viol, "min_rho_over_absn_observed_over_v": worst},
    "saturation_check_rho_minus_mstar_n": rho_sat - mstar_closed * n_sat,
    "D_sector_nonnegative": bool(mD0_over_v ** 2 >= 0 and lD >= 0),
}

# ---------------- R2 flat-FRW floor
wb, sig_wb = 0.02233, 0.00036       # Mossa et al. 2020, Nature 587, 210 (LUNA D(p,gamma) + D/H), standard BBN
mb = 0.938                            # GeV (protocol convention); old 1% upward allowance
mstar_GeV = mstar_closed * v


def Hfloor(wb_, mb_):
    wh = mstar_GeV * r / mb_ * wb_
    return 100.0 * math.sqrt(wb_ + wh), wh


H_nom, wh_nom = Hfloor(wb, mb)
H_cons, wh_cons = Hfloor(wb - 3 * sig_wb, mb * 1.01)
H_amu, _ = Hfloor(wb - 3 * sig_wb, 0.9315 * 1.01)
res["R2_frw_floor"] = {"omega_h_min_nominal": wh_nom, "Hfloor_nominal": H_nom,
                       "omega_h_min_conservative": wh_cons, "Hfloor_conservative": H_cons,
                       "Hfloor_conservative_with_mb_0.9315_GeV": H_amu,
                       "dHfloor_from_one_sigma_omega_b": 0.5 * H_nom * sig_wb / wb}


# ---------------- R3 distance envelope
def DA_max(z, H, omega_L_ratio=0.0):
    """Upper bound on D_A in flat FRW with H(z)^2 >= H^2 [(1+z)^3 + omega_L_ratio], omega_L_ratio = omega_L/omega_floor (may be < 0)."""
    if omega_L_ratio == 0.0:
        return 2 * c_kms / H / (1 + z) * (1 - (1 + z) ** -0.5)
    val, _ = quad(lambda x: 1.0 / math.sqrt((1 + x) ** 3 + omega_L_ratio), 0, z, epsabs=0, epsrel=1e-12)
    return c_kms / H * val / (1 + z)


def zcos(v_obs, v_pec, convention="optical"):
    z_obs = v_obs / c_kms if convention == "optical" else (v_obs / c_kms) / (1 - v_obs / c_kms)
    b = v_pec / c_kms
    z_pec = math.sqrt((1 + b) / (1 - b)) - 1
    return (1 + z_obs) / (1 + z_pec) - 1


tab = json.loads((HERE.parent / "data" / "pesce2020_table1_claude_transcription.json").read_text())["rows"]
rows = []
for t in tab:
    Dlow = t["D_med"] - 3 * t["D_minus"]
    vhi = t["v"] + 3 * t["sigma_v"]
    out = {"galaxy": t["galaxy"], "D_med": t["D_med"], "D_low_corner": Dlow}
    for lab, vp in (("primary_vpec_-750", -750.0), ("stress_vpec_-1500", -1500.0)):
        z = zcos(vhi, vp)
        out[lab] = {"z_cos": z, "DA_max": DA_max(z, H_cons), "mismatch": bool(Dlow > DA_max(z, H_cons))}
    z_nom = t["v"] / c_kms
    out["nominal_ratio_Dmed_over_DAmax_nominal_floor"] = t["D_med"] / DA_max(z_nom, H_nom)
    z_radio = zcos(vhi, -1500.0, "radio")
    out["stress_if_velocity_were_radio_convention"] = {"DA_max": DA_max(z_radio, H_cons), "mismatch": bool(Dlow > DA_max(z_radio, H_cons))}
    rows.append(out)


def product_bound(row_key, vp):
    """Largest mX*r (GeV, at fixed couplings so m* = mu mX) for which every row passes, at the conservative corners."""
    def worst_margin(mXr):
        wh = (mstar_closed / g) * mXr / (mb * 1.01) * (wb - 3 * sig_wb)
        H = 100 * math.sqrt(wb - 3 * sig_wb + wh)
        return min(DA_max(zcos(t["v"] + 3 * t["sigma_v"], vp), H) - (t["D_med"] - 3 * t["D_minus"]) for t in tab)
    return brentq(worst_margin, 1.0, 1000.0, xtol=1e-10)


res["R3_distance_envelope"] = {"Hfloor_used": H_cons, "rows": rows,
                               "mismatching_rows_primary": [x["galaxy"] for x in rows if x["primary_vpec_-750"]["mismatch"]],
                               "mismatching_rows_stress": [x["galaxy"] for x in rows if x["stress_vpec_-1500"]["mismatch"]],
                               "necessary_bound_mX_r_GeV_primary": product_bound("primary", -750.0),
                               "necessary_bound_mX_r_GeV_stress": product_bound("stress", -1500.0)}

# ---------------- R4 vacuum-energy scope
GeV4_per_rho_crit100 = 1.0 / (1.87834e-26 * 100 ** 2 / (3.0856775814913673e22 ** 2) * 0)  if False else None
# rho_crit,100 = 3 H100^2 / (8 pi G) in GeV^4: H100 = 100 km/s/Mpc = 2.1331e-42 GeV ; M_P(reduced) = 2.435e18 GeV
H100_GeV = 100.0 / 3.0856775814913673e19 * 6.582119569e-25   # (km/s/Mpc -> 1/s) * hbar[GeV s]
MPred = 2.435323e18
rho_c100 = 3 * H100_GeV ** 2 * MPred ** 2
omega_vac_new = -1.22652632365e-6 * v ** 4 / rho_c100
crit = {}
for name in ("NGC 5765b", "UGC 3789"):
    t = next(x for x in tab if x["galaxy"] == name)
    Dlow = t["D_med"] - 3 * t["D_minus"]
    for lab, vp in (("primary", -750.0), ("stress", -1500.0)):
        z = zcos(t["v"] + 3 * t["sigma_v"], vp)
        f = lambda xL: DA_max(z, H_cons, xL) - Dlow
        if f(0.0) >= 0:
            crit[f"{name} {lab}"] = None
        else:
            xr = brentq(f, -0.999999, 0.0, xtol=1e-12)
            wfloor = (H_cons / 100) ** 2
            crit[f"{name} {lab}"] = {"omega_L_over_omega_floor": xr, "omega_L_critical": xr * wfloor}
res["R4_vacuum_energy_scope"] = {
    "rho_crit_100_GeV4": rho_c100,
    "benchmark_one_loop_new_scalar_vacuum_energy_GeV4": -1.22652632365e-6 * v ** 4,
    "same_in_omega_units": omega_vac_new,
    "critical_total_omega_Lambda_below_which_the_row_no_longer_mismatches": crit,
    "note": "Mismatch survives only if the total effective vacuum density is above omega_L_critical (negative values allowed down to it).",
}


# ---------------- R5 one-loop Coleman-Weinberg estimate of the homogeneous floor (units v = 1, MSbar, mu = mX = g v)
mu2 = g ** 2


def masses2(phi, s):
    """Field-dependent mass^2 eigenvalues of (phi, q1, q2, D) at background X = s (real), D = 0."""
    q1 = math.sqrt(2) * s
    m_pp = lphi * (3 * phi ** 2 - 1) + 2 * g ** 2 * s ** 2
    m_p1 = 2 * g ** 2 * phi * q1
    m_11 = g ** 2 * phi ** 2 + 6 * lX * s ** 2
    m_22 = g ** 2 * phi ** 2 + 2 * lX * s ** 2
    m_DD = mD0_over_v ** 2 + hD ** 2 * phi ** 2
    M = np.array([[m_pp, m_p1], [m_p1, m_11]])
    ev = np.linalg.eigvalsh(M)
    return [ev[0], ev[1], m_22, m_DD]


def V1(phi, s):
    tot, neg = 0.0, False
    for m2 in masses2(phi, s):
        if m2 < 0:
            neg = True
        a = abs(m2)
        if a > 0:
            tot += a ** 2 * (math.log(a / mu2) - 1.5)  # real part for m2 < 0
    return tot / (64 * math.pi ** 2), neg


def Vtree(phi, s):
    return U(phi, math.sqrt(2) * s, 0.0, 0.0)


def Veff(phi, s):
    return Vtree(phi, s) + V1(phi, s)[0]


ovac = minimize(lambda x: Veff(x[0], 0.0), [1.0], method="Nelder-Mead", options={"xatol": 1e-14, "fatol": 1e-20})
phi_vac = float(ovac.x[0]); V_vac = float(ovac.fun)
best1 = (1e9, None)
for p0 in np.linspace(0, 1.1, 12):
    for ls0 in np.linspace(-3, 1.5, 10):
        f = lambda x: (Veff(x[0], math.exp(x[1])) - V_vac) / math.exp(2 * x[1])
        o = minimize(f, [p0, ls0], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-18, "maxiter": 40000})
        if o.fun < best1[0]:
            best1 = (o.fun, o.x)
mstar_1l = math.sqrt(best1[0])
phi_at, s_at = float(best1[1][0]), float(math.exp(best1[1][1]))
neg_at = V1(phi_at, s_at)[1]
res["R5_one_loop_floor_estimate"] = {
    "one_loop_vacuum_phi_over_v": phi_vac, "archived_value": 1.00012193740,
    "vacuum_energy_one_loop_over_v4": V_vac,
    "archived_new_scalar_vacuum_energy_over_v4": -1.22652632365e-6,
    "mstar_tree_over_v": mstar_closed, "mstar_one_loop_over_v": mstar_1l,
    "relative_change": mstar_1l / mstar_closed - 1,
    "location_phi_over_v": phi_at, "location_s_over_v": s_at,
    "negative_mass_squared_at_location": neg_at,
    "note": "Derivative-expansion estimate: one-loop effective potential, vacuum energy subtracted (it gravitates through the separately "
            "matched Lambda), wavefunction renormalization neglected. Not a quantum energy inequality.",
}
# fractional reduction of m* that would remove each mismatch
need = {}
for lab, vp in (("primary", -750.0), ("stress", -1500.0)):
    bound = res["R3_distance_envelope"][f"necessary_bound_mX_r_GeV_{lab}"]
    need[lab] = 1 - bound / (mX * r)
res["R5_one_loop_floor_estimate"]["fractional_reduction_of_mstar_needed_to_remove_all_mismatches"] = need

# ---------------- R6 order-of-magnitude preview of the proposed next calculation: Z2 domain walls
# Leading high-T thermal mass of phi from bosonic loops: V_T ~ (T^2/24) * sum_i m_i^2(phi), coefficient of phi^2 = (3 lphi + 2 g^2 + hD^2)/24
Tc_over_v = math.sqrt(12 * lphi / (3 * lphi + 2 * g ** 2 + hD ** 2))
sigma = (2 * math.sqrt(2) / 3) * math.sqrt(lphi) * v ** 3           # GeV^3, kink of (lphi/4)(phi^2-v^2)^2 with X = D = 0
MPred_ = 2.435323e18
# a scaling wall network has rho_wall ~ A sigma H with A ~ 1; it dominates radiation (rho ~ 3 H^2 MP^2) when H ~ A sigma / (3 MP^2)
H_dom = sigma / (3 * MPred_ ** 2)
gstar = 3.36                                                        # photons + neutrinos after e+e- annihilation (illustrative)
T_dom = (H_dom * MPred_ * math.sqrt(90 / (math.pi ** 2 * gstar))) ** 0.5
G_inv = (1.220890e19) ** 2                                          # 1/G in GeV^2
t0_GeV_inv = 4.35e17 / 6.582119569e-25
sigma_ZKO = 1e-5 * G_inv / t0_GeV_inv                               # G sigma t0 < 1e-5 (CMB isotropy)
res["R6_domain_wall_preview"] = {
    "status": "order-of-magnitude preview of the PROPOSED next calculation, not a result",
    "Tc_over_v_leading_high_T": Tc_over_v, "Tc_GeV": Tc_over_v * v,
    "wall_tension_GeV3": sigma, "wall_tension_cuberoot_GeV": sigma ** (1 / 3),
    "domination_temperature_GeV_if_scaling_with_A_1": T_dom,
    "ZKO_CMB_bound_sigma_GeV3": sigma_ZKO, "ZKO_bound_cuberoot_MeV": 1e3 * sigma_ZKO ** (1 / 3),
    "ratio_sigma_over_bound": sigma / sigma_ZKO,
    "note": "Exact Z2 (phi -> -phi) of the unchanged action: walls are stable. If the hidden sector was ever hotter than ~Tc, a network forms; "
            "the estimate says it would dominate well before recombination and violate CMB isotropy by many orders. Quantities to be "
            "computed properly: Tc with the full one-loop thermal potential and X/D dynamics; whether symmetric-phase core interiors "
            "(phi ~ 0 inside the charged objects) can seed -v domains on fragmentation; the network's scaling constant.",
}

res["provenance"] = {"script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "python": sys.version.split()[0],
                     "numpy": np.__version__, "scipy": scipy.__version__, "wall_seconds": round(time.time() - t0, 2)}


def rnd(o):
    if isinstance(o, float):
        return float(f"{o:.12g}")
    if isinstance(o, dict):
        return {k: rnd(x) for k, x in o.items()}
    if isinstance(o, list):
        return [rnd(x) for x in o]
    return o


(OUT / "review_results.json").write_text(json.dumps(rnd(res), indent=1))
R1, R2, R3, R4, R5 = (res[k] for k in ("R1_classical_floor", "R2_frw_floor", "R3_distance_envelope", "R4_vacuum_energy_scope", "R5_one_loop_floor_estimate"))
print("R1 m*/v closed %.12f direct %.12f; mu = %.12f; m* = %.8f GeV; random violations %d of %d; saturation residual %.2e" % (
    R1["mstar_over_v_closed_form"], R1["mstar_over_v_direct_minimization"], R1["mu_mstar_over_mX"], R1["mstar_GeV"], R1["pointwise_random_test"]["violations"], 4 * N, R1["saturation_check_rho_minus_mstar_n"]))
print("R2 omega_h,min %.10f Hfloor %.8f; conservative omega_h %.10f Hfloor %.8f; with mb=0.9315: %.4f" % (R2["omega_h_min_nominal"], R2["Hfloor_nominal"], R2["omega_h_min_conservative"], R2["Hfloor_conservative"], R2["Hfloor_conservative_with_mb_0.9315_GeV"]))
for x in R3["rows"]:
    print("R3 %-13s Dlow %7.2f | DAmax primary %7.2f %s | stress %7.2f %s | radio-stress %7.2f | nominal Dmed/DAmax %.3f" % (
        x["galaxy"], x["D_low_corner"], x["primary_vpec_-750"]["DA_max"], "MISMATCH" if x["primary_vpec_-750"]["mismatch"] else "ok",
        x["stress_vpec_-1500"]["DA_max"], "MISMATCH" if x["stress_vpec_-1500"]["mismatch"] else "ok", x["stress_if_velocity_were_radio_convention"]["DA_max"], x["nominal_ratio_Dmed_over_DAmax_nominal_floor"]))
print("R3 necessary bound mX*r: primary %.4f GeV, stress %.4f GeV" % (R3["necessary_bound_mX_r_GeV_primary"], R3["necessary_bound_mX_r_GeV_stress"]))
print("R4 rho_c100 %.4e GeV^4; benchmark one-loop vacuum %.4e GeV^4 = %.3e in omega units; critical omega_L: %s" % (R4["rho_crit_100_GeV4"], R4["benchmark_one_loop_new_scalar_vacuum_energy_GeV4"], R4["same_in_omega_units"], {k: (None if v_ is None else round(v_["omega_L_critical"], 5)) for k, v_ in R4["critical_total_omega_Lambda_below_which_the_row_no_longer_mismatches"].items()}))
print("R5 one-loop vacuum phi/v %.11f (archived 1.00012193740); V_vac/v^4 %.6e; m* tree %.10f one-loop %.10f rel %.3e at phi %.4f s %.4f negm2 %s; needed reduction %s" % (
    R5["one_loop_vacuum_phi_over_v"], R5["vacuum_energy_one_loop_over_v4"], R5["mstar_tree_over_v"], R5["mstar_one_loop_over_v"], R5["relative_change"], R5["location_phi_over_v"], R5["location_s_over_v"], R5["negative_mass_squared_at_location"], {k: round(v_, 4) for k, v_ in R5["fractional_reduction_of_mstar_needed_to_remove_all_mismatches"].items()}))
R6 = res["R6_domain_wall_preview"]
print("R6 preview: Tc ~ %.3f v = %.0f GeV; sigma^(1/3) = %.1f GeV; walls would dominate at T ~ %.2e GeV; ZKO bound sigma^(1/3) < %.2f MeV; ratio %.1e" % (
    R6["Tc_over_v_leading_high_T"], R6["Tc_GeV"], R6["wall_tension_cuberoot_GeV"], R6["domination_temperature_GeV_if_scaling_with_A_1"], R6["ZKO_bound_cuberoot_MeV"], R6["ratio_sigma_over_bound"]))
print("wall %.1fs" % (time.time() - t0))
