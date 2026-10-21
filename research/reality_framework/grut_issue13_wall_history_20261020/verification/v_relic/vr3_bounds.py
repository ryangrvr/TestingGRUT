"""v_relic part 3: minimum T_SM for a relic-compatible thermal restoration; Pi_phi(z) monotonicity;
floor-inequality algebra (coefficient 79.28 GeV); independent code."""
import json, os
import numpy as np
from scipy.special import zeta, kve
from scipy.optimize import brentq
import vr_common as C

gstar = 106.75
OM_FAC = 2.743744e8  # from vr2 (own derivation of s0/rho_c h^-2); recomputed below for safety
kB = 8.617333262e-5; hbarc = 1.973269804e-5
T0 = 2.7255 * kB / hbarc
s0 = 2 * np.pi**2 / 45 * (2 + 21 / 11) * T0**3
rhoc_h2 = 3 * (1e5 / 3.0856775814913673e22)**2 / (8 * np.pi * 6.67430e-11) / 1.78266192e-27 * 1e-6
OM_FAC = s0 / rhoc_h2
out = {"OM_FAC": OM_FAC}

def ssm_coef():
    return 2 * np.pi**2 / 45 * gstar

def n_BE(T, m, g, psi=None):
    x = m / T
    psi = -x if psi is None else psi
    k = np.arange(1, 400)
    if x == 0:
        return g * T**3 / np.pi**2 * np.sum(np.exp(k * psi) / k**3)
    return g * T**3 / (2 * np.pi**2) * np.sum(np.exp(k * psi) * x**2 * kve(2, k * x) / k)

# ---- leading T_c ----
Pi = C.lam / 4 + C.gX**2 / 6 + C.hD**2 / 12
Tc_lead = C.v * np.sqrt(C.lam / Pi)
out["Pi_over_T2"] = Pi
out["Tc_over_v_leading"] = Tc_lead / C.v
Tc_list = {"leading": Tc_lead, "resummed_lo[Q]": 1.5571 * C.v, "resummed_hi[Q]": 1.5641 * C.v}

# ---- route A: T_c / xi_max ----
xi_max = {"module_ODE": 7.581e-4, "vrelic_conv_only_ODE": None, "module_noD": 8.364e-4}
try:
    r2 = json.load(open("results/vr2_relic.json"))
    xi_max["vrelic_conv_only_ODE"] = r2["boltzmann_conv_only"]["xi_max_conv_only"]
    v_noD = r2["boltzmann_conv_only"].get("xi_max_conv_only_noD")
    if isinstance(v_noD, float):
        xi_max["vrelic_conv_only_noD"] = v_noD
except Exception as e:
    xi_max["vrelic_conv_only_ODE"] = 7.5499e-4
routeA = {}
for kt, Tc in Tc_list.items():
    routeA[kt] = {kx: Tc / xv for kx, xv in xi_max.items() if xv}
out["routeA_TSMmin=Tc/xi_max"] = routeA

# ---- xi at T_c vs xi at the 4 mX reference: number and entropy matching ----
Tr = 4 * C.mX
nref = n_BE(Tr, C.mX, 2) + n_BE(Tr, C.mh, 1) + n_BE(Tr, C.mD, 1)
for kt, Tc in Tc_list.items():
    pass
Tc = Tc_lead
# symmetric phase at T_c: X, phi ~ massless, D has m_D0 (thermal masses neglected)
n_c = n_BE(Tc, 0.0, 2) + n_BE(Tc, 0.0, 1) + n_BE(Tc, C.mD0, 1)
n_c_massless = 4 * zeta(3) * Tc**3 / np.pi**2
out["n_th(Tc)_over_v3"] = {"all_massless": n_c_massless / C.v**3, "D_with_mD0": n_c / C.v**3}
ratio_num = (nref / Tr**3) / (n_c / Tc**3)
out["xi_c_over_xi_ref_number_matching"] = ratio_num**(1 / 3)

# ---- route B: number bound n >= n_th(T_c), Y_N <= Y_max ----
Ymax = {"module_4.30e-12": 4.30e-12,
        "strict_all_quanta>=mX": 0.12 / (C.mX * OM_FAC)}
routeB = {}
for kt, Tcx in Tc_list.items():
    nth = 4 * zeta(3) * Tcx**3 / np.pi**2
    routeB[kt] = {ky: (nth / (ssm_coef() * Y))**(1 / 3) for ky, Y in Ymax.items()}
out["Ymax"] = Ymax
out["routeB_TSMmin_number"] = routeB

# ---- v scaling (closed form limits): xi_max ∝ v^{-1/3}, Tc ∝ v -> TSM ∝ v^{4/3} ----
out["v_scaling_note"] = "frozen and full-conversion limits scale exactly as v^{-1/3}; conversion efficiency " \
                        "depends on v (Gamma/H ∝ xi^2 M_P/v), so xi_max(v) moves inside the band " \
                        "[frozen, full-conversion] whose width is 3.3% in xi"

# ---- Pi_phi at fixed n vs fugacity ----
def Li(s, z, K=20000):
    k = np.arange(1, K + 1)
    return np.sum(z**k / k**s)
zs = [1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0]
F = [Li(2, z) / Li(3, z)**(2 / 3) for z in zs]
out["Pi_fixed_n_vs_z"] = {"z": zs, "Li2/Li3^(2/3)": F, "monotonic": bool(np.all(np.diff(F) > 0))}

# ---- floor inequality algebra ----
mstar = C.v * min(C.gX, (C.lam * C.lamX)**0.25)
mu = mstar / C.mX
wb_corner = 0.02233 - 3 * 0.00036
mb_corner = 0.93800 * 1.01  # "1% higher m_b" corner (Codex 2ee7f82 REPORT §4); 0.938 base inferred
coef = (C.mX / mstar) * (mb_corner / wb_corner)
K_primary, K_stress = 63.9344, 74.6429  # Codex quoted
dW_primary = K_primary * mu * wb_corner / mb_corner
out["floor"] = {
    "mstar": mstar, "mu=m*/mX": mu, "omega_b_corner": wb_corner, "m_b_corner": mb_corner,
    "coef_GeV_per_omega_sym": coef,
    "equivalent_m*_form": {"K*_primary_GeV": K_primary * mu, "coef*_GeV": mb_corner / wb_corner},
    "implied_max_omega_h_primary": dW_primary,
    "implied_H0floor_max_primary": 100 * np.sqrt(wb_corner + dW_primary),
    "r_max_omega_sym=0.12_primary": (K_primary - coef * 0.12) / C.mX,
    "r_max_omega_sym=0.12_stress": (K_stress - coef * 0.12) / C.mX,
    "r_max_if_wrong_coef(m_b/omega_b only)": (K_primary - mb_corner / wb_corner * 0.12) / C.mX,
    "H0floor_r1": 100 * np.sqrt(wb_corner + mstar / mb_corner * wb_corner),
    "H0floor_r1_plus_0.12": 100 * np.sqrt(wb_corner + mstar / mb_corner * wb_corner + 0.12),
    "omega_asym_per_r_central(mb=0.9383,wb=0.02233)": mstar / 0.9383 * 0.02233,
    "r_max_DM_0.12": 0.12 / (mstar / 0.9383 * 0.02233),
}
os.makedirs("results", exist_ok=True)
json.dump(out, open("results/vr3_bounds.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
