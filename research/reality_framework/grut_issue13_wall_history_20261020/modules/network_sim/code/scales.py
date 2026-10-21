"""Unit conversions linking the dimensionless simulations to the unchanged benchmark (derived numbers; the
cosmological inputs g_* and the formation temperature are labelled illustrative/free).
Writes results/scales.json."""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")

v = 1000.0                 # GeV, benchmark
lam = 0.01
m_phi = math.sqrt(2 * lam) * v
sigma = 2 * math.sqrt(2 * lam) / 3 * v ** 3          # GeV^3
delta = 2.0 / m_phi                                  # GeV^-1, tanh half-width
M_P = 2.435323e18                                    # reduced Planck mass, GeV
M_Pl = math.sqrt(8 * math.pi) * M_P                  # non-reduced
hbar_s = 6.582119569e-25                             # GeV s
hbarc_m = 1.973269804e-16                            # GeV m


def H_rad(T, gstar):
    """Hubble rate in a radiation era with g_* (visible + any hidden radiation; free input)."""
    return math.sqrt(math.pi ** 2 * gstar / 90.0) * T * T / M_P


out = {"benchmark": {"v_GeV": v, "lambda_phi": lam, "m_phi_GeV": m_phi, "sigma_GeV3": sigma,
                     "sigma_cuberoot_GeV": sigma ** (1 / 3), "delta_GeV_inv": delta, "delta_m": delta * hbarc_m,
                     "status": "derived from the unchanged action (tree level, X = D = 0)"}}
# horizon / wall-width hierarchy at an illustrative formation temperature (NOT a transition result)
rows = []
for T in (1.547646465 * v, v, 100.0):
    for g in (106.75, 110.75):
        H = H_rad(T, g)
        rows.append({"T_GeV": T, "gstar": g, "H_GeV": H, "Hinv_over_delta": 1 / (H * delta)})
out["hierarchy_at_formation_illustrative"] = {
    "rows": rows,
    "note": "T values are illustrative (1.5476 v is the leading curvature-zero CHECK of BRIEF.md, not a transition "
            "temperature); g_* = 106.75 (SM) or 110.75 (SM + 4 hidden real scalar dof phi, X, D in equilibrium) are free inputs."}
# simulated hierarchies (physical-equation runs, a = tau in m_phi units): H^-1/delta = tau^2/2
out["simulated_hierarchy"] = {"phys2d_tau_end_32": 32 ** 2 / 2, "phys3d_tau_end_12": 12 ** 2 / 2,
                              "phys2d_tau_i_2": 2.0, "note": "H^-1 / delta = a tau / 2 = tau^2 / 2 for a = tau"}
# wall domination time of a scaling network rho_w = A sigma / t in radiation era rho = 3 M_P^2 / (4 t^2)
dom = []
for A in (0.5, 0.8, 1.0):
    t_w = 3 * M_P ** 2 / (4 * A * sigma)    # GeV^-1
    # visible temperature at t_w for g_* = 3.36 (photons + 3 nu, after e+e-; illustrative; hidden radiation free)
    g = 3.36
    T = math.sqrt(1.0 / (2 * t_w) / (math.sqrt(math.pi ** 2 * g / 90.0) / M_P))
    dom.append({"A": A, "t_w_GeV_inv": t_w, "t_w_s": t_w * hbar_s, "T_vis_at_t_w_keV_illustrative": T * 1e6})
out["wall_domination_if_scaling_persists"] = {
    "rows": dom, "formula": "t_w = 3 M_P^2 / (4 A sigma) (radiation era, no backreaction until then)",
    "note": "Conditional: holds only if an unbiased scaling network forms and persists with the tree tension; it is "
            "the validity limit of the no-backreaction simulations, not an observable prediction."}
with open(os.path.join(RES, "scales.json"), "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
