"""Item 1: which interactions couple the hidden scalars to the Standard Model inside the unchanged action.

(a) One-loop scalar beta tensor for 8 real scalars (4 Higgs components + phi, q1, q2, D):
    16 pi^2 beta_abcd = lambda_abef lambda_efcd + lambda_acef lambda_efbd + lambda_adef lambda_efbc
    (+ gauge/Yukawa anomalous-dimension terms, which are diagonal: gamma_a lambda_abcd + perms, and
    pure-gauge/Yukawa terms that only carry Higgs indices).
    Checks: mixed Higgs-hidden entries vanish identically when the three Higgs portals vanish;
    single-singlet portal beta reproduces Lebedev arXiv:1203.0156 eq.(15); hidden betas reproduce
    the archived rg_beta_functions.txt.
(b) Running of the nonminimal couplings xi(mu) from zero at mu_0 and their size during inflation,
    matter-like reheating and the radiation era (trace anomaly of the SM plasma).
(c) Dimensional estimates of graviton-induced local portals.
"""
import itertools
import numpy as np
from common import (M_P, V_BENCH, LAM_PHI, LAM_X, LAM_D, C_PORT, D_PORT, MU0, G_STAR_SM,
                    beta_xi_coeffs, xi_running, dump, rho_rad)

NF = 8  # h1..h4, phi, q1, q2, D
IH = [0, 1, 2, 3]
IPHI, IQ1, IQ2, ID = 4, 5, 6, 7


def quartic(x, p):
    h2 = x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2           # = 2 H^dag H
    phi, q1, q2, D = x[4], x[5], x[6], x[7]
    R2 = q1**2 + q2**2                                   # = 2|X|^2
    V = (p["lh"] / 4) * h2**2
    V += (p["lp"] / 4) * phi**4 + (p["lx"] / 4) * R2**2 + (p["ld"] / 4) * D**4
    V += (p["c"] / 2) * phi**2 * R2 + (p["d"] / 2) * phi**2 * D**2 + (p["al"] / 2) * R2 * D**2
    # Higgs portals (zero at mu_0 by matching); Lebedev-like normalisation lambda/4 h^2 s^2
    V += (p["hp"] / 4) * h2 * phi**2 + (p["hx"] / 4) * h2 * R2 + (p["hd"] / 4) * h2 * D**2
    return V


def lam_tensor(p):
    """Exact 4th-derivative tensor via the polarisation identity (Q homogeneous quartic)."""
    L = np.zeros((NF,) * 4)
    E = np.eye(NF)
    signs = list(itertools.product([1, -1], repeat=4))
    for a, b, c, d in itertools.combinations_with_replacement(range(NF), 4):
        s = 0.0
        for e in signs:
            s += e[0] * e[1] * e[2] * e[3] * quartic(e[0] * E[a] + e[1] * E[b] + e[2] * E[c] + e[3] * E[d], p)
        val = s / 16.0
        for perm in set(itertools.permutations((a, b, c, d))):
            L[perm] = val
    return L


def beta_tensor(L):
    t1 = np.einsum("abef,efcd->abcd", L, L)
    t2 = np.einsum("acef,efbd->abcd", L, L)
    t3 = np.einsum("adef,efbc->abcd", L, L)
    return (t1 + t2 + t3) / (16 * np.pi**2)


def part_a():
    out = {}
    p0 = dict(lh=0.13, lp=LAM_PHI, lx=LAM_X, ld=LAM_D, c=C_PORT, d=D_PORT, al=1.3e-8, hp=0.0, hx=0.0, hd=0.0)
    L = lam_tensor(p0)
    B = beta_tensor(L)
    hidden = [IPHI, IQ1, IQ2, ID]
    mixed_max = 0.0
    for idx in itertools.product(range(NF), repeat=4):
        nh = sum(1 for i in idx if i in IH)
        if 0 < nh < 4:
            mixed_max = max(mixed_max, abs(B[idx]))
    nmix_lam = 0
    for idx in itertools.product(range(NF), repeat=4):
        nh = sum(1 for i in idx if i in IH)
        if 0 < nh < 4 and abs(L[idx]) > 1e-14:
            nmix_lam += 1
    out["mixed_lambda_entries_above_1e-14_at_zero_portals"] = nmix_lam
    out["mixed_beta_max_at_zero_portals"] = mixed_max
    out["mixed_beta_below_roundoff_1e-18"] = bool(mixed_max < 1e-18)
    out["mixed_beta_note"] = "nonzero residue is floating roundoff of the polarisation identity; structurally zero"
    # hidden betas vs archived rg_beta_functions.txt (B = 16 pi^2 beta)
    lp, lx, ld, c, d, al = p0["lp"], p0["lx"], p0["ld"], p0["c"], p0["d"], p0["al"]
    k = 16 * np.pi**2
    checks = {
        "B_lp": (k * B[IPHI, IPHI, IPHI, IPHI] / 6, 18 * lp**2 + 4 * c**2 + 2 * d**2),
        "B_lx": (k * B[IQ1, IQ1, IQ1, IQ1] / 6, 20 * lx**2 + 2 * c**2 + 2 * al**2),
        "B_ld": (k * B[ID, ID, ID, ID] / 6, 18 * ld**2 + 2 * d**2 + 4 * al**2),
        "B_c": (k * B[IPHI, IPHI, IQ1, IQ1] / 2, 6 * lp * c + 8 * lx * c + 8 * c**2 + 2 * d * al),
        "B_d": (k * B[IPHI, IPHI, ID, ID] / 2, 6 * lp * d + 6 * ld * d + 4 * c * al + 8 * d**2),
        "B_alpha": (k * B[IQ1, IQ1, ID, ID] / 2, 8 * lx * al + 6 * ld * al + 2 * c * d + 8 * al**2),
    }
    out["hidden_beta_vs_archive"] = {kk: {"tensor": v[0], "archive_formula": v[1],
                                          "rel_diff": abs(v[0] - v[1]) / max(abs(v[1]), 1e-300)}
                                     for kk, v in checks.items()}
    out["hidden_beta_archive_pass"] = bool(all(abs(v[0] - v[1]) <= 1e-12 * max(abs(v[1]), 1e-30) for v in checks.values()))
    # Lebedev single-singlet check: V = lh/4 h^4 + lhs/4 h^2 s^2 + ls/4 s^4  (s = phi, others decoupled)
    pl = dict(lh=0.13, lp=0.02, lx=0.0, ld=0.0, c=0.0, d=0.0, al=0.0, hp=0.05, hx=0.0, hd=0.0)
    Bl = beta_tensor(lam_tensor(pl)) * k
    lebedev = 4 * pl["hp"]**2 + 12 * pl["lh"] * pl["hp"] + 6 * pl["lp"] * pl["hp"]
    out["lebedev_portal_check"] = {"tensor_16pi2_beta_lhs": Bl[0, 0, IPHI, IPHI], "lebedev_eq15_scalar_part": lebedev,
                                   "rel_diff": abs(Bl[0, 0, IPHI, IPHI] - lebedev) / lebedev,
                                   "note": "gauge/Yukawa part of Lebedev eq.(15), 6 yt^2 - (3/2)(3g^2+g'^2), multiplies lambda_hs"}
    # linearised portal running at zero portals: derivative of mixed betas wrt each portal
    jac = np.zeros((3, 3))
    names = ["hp", "hx", "hd"]
    ent = [(0, 0, IPHI, IPHI), (0, 0, IQ1, IQ1), (0, 0, ID, ID)]
    eps = 1e-6
    for j, nm in enumerate(names):
        pp = dict(p0)
        pp[nm] = eps
        Bp = beta_tensor(lam_tensor(pp)) * k
        for i, e in enumerate(ent):
            jac[i, j] = Bp[e] / eps
    out["portal_anomalous_matrix_16pi2_at_mu0_scalar_part"] = jac
    out["portal_anomalous_matrix_note"] = ("rows: d(lambda_H phi, lambda_HX, lambda_HD)/dlnmu x 16pi^2; columns: same portals; "
                                           "beta(portals)=M.portals+O(portal^2); zero is a fixed subspace")
    return out


def trace_anomaly_ratio(T):
    """(rho - 3p)/T^4 for the SM plasma from the O(g^2) pressure of each gauge factor,
    p2 = -(g^2 T^4/144)[d_A C_A + (5/2) sum_Dirac d_R C_R]; rho-3p = T d(p/T^4)/dT T^4.
    Couplings: one-loop SM running from M_Z, mu = 2 pi T (scale choice; labelled estimate).
    Scalars, Yukawa and quartic contributions omitted (estimated <~ 50%)."""
    MZ = 91.1876
    a3, a2, a1 = 0.1179, 0.0338, 0.0102          # alpha_s, alpha_2, alpha_Y (non-GUT) at M_Z, rounded
    b = {"3": 7.0, "2": 19.0 / 6.0, "1": -41.0 / 6.0}   # 16pi^2 dg/dlnmu = -b g^3
    groups = {"3": (a3, 8, 3, 24.0), "2": (a2, 3, 2, 9.0), "1": (a1, 1, 0, 5.0)}
    mu = 2 * np.pi * T
    total = 0.0
    parts = {}
    for gname, (a0, dA, CA, SF) in groups.items():
        inv = 1.0 / a0 + b[gname] / (2 * np.pi) * np.log(mu / MZ)
        al = 1.0 / inv
        g4 = (4 * np.pi * al)**2
        val = (2.0 / 144.0) * (dA * CA + 2.5 * SF) * b[gname] * g4 / (16 * np.pi**2)
        parts[gname] = val
        total += val
    return total, parts


def part_b():
    out = {}
    A = beta_xi_coeffs()
    b0 = A @ (np.zeros(3) - 1.0 / 6.0) / (16 * np.pi**2)
    out["beta_xi_at_mu0_per_efold"] = {"xi_phi": b0[0], "xi_X": b0[1], "xi_D": b0[2]}
    out["sign_note"] = ("beta_xi < 0 at xi=0: xi(mu) < 0 for mu > mu_0 (anti-restoring curvature mass), "
                        "xi(mu) > 0 for mu < mu_0")
    mus = [10.0, 100.0, 1e3, 1e4, 1e6, 1e10, 1e13, 1e16]
    out["xi_of_mu"] = {f"{m:.0e}": dict(zip(["xi_phi", "xi_X", "xi_D"], xi_running(m))) for m in mus}
    lv2 = LAM_PHI * V_BENCH**2
    # radiation era
    rows = []
    Ts = np.logspace(2, 16, 57)
    Tstar_abs = None
    for T in Ts:
        tr, parts = trace_anomaly_ratio(T)
        rho = rho_rad(T, G_STAR_SM)
        R = tr * T**4 / M_P**2
        xi = xi_running(T)[0]
        rows.append({"T": T, "rho_minus_3p_over_T4": tr, "one_minus_3w": tr * T**4 / rho, "R_GeV2": R,
                     "xi_phi_mu_T": xi, "xi_R_over_lambda_v2": xi * R / lv2})
        if Tstar_abs is None and abs(xi * R) > lv2:
            Tstar_abs = T
    out["radiation_era_rows"] = rows[::4]
    # find crossing |xi R| = lambda v^2 precisely
    f = lambda lT: np.log(abs(xi_running(10**lT)[0] * trace_anomaly_ratio(10**lT)[0] * 10**(4 * lT) / M_P**2) / lv2)
    from scipy.optimize import brentq
    lTs = brentq(f, 6, 18)
    out["T_where_abs_xiR_equals_lambda_v2_radiation"] = 10**lTs
    out["T_where_abs_xiR_equals_lambda_v2_if_trace_x10"] = 10**brentq(lambda lT: f(lT) + np.log(10.0), 6, 18)
    out["radiation_verdict"] = ("xi_phi(mu~T) < 0, so xi R phi^2 with R=(rho-3p)/M_P^2 > 0 deepens the broken minimum; "
                                "its magnitude exceeds lambda_phi v^2 only above T_* (anti-restoring there).")
    # inflation and matter-like reheating (outside-action preparation data: H)
    infl = []
    for H in [1e2, 3e2, 1e3, 3e3, 1e4, 1e6, 1e10, 1e13]:
        xi = xi_running(H)
        infl.append({"H": H, "xi_phi(H)": xi[0], "xi_X(H)": xi[1],
                     "xi_phi*12H^2/(lambda v^2) (de Sitter)": xi[0] * 12 * H**2 / lv2,
                     "xi_phi*3H^2/(lambda v^2) (matter-like)": xi[0] * 3 * H**2 / lv2,
                     "m_X^2(phi=v)+xi_X 12 H^2 [GeV^2]": C_PORT * V_BENCH**2 + xi[1] * 12 * H**2})
    out["inflation_reheating_rows"] = infl
    # H above which X becomes tachyonic at phi=+v during de Sitter (xi_X(mu=H)<0)
    g = lambda lH: C_PORT * V_BENCH**2 + xi_running(10**lH)[1] * 12 * 10**(2 * lH)
    out["H_X_tachyonic_at_phi_v_dS"] = 10**brentq(g, 2.5, 8)
    return out


def part_c():
    """Graviton-induced local portal estimates (labelled dimensional estimates, not computed diagrams)."""
    mphi = np.sqrt(2 * LAM_PHI) * V_BENCH
    mh = 125.25
    out = {}
    # (i) one-loop graviton diagrams with one graviton line between sectors: local H^dag H phi^2
    #     counterterm needs two powers of a mass (Rodigast-Schuster 0908.2422: gravitational running only for massive fields)
    est1 = max(mphi, mh)**2 / (16 * np.pi**2 * M_P**2)
    out["delta_lambda_portal_one_loop_graviton_estimate"] = est1
    # (ii) classical: eliminate R with Einstein eq. in xi_H R H^dag H -> xi_H T_hidden H^dag H / M_P^2
    xiH_bound = 0.1
    est2 = xiH_bound * mphi**2 / M_P**2
    out["xi_H_bound_used"] = xiH_bound
    out["effective_portal_from_xi_H_estimate"] = est2
    lam_eff = max(est1, est2)
    # thermal mass from such a portal: Pi ~ lam_eff * (4 Higgs dof) T^2/12 -> T needed for Pi = lambda v^2
    Tneed = np.sqrt(LAM_PHI * V_BENCH**2 / (lam_eff * 4 / 12.0))
    out["T_needed_for_portal_thermal_mass_to_restore"] = Tneed
    out["T_needed_over_M_P"] = Tneed / M_P
    out["verdict"] = "local graviton-induced portals are negligible below M_P; sectors couple through the metric only"
    return out


def run():
    res = {"item": "1 couplings between sectors", "a_portal_beta": part_a(), "b_nonminimal": part_b(),
           "c_graviton_portal_estimates": part_c()}
    res["coupling_verdict"] = ("Inside the unchanged action with zero Higgs portals at mu_0 the portal subspace is "
                               "a fixed subspace of the non-gravitational RG (exact tensor check: mixed betas = 0); "
                               "graviton-induced local portals are <~1e-34; xi(mu) R phi^2 is anti-restoring for mu>mu_0. "
                               "The hidden and visible sectors interact only through the metric (gravity), incl. "
                               "s-channel graviton exchange and the shared FRW background.")
    dump("c1_couplings.json", res)
    return res


if __name__ == "__main__":
    r = run()
    a = r["a_portal_beta"]
    print("mixed beta max:", a["mixed_beta_max_at_zero_portals"], "archive pass:", a["hidden_beta_archive_pass"])
    print("lebedev:", a["lebedev_portal_check"])
    b = r["b_nonminimal"]
    print("beta_xi:", b["beta_xi_at_mu0_per_efold"])
    print("T* rad:", b["T_where_abs_xiR_equals_lambda_v2_radiation"], b["T_where_abs_xiR_equals_lambda_v2_if_trace_x10"])
    print("H_X_tach:", b["H_X_tachyonic_at_phi_v_dS"])
    for row in b["radiation_era_rows"]:
        print(row)
    print(r["c_graviton_portal_estimates"])
