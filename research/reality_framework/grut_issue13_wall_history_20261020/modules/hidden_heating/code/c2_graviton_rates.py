"""Item 2a: s-channel graviton exchange SM SM -> h_munu -> (hidden real scalar pair), computed here.

Amplitude M = T^i_{mu nu} P^{mu nu rho sigma} T^f_{rho sigma} / (M_P^2 s),
P = (eta eta + eta eta - eta eta)/2, vertex -i h_{mu nu} T^{mu nu}/M_P (g = eta + 2h/M_P),
normalisation fixed by the Newtonian limit (checked below).  Mostly-minus metric internally.
Initial states: real scalar (4 Higgs components, xi_H = 0 at mu_0), left-handed Weyl fermion pair
(45 Weyl fields), massless vector pair (12 real vector fields).  Final: one massless real scalar with
nonminimal coupling xi (xi = 0 at mu_0); the 4 hidden real dof (phi, q1, q2, D) each get this rate.
Checks: Dirac completeness, tracelessness of fermion/vector T, comparison with the angular shapes of
Clery-Mambrini-Olive-Verner arXiv:2112.15214 eqs.(83)-(85) and with their total beta_0 (eq.86).
Outputs: reaction (pair) rate and energy-transfer rate coefficients with Bose/Fermi statistics, and a
Monte-Carlo <T/E> of the produced hidden quanta (for the nonthermal induced mass, item 5).
"""
import numpy as np
from numpy.polynomial import legendre
from common import J3, J4, N_HIGGS_REAL, N_WEYL, N_VECTOR, dump

ETA = np.diag([1.0, -1.0, -1.0, -1.0])
s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [s0, sx, sy, sz]
Z2 = np.zeros((2, 2), dtype=complex)
GAM = [np.block([[Z2, s0], [s0, Z2]])] + [np.block([[Z2, s], [-s, Z2]]) for s in (sx, sy, sz)]
G5 = np.block([[-s0, Z2], [Z2, s0]])
PL = (np.eye(4) - G5) / 2


def low(p):
    return ETA @ p


def dot(a, b):
    return a @ ETA @ b


def contract(A, B):
    return np.einsum("mn,ma,nb,ab->", A, ETA, ETA, B)


def trace(A):
    return np.einsum("mn,mn->", ETA, A)


def T_scalar(p1, p2):
    a, b = low(p1), low(p2)
    return np.outer(a, b) + np.outer(b, a) - ETA * dot(p1, p2)


def T_final(k1, k2, xi):
    a, b = low(k1), low(k2)
    q = a + b
    Q2 = dot(k1 + k2, k1 + k2)
    return np.outer(a, b) + np.outer(b, a) - ETA * dot(k1, k2) - 2 * xi * (np.outer(q, q) - ETA * Q2)


def fstr(p, e):
    a, b = low(p), low(e)
    return np.outer(a, b) - np.outer(b, a)


def T_vector(p1, e1, p2, e2):
    f1, f2 = fstr(p1, e1), fstr(p2, e2)
    m1 = np.einsum("ma,ab,nb->mn", f1, ETA, f2)    # f1_{m a} f2_n^a
    m2 = np.einsum("ma,ab,nb->mn", f2, ETA, f1)
    ff = np.einsum("ab,ac,bd,cd->", f1, ETA, ETA, f2)
    return -(m1 + m2) + 0.5 * ETA * ff


def sqrtm_psd(M):
    w, U = np.linalg.eigh(M)
    w = np.clip(w, 0, None)
    return U @ np.diag(np.sqrt(w)) @ U.conj().T


def spinors(p):
    ps = p[0] * s0 - p[1] * sx - p[2] * sy - p[3] * sz
    psb = p[0] * s0 + p[1] * sx + p[2] * sy + p[3] * sz
    A, B = sqrtm_psd(ps), sqrtm_psd(psb)
    us, vs = [], []
    for xi in (np.array([1, 0], complex), np.array([0, 1], complex)):
        us.append(np.concatenate([A @ xi, B @ xi]))
        vs.append(np.concatenate([A @ xi, -B @ xi]))
    return us, vs


def slash(p):
    pl = low(p)
    return sum(GAM[m] * pl[m] for m in range(4))


def kinematics(cth, phi=0.3):
    E = 0.5  # s = 1
    p1 = np.array([E, 0, 0, E])
    p2 = np.array([E, 0, 0, -E])
    sth = np.sqrt(1 - cth**2)
    k1 = np.array([E, E * sth * np.cos(phi), E * sth * np.sin(phi), E * cth])
    k2 = np.array([E, -k1[1], -k1[2], -k1[3]])
    return p1, p2, k1, k2


def msq_scalar(cth, xi=0.0):
    p1, p2, k1, k2 = kinematics(cth)
    Ti, Tf = T_scalar(p1, p2), T_final(k1, k2, xi)
    s = dot(p1 + p2, p1 + p2)
    M = (contract(Ti, Tf) - 0.5 * trace(Ti) * trace(Tf)) / s
    return abs(M)**2


def msq_weyl(cth, xi=0.0, return_checks=False):
    p1, p2, k1, k2 = kinematics(cth)
    Tf = T_final(k1, k2, xi)
    s = dot(p1 + p2, p1 + p2)
    us, _ = spinors(p1)
    _, vs = spinors(p2)
    d = low(p1 - p2)
    tot = 0.0
    trmax = 0.0
    for u in us:
        for v in vs:
            vb = v.conj() @ GAM[0]
            J = np.array([vb @ GAM[m] @ PL @ u for m in range(4)])   # J^m (upper index)
            Jl = ETA @ J
            TF = 0.25 * (np.outer(Jl, d) + np.outer(d, Jl))
            trmax = max(trmax, abs(trace(TF)))
            M = (contract(TF, Tf) - 0.5 * trace(TF) * trace(Tf)) / s
            tot += abs(M)**2
    if return_checks:
        return tot, trmax
    return tot


def msq_vector(cth, xi=0.0, return_checks=False):
    p1, p2, k1, k2 = kinematics(cth)
    Tf = T_final(k1, k2, xi)
    s = dot(p1 + p2, p1 + p2)
    pols = [np.array([0, 1.0, 0, 0]), np.array([0, 0, 1.0, 0])]
    tot, trmax = 0.0, 0.0
    for e1 in pols:
        for e2 in pols:
            TV = T_vector(p1, e1, p2, e2)
            trmax = max(trmax, abs(trace(TV)))
            M = (contract(TV, Tf) - 0.5 * trace(TV) * trace(Tf)) / s
            tot += abs(M)**2
    if return_checks:
        return tot, trmax
    return tot


def clery_shapes(cth):
    s = 1.0
    t = -0.5 * s * (1 - cth)
    return {"00": t**2 * (s + t)**2 / (4 * s**2),
            "h0": -t * (s + t) * (s + 2 * t)**2 / (4 * s**2),
            "10": t**2 * (s + t)**2 / (2 * s**2)}


def kappa(fn, xi=0.0, n=40):
    """kappa = sigma_hat M_P^4 / s^2, sigma_hat = (1/2) int dPhi_2 sum|M|^2 = (1/32pi) int dc sum|M|^2."""
    x, w = legendre.leggauss(n)
    return sum(wi * fn(xi_, xi) for xi_, wi in zip(x, w)) / (32 * np.pi)


def newton_check():
    """Static limit: M = T1 P T2 /(M_P^2 q^2) with T^00 = 2m^2 gives V(r) = -G m1 m2 / r if P_0000 = 1/2."""
    P0000 = 0.5 * (1 + 1 - 1)
    m1 = m2 = 1.0
    MP = 1.0
    M_over_q2 = (2 * m1**2) * P0000 * (2 * m2**2) / MP**2
    Vq = -M_over_q2 / (2 * m1 * 2 * m2)           # V(q) * q^2
    G = 1.0 / (8 * np.pi * MP**2)
    return {"V(q) q^2": Vq, "-4 pi G m1 m2": -4 * np.pi * G * m1 * m2, "pass": bool(abs(Vq + 4 * np.pi * G) < 1e-14)}


def completeness_check():
    p = np.array([2.0, 0.6, -0.8, np.sqrt(4 - 1.0)])
    us, vs = spinors(p)
    G0 = GAM[0]
    su = sum(np.outer(u, u.conj() @ G0) for u in us)
    sv = sum(np.outer(v, v.conj() @ G0) for v in vs)
    return {"u_ubar_minus_pslash": float(np.abs(su - slash(p)).max()),
            "v_vbar_minus_pslash": float(np.abs(sv - slash(p)).max())}


# ------------------------------- Monte Carlo <T/E> ---------------------------------
def _sample_E(rng, kind, n):
    E = np.linspace(1e-6, 60, 200001)
    if kind == "B":
        pdf = E**3 / np.expm1(E)
    else:
        pdf = E**3 / (np.exp(E) + 1)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, E)


def _angular_poly(fn):
    x = np.linspace(-1, 1, 9)
    y = np.array([fn(c, 0.0) for c in x])
    return np.polyfit(x, y, 6)


def mc_inverse_energy(rng, kind, fn, n=400000):
    E1 = _sample_E(rng, kind, n)
    E2 = _sample_E(rng, kind, n)
    u = rng.random(n)
    c12 = 1 - 2 * (1 - u)**(1.0 / 3.0)
    s12 = np.sqrt(np.clip(1 - c12**2, 0, None))
    p1 = np.stack([E1, 0 * E1, 0 * E1, E1], 1)
    p2 = np.stack([E2, E2 * s12, 0 * E2, E2 * c12], 1)
    P = p1 + p2
    s = 2 * E1 * E2 * (1 - c12)
    beta = P[:, 1:] / P[:, :1]
    b2 = np.sum(beta**2, 1)
    gam = 1 / np.sqrt(1 - b2)

    def boost(vec, bvec, g, sign):
        # boost 4-vectors by velocity sign*bvec
        bv = sign * bvec
        bb = np.sum(bv**2, 1)
        bp = np.sum(bv * vec[:, 1:], 1)
        g2 = np.where(bb > 0, (g - 1) / np.where(bb > 0, bb, 1), 0)
        E = g * (vec[:, 0] + bp)
        sp = vec[:, 1:] + (g2 * bp)[:, None] * bv + (g * vec[:, 0])[:, None] * bv
        return np.concatenate([E[:, None], sp], 1)

    p1s = boost(p1, beta, gam, -1.0)
    n1 = p1s[:, 1:] / np.linalg.norm(p1s[:, 1:], axis=1)[:, None]
    coef = _angular_poly(fn)
    fmax = np.max(np.polyval(coef, np.linspace(-1, 1, 2001))) * 1.05
    cst = np.empty(n)
    todo = np.arange(n)
    while todo.size:
        c = 2 * rng.random(todo.size) - 1
        acc = rng.random(todo.size) * fmax < np.polyval(coef, c)
        cst[todo[acc]] = c[acc]
        todo = todo[~acc]
    phi = 2 * np.pi * rng.random(n)
    tmp = np.where(np.abs(n1[:, 0:1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    e1 = np.cross(n1, tmp)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(n1, e1)
    sn = np.sqrt(np.clip(1 - cst**2, 0, None))
    kd = cst[:, None] * n1 + (sn * np.cos(phi))[:, None] * e1 + (sn * np.sin(phi))[:, None] * e2
    Ek = np.sqrt(s) / 2
    k1s = np.concatenate([Ek[:, None], Ek[:, None] * kd], 1)
    k2s = np.concatenate([Ek[:, None], -Ek[:, None] * kd], 1)
    k1 = boost(k1s, beta, gam, +1.0)
    k2 = boost(k2s, beta, gam, +1.0)
    Eall = np.concatenate([k1[:, 0], k2[:, 0]])
    econs = np.max(np.abs(k1[:, 0] + k2[:, 0] - (E1 + E2)) / (E1 + E2))
    return {"mean_T_over_E": float(np.mean(1 / Eall)), "mean_E_over_T": float(np.mean(Eall)),
            "sem_T_over_E": float(np.std(1 / Eall) / np.sqrt(Eall.size)),
            "pair_energy_mean_over_T": float(np.mean(E1 + E2)), "energy_conservation_max_rel": float(econs)}


def run(seed=20261020):
    out = {"newton_normalisation_check": newton_check(), "dirac_completeness": completeness_check()}
    cs = np.linspace(-0.95, 0.95, 8)   # avoids c=0 where the Weyl shape (s+2t)^2 vanishes
    shapes = {"scalar": [], "weyl": [], "vector": []}
    trf, trv = 0.0, 0.0
    for c in cs:
        cl = clery_shapes(c)
        mw, tw = msq_weyl(c, return_checks=True)
        mv, tv = msq_vector(c, return_checks=True)
        trf, trv = max(trf, tw), max(trv, tv)
        shapes["scalar"].append(msq_scalar(c) / cl["00"])
        shapes["weyl"].append(mw / cl["h0"])
        shapes["vector"].append(mv / cl["10"])
    out["trace_TF_max"] = trf
    out["trace_TV_max"] = trv
    out["ratio_ours_to_clery_eqs83_85"] = {k: {"min": min(v), "max": max(v)} for k, v in shapes.items()}
    # analytic scalar check: sum|M|^2 = t^2 u^2 / s^2 (s=1)
    t = -0.5 * (1 - 0.37)
    u = -0.5 * (1 + 0.37)
    out["scalar_analytic_check_rel"] = abs(msq_scalar(0.37) - t**2 * u**2) / (t**2 * u**2)
    # kappa per initial species (sum over its spins), final: one real scalar with 1/2 symmetry factor
    kap = {"higgs_real": kappa(msq_scalar), "weyl_pair": kappa(msq_weyl), "vector": kappa(msq_vector)}
    out["kappa"] = kap
    out["kappa_analytic_scalar"] = 1.0 / (480 * np.pi)
    species = [("higgs_real", N_HIGGS_REAL, 0.5, "B"), ("weyl_pair", N_WEYL, 1.0, "F"), ("vector", N_VECTOR, 0.5, "B")]
    pair_rate = {}
    energy_rate = {}
    pair_rate_MB = {}
    for nm, N, S, st in species:
        pair_rate[nm] = N * S * kap[nm] * J3[st]**2 / (3 * np.pi**4)
        energy_rate[nm] = N * S * kap[nm] * 2 * J3[st] * J4[st] / (3 * np.pi**4)
        pair_rate_MB[nm] = N * S * kap[nm] * J3["MB"]**2 / (3 * np.pi**4)
    b_pair = sum(pair_rate.values())
    b_E1 = sum(energy_rate.values())
    out["per_hidden_real_dof"] = {"pair_rate_coeff_T8_over_MP4": b_pair, "pair_rate_by_species": pair_rate,
                                  "energy_rate_coeff_T9_over_MP4": b_E1, "energy_rate_by_species": energy_rate,
                                  "pair_rate_coeff_maxwell_boltzmann": sum(pair_rate_MB.values())}
    # central-charge check: for conformal sectors S_i kappa_i ∝ c_Weyl-anomaly, ratios real scalar : Weyl : vector = 1 : 3 : 12
    k_conf = kappa(msq_scalar, 1.0 / 6.0)
    out["central_charge_check"] = {"S*kappa conformal real scalar": 0.5 * k_conf,
                                   "Weyl/scalar": kap["weyl_pair"] / (0.5 * k_conf),
                                   "vector/scalar": 0.5 * kap["vector"] / (0.5 * k_conf),
                                   "expected (c = (N_s + 3 N_Weyl + 12 N_v)/120)": [3.0, 12.0]}
    beta0_clery = 3997 * np.pi**3 / 41472000
    out["ours_with_weyl_channel_doubled"] = b_pair + pair_rate["weyl_pair"]
    out["ours_with_weyl_doubled_minus_clery_beta0"] = b_pair + pair_rate["weyl_pair"] - beta0_clery
    out["clery_beta0_eq86"] = beta0_clery
    out["ratio_our_pair_rate_to_clery_beta0"] = b_pair / beta0_clery
    out["interpretation"] = ("ratio 1 means Clery et al. R_0^T equals the reaction (pair) rate per real scalar; "
                             "the particle-number rate is 2x that")
    # xi dependence of the final-state coupling (Higgs initial only; conformal initial states are xi-blind)
    xis = [-0.01, -0.0034, 0.0, 1.0 / 6.0, 1.0]
    out["xi_dependence_higgs_channel_kappa_ratio"] = {f"{x:.4g}": kappa(msq_scalar, x) / kap["higgs_real"] for x in xis}
    w_h = energy_rate["higgs_real"] / b_E1
    out["xi_dependence_total_energy_rate_ratio"] = {f"{x:.4g}": 1 - w_h + w_h * kappa(msq_scalar, x) / kap["higgs_real"]
                                                    for x in xis}
    NH = 4
    out["beta_E_4dof"] = NH * b_E1
    out["beta_pair_4dof"] = NH * b_pair
    out["beta_particle_number_4dof"] = 2 * NH * b_pair
    out["mean_pair_energy_over_T"] = b_E1 / b_pair
    # Monte Carlo <T/E>
    rng = np.random.default_rng(seed)
    mc = {}
    for nm, N, S, st in species:
        fn = {"higgs_real": msq_scalar, "weyl_pair": msq_weyl, "vector": msq_vector}[nm]
        mc[nm] = mc_inverse_energy(rng, st, fn)
    out["mc"] = mc
    wts = {nm: pair_rate[nm] / b_pair for nm, *_ in species}
    TE = sum(wts[nm] * mc[nm]["mean_T_over_E"] for nm in wts)
    out["mc_mean_T_over_E_weighted"] = TE
    out["mc_mean_E_over_T_weighted"] = sum(wts[nm] * mc[nm]["mean_E_over_T"] for nm in wts)
    out["analytic_mean_E_over_T_per_particle"] = 0.5 * b_E1 / b_pair
    # chi2 injection: d<chi^2>/dt per hidden dof = 2 * pair_rate * <1/E> = gamma_chi T^7 / M_P^4
    out["gamma_chi_per_dof_T7_over_MP4"] = 2 * b_pair * TE
    dump("c2_graviton_rates.json", out)
    return out


if __name__ == "__main__":
    r = run()
    for k in ["newton_normalisation_check", "dirac_completeness", "trace_TF_max", "trace_TV_max",
              "ratio_ours_to_clery_eqs83_85", "scalar_analytic_check_rel", "kappa", "kappa_analytic_scalar",
              "per_hidden_real_dof", "clery_beta0_eq86", "ratio_our_pair_rate_to_clery_beta0",
              "xi_dependence_higgs_channel_kappa_ratio", "xi_dependence_total_energy_rate_ratio", "beta_E_4dof",
              "beta_pair_4dof", "mean_pair_energy_over_T", "mc", "mc_mean_T_over_E_weighted",
              "mc_mean_E_over_T_weighted", "analytic_mean_E_over_T_per_particle", "gamma_chi_per_dof_T7_over_MP4"]:
        print(k, r[k])
