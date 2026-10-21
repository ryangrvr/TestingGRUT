"""cosmo_v.py -- verifier's independent recomputation of the coarsening factors, T_par, T_safe and thresholds.

Independent of modules/bias_threshold/code (not imported, not copied):
  * own g_rho(T), g_s(T): ideal-gas integrals for every SM species (FD/BE, masses), a QCD crossover modelled as a
    tanh switch between a quark-gluon gas and a pion gas at T_QCD = 0.155 GeV (width 0.02 GeV), neutrino decoupling
    and e+e- annihilation by entropy conservation (T_nu = (4/11)^(1/3) T after annihilation; N_eff = 3 or 3.044).
    Precision: a few % in g around the QCD and e+e- epochs, which enters K as g^(1/3) and T_safe as g^(-1/2).
  * own radiation-era formulas: rho_w = A sigma / t, t = 1/(2H), parity rho_w = rho_rad, relic-safe
    f kappa rho_w/s <= Y_max with Y_max from Omega h^2 = 0.12 (rho_c/h^2 = 1.05368e-5 GeV cm^-3, s0 = 2891.2 cm^-3).
  * thresholds eps(K) = erf(x10 K^-gamma), K10(m) = (x10/erfinv m)^(1/gamma).
Writes results/cosmo_v.json.
"""
import json
import math
import os

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erf, erfinv

HERE = os.path.dirname(os.path.abspath(__file__))
MODS = os.path.join(HERE, "..", "..", "modules")

V = 1000.0
LAM = 0.01
SIGMA = 2 * math.sqrt(2 * LAM) / 3 * V ** 3
M_P = 2.435323e18
RHO_C_H2 = 1.05368e-5      # GeV cm^-3  (PDG: rho_c = 1.05368e-5 h^2 GeV cm^-3)
S0 = 2891.2                # cm^-3      (T0 = 2.7255 K)
YMAX = 0.12 * RHO_C_H2 / S0

# ---------------------------------------------------------------- ideal-gas degrees of freedom
# (name, mass GeV, internal dof, fermion?)
LEPTONS = [("e", 0.000511, 4, True), ("mu", 0.10566, 4, True), ("tau", 1.777, 4, True)]
QUARKS = [("u", 0.0022, 12, True), ("d", 0.0047, 12, True), ("s", 0.095, 12, True), ("c", 1.27, 12, True),
          ("b", 4.18, 12, True), ("t", 172.7, 12, True)]
BOSONS_EW = [("W", 80.4, 6, False), ("Z", 91.19, 3, False), ("h", 125.1, 1, False)]
GLUON = [("g", 0.0, 16, False)]
PIONS = [("pi0", 0.135, 1, False), ("pic", 0.1396, 2, False)]


def _rho_s(m, T, fermion):
    """rho/T^4 and s/T^3 per internal dof (zero chemical potential)."""
    x = m / T
    sgn = 1.0 if fermion else -1.0
    if x > 60:
        return 0.0, 0.0
    def occ(u):
        e = math.sqrt(u * u + x * x)
        return 1.0 / (math.exp(e) + sgn) if e < 700 else 0.0
    rho = quad(lambda u: u * u * math.sqrt(u * u + x * x) * occ(u), 0, 80 + x, limit=200)[0] / (2 * math.pi ** 2)
    p = quad(lambda u: u ** 4 / math.sqrt(u * u + x * x) * occ(u) if (u > 0 or x > 0) else 0.0, 0, 80 + x,
             limit=200)[0] / (6 * math.pi ** 2)
    return rho, rho + p


def _g_of(species, T):
    gr = gs = 0.0
    for _, m, g, ferm in species:
        r, s = _rho_s(m, T, ferm)
        gr += g * r * 30 / math.pi ** 2
        gs += g * s * 45 / (2 * math.pi ** 2)
    return gr, gs


def g_both(T, neff=3.0):
    """Returns (g_rho, g_s) of the visible plasma at photon temperature T."""
    TQ, W = 0.155, 0.02
    wq = 0.5 * (1 + math.tanh((T - TQ) / W))
    # photons
    gr, gs = 2.0, 2.0
    # charged leptons + EW bosons (in equilibrium with photons)
    a, b = _g_of(LEPTONS[1:] + BOSONS_EW, T)
    gr += a
    gs += b
    # QCD sector: tanh switch between quark-gluon gas and pion gas
    if wq > 1e-6:
        a, b = _g_of(QUARKS + GLUON, T)
        gr += wq * a
        gs += wq * b
    if wq < 1 - 1e-6:
        a, b = _g_of(PIONS, T)
        gr += (1 - wq) * a
        gs += (1 - wq) * b
    # electrons and neutrinos: coupled above T_dec = 2 MeV; below, neutrinos decoupled with T_nu from
    # entropy conservation of the (gamma, e+-) plasma.
    re_, se_ = _rho_s(0.000511, T, True)
    ge_r, ge_s = 4 * re_ * 30 / math.pi ** 2, 4 * se_ * 45 / (2 * math.pi ** 2)
    if T > 0.002:
        gr += ge_r + 7 / 8 * 2 * 3
        gs += ge_s + 7 / 8 * 2 * 3
    else:
        # (2 + g_e(T)) (aT)^3 = const = (2 + 3.5) (a T_dec)^3 ;  T_nu = T_dec a_dec/a
        tnu_over_t = ((2 + ge_s) / (2 + 3.5)) ** (1 / 3)
        gr += ge_r + 7 / 8 * 2 * neff * tnu_over_t ** 4
        gs += ge_s + 7 / 8 * 2 * neff * tnu_over_t ** 3
    return gr, gs


_cache = {}


def g(T, neff=3.0):
    key = (round(math.log(T), 6), neff)
    if key not in _cache:
        _cache[key] = g_both(T, neff)
    return _cache[key]


def H(T, neff=3.0):
    return math.sqrt(math.pi ** 2 / 90 * g(T, neff)[0]) * T * T / M_P


def t_of(T, neff=3.0):
    return 1 / (2 * H(T, neff))


def rho_rad(T, neff=3.0):
    return math.pi ** 2 / 30 * g(T, neff)[0] * T ** 4


def s_of(T, neff=3.0):
    return 2 * math.pi ** 2 / 45 * g(T, neff)[1] * T ** 3


def T_par(A, neff=3.0):
    f = lambda lt: math.log(A * SIGMA / t_of(math.exp(lt), neff)) - math.log(rho_rad(math.exp(lt), neff))
    return math.exp(brentq(f, math.log(1e-7), math.log(1e-3)))


def T_safe(A, fk, Ymax=YMAX, neff=3.0):
    f = lambda lt: math.log(fk * A * SIGMA / t_of(math.exp(lt), neff) / s_of(math.exp(lt), neff)) - math.log(Ymax)
    return math.exp(brentq(f, math.log(1e-6), math.log(10.0)))


def K_of(T, T_f, L_f, A, sQ=1.0, neff=3.0):
    L_T = t_of(T, neff) / (sQ * A)
    return L_T * (g(T, neff)[1] / g(T_f, neff)[1]) ** (1 / 3) * (T / T_f) / L_f


def T_of_H(Hv):
    return math.exp(brentq(lambda lt: math.log(H(math.exp(lt))) - math.log(Hv), math.log(1e-6), math.log(1e18)))


def eps(K, gamma, x10):
    return float(erf(x10 * K ** (-gamma)))


def K10(m, gamma, x10):
    return (x10 / float(erfinv(m))) ** (1 / gamma)


def main():
    out = {"Ymax_GeV": YMAX}
    # g checks
    out["g_checks"] = {f"{T:g}": g(T) for T in (3.11e6, 1e3, 1.0, 0.2, 0.1, 0.0306, 0.003, 1e-5, 6e-6)}
    # module's 3D calibration points (y_f, K10) for x10(gamma); my own refit is in fit_ns_v.json
    fk = json.load(open(os.path.join(MODS, "bias_threshold", "results", "fit_K10.json")))
    pts = fk["prs3d_bias_N128_L128"]["formation_fit"]["points"]
    x10_for = lambda gm: float(np.exp(np.median([math.log(p["y_f"]) + gm * math.log(p["K10"]) for p in pts])))
    # KZ inputs (H-star)
    kz = json.load(open(os.path.join(MODS, "thermal_potential", "results", "kibble_zurek.json")))
    row = min(kz["rows"], key=lambda r: abs(math.log(r["xi"] / 5e-4)))
    T_f = row["T_SM_GeV"]
    xi_main = [row["z1_ising_ginzburg_matched"]["xi_KZ_GeV^-1"], row["z2_ising_ginzburg_matched"]["xi_KZ_GeV^-1"]]
    xi_c = math.sqrt(xi_main[0] * xi_main[1])
    out["H_star_inputs"] = {"T_f": T_f, "H_f_kz_json": row["H_GeV"], "H_f_mine": H(T_f), "xi_main": xi_main, "xi_c": xi_c}
    res = {}
    for A in (0.7, 0.8, 0.9):
        res[f"T_par_keV_A{A}"] = T_par(A) * 1e6
        res[f"T_par_keV_A{A}_Neff3.044"] = T_par(A, 3.044) * 1e6
    for fkap in (1.0, 2.0, 0.1, 0.01, 0.02):
        res[f"T_safe_MeV_fk{fkap}"] = T_safe(0.8, fkap) * 1e3
    res["T_safe_MeV_band_A0.7-0.9_fk0.01-2"] = [T_safe(0.7, 0.01) * 1e3, T_safe(0.9, 2.0) * 1e3]
    # closed form check of T_safe at f kappa = 1 using my g at that T
    Ts = T_safe(0.8, 1.0)
    gr, gs = g(Ts)
    res["T_safe_closed_form_MeV"] = 0.8 * SIGMA * (45 / math.pi ** 2) * math.sqrt(math.pi ** 2 * gr / 90) / (gs * M_P * YMAX) * 1e3
    res["g_at_T_safe"] = [gr, gs]
    Tp = T_par(0.8)
    Kp = K_of(Tp, T_f, xi_c, 0.8)
    Ks = K_of(Ts, T_f, xi_c, 0.8)
    Kh = K_of(T_f, T_f, xi_c, 0.8)
    res["log10_K_par"] = math.log10(Kp)
    res["log10_K_par_backreaction_6.37keV"] = math.log10(K_of(6.37e-6, T_f, xi_c, 0.8))
    res["log10_K_safe"] = math.log10(Ks)
    res["log10_K_hor"] = math.log10(Kh)
    # band over A, xi0 (z1/z2 GM), g_geo 1-3, f kappa 0.01-2
    lkp, lks = [], []
    for A in (0.7, 0.8, 0.9):
        for xi0 in xi_main:
            for gg in (1.0, 3.0):
                for Tpp in (T_par(A), 6.37e-6):
                    lkp.append(math.log10(K_of(Tpp, T_f, gg * xi0, A)))
                for fkap in (0.01, 0.02, 1.0, 2.0):
                    lks.append(math.log10(K_of(T_safe(A, fkap), T_f, gg * xi0, A)))
    res["log10_K_par_band"] = [min(lkp), max(lkp)]
    res["log10_K_safe_band"] = [min(lks), max(lks)]
    # thresholds at central inputs vs gamma
    thr = {}
    for gm in (1.15, 1.19, 1.30, 1.33, 1.35, 1.39, 1.406, 1.46, 1.50, 0.84):
        x10 = x10_for(gm) if gm > 1.0 else 1.55
        thr[f"{gm}"] = {"x10": x10, "eps_persist": eps(Kp, gm, x10), "eps_safe": eps(Ks, gm, x10),
                        "eps_persist_K0=2.3K10": eps(Kp / 2.3, gm, x10)}
    res["thresholds_vs_gamma"] = thr
    # extreme band (my own): gamma 1.15-1.50, all nuisance inputs above
    ep, es = [], []
    for gm in (1.15, 1.50):
        x10 = x10_for(gm)
        for lk in lkp:
            for k0 in (1.0, 2.3):
                ep.append(math.log10(eps(10 ** lk / k0, gm, x10)))
        for lk in lks:
            es.append(math.log10(eps(10 ** lk, gm, x10)))
    res["log10_eps_persist_band_gamma1.15-1.50"] = [min(ep), max(ep)]
    res["log10_eps_safe_band_gamma1.15-1.50"] = [min(es), max(es)]
    # candidates
    gm, x10 = 1.406, x10_for(1.406)
    cands = {"thermal 1e-3": math.erf(1e-3 / math.sqrt(2)), "trap 0.4999": 2e-4, "trap 0.4997": 6e-4,
             "forced 0.71": 0.42, "relaxed 0.23": 0.54}
    res["candidates_log10_K10"] = {k: {"m": v, "gamma1.406": math.log10(K10(v, gm, x10)),
                                       "gamma1.15": math.log10(K10(v, 1.15, x10_for(1.15))),
                                       "gamma0.84": math.log10(K10(v, 0.84, 1.55))} for k, v in cands.items()}
    m_rem = cands["thermal 1e-3"]
    res["damping_efolds"] = {"to_eps_persist": math.log(m_rem / thr["1.406"]["eps_persist"]),
                             "to_eps_safe": math.log(m_rem / thr["1.406"]["eps_safe"]),
                             "to_eps_persist_gamma1.19": math.log(m_rem / thr["1.19"]["eps_persist"]),
                             "to_eps_safe_gamma1.19": math.log(m_rem / thr["1.19"]["eps_safe"])}
    res["statistical_floor_K_par^-1.5"] = Kp ** -1.5
    res["log10_eps_persist_minus_log10_floor"] = math.log10(thr["1.406"]["eps_persist"]) - math.log10(Kp ** -1.5)
    out["H_star"] = res
    # H-star-Q (T_f = 5.8e5 GeV, xi0 = sqrt(21*960) ell)
    TfQ = 0.29 * V / 5e-4
    xiQ = math.sqrt(21 * 960) / 100.0
    out["H_star_Q"] = {"T_f": TfQ, "xi0": xiQ, "log10_K_par": math.log10(K_of(Tp, TfQ, xiQ, 0.8)),
                       "log10_K_safe": math.log10(K_of(Ts, TfQ, xiQ, 0.8)),
                       "log10_K_hor": math.log10(K_of(TfQ, TfQ, xiQ, 0.8))}
    out["H_star_Q"]["eps_persist"] = eps(10 ** out["H_star_Q"]["log10_K_par"], 1.406, x10_for(1.406))
    out["H_star_Q"]["eps_safe"] = eps(10 ** out["H_star_Q"]["log10_K_safe"], 1.406, x10_for(1.406))
    # inflationary H_inf = 1e8 (central xi0 = 1/m_phi)
    infl = {}
    for Hinf in (2.4e4, 1e8, 1e13):
        Htr = math.sqrt(LAM) * V ** 2 / (1.5 * Hinf)
        Ttr = T_of_H(Htr)
        xi0 = 1 / (math.sqrt(2 * LAM) * V)
        Kp_i = K_of(Tp, Ttr, xi0, 0.8)
        infl[f"{Hinf:g}"] = {"H_trap": Htr, "T_trap": Ttr, "log10_K_par": math.log10(Kp_i),
                             "eps_persist": eps(Kp_i, 1.406, x10_for(1.406)),
                             "eps_safe": eps(K_of(Ts, Ttr, xi0, 0.8), 1.406, x10_for(1.406)),
                             "xi0_hubble_patch": 1 / math.sqrt(Hinf * Htr)}
    out["inflationary"] = infl
    out["x10_for"] = {f"{gm}": x10_for(gm) for gm in (1.15, 1.19, 1.30, 1.406, 1.50)}
    json.dump(out, open(os.path.join(HERE, "results", "cosmo_v.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    o = main()
    print(json.dumps(o, indent=1))
