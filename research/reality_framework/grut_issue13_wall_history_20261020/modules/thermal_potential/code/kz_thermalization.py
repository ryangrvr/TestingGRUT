"""Items 4 and 5: Kibble-Zurek input and hidden self-thermalisation for a hidden bath at T_h = xi T_SM, visible sector
radiation dominated with g* = 106.75.  Uses Tc and kappa = d m_phi^2(0,T)/dT^2 from results/transition.json (DR scheme,
running couplings, mu0) and the DR potential for the m_phi(T) table.  Writes results/kibble_zurek.json and
results/thermalization.json."""
import json
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import zeta, spence

from common import BENCH, MU0, Veff, JBSeries, hubble, MP_RED, PI, V_GEV, to_builtin

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent.parent / "results"
HBARC_GEV_CM = 1.973269804e-14
NU_ISING = 1.0 / (3.0 - 1.412625)      # from Delta_epsilon (Kos, Poland, Simmons-Duffin, Vichi 2016)
G_STAR = 106.75
G_HID = 4.0                            # phi, two real X modes, D
XI_GRID = np.geomspace(1e-6, 1.0, 61)


def kz_naive(xi0, tau0, tauQ, nu, z):
    eps = (tau0 / tauQ) ** (1.0 / (1.0 + nu * z))
    return xi0 * eps ** (-nu), eps


def kz_matched(xi0, tau0, tauQ, nu, z, epsG):
    """Mean-field xi = xi0 eps^-1/2 for eps >= epsG, Ising xi = xi0 epsG^-1/2 (eps/epsG)^-nu below; tau = tau0 (xi/xi0)^z."""
    def corr(e):
        return xi0 * e ** -0.5 if e >= epsG else xi0 * epsG ** -0.5 * (e / epsG) ** (-nu)

    f = lambda le: np.log(tau0 * (corr(np.exp(le)) / xi0) ** z) - np.log(np.exp(le) * tauQ)
    le = brentq(f, np.log(1e-300), np.log(10.0), xtol=1e-12)
    e = np.exp(le)
    return corr(e), e


def main():
    t0 = time.time()
    tr = json.loads((OUT / "transition.json").read_text())
    dr0 = [r for r in tr["runs"] if r["scheme"] == "DR" and r["rg"] == "running" and r["scale_label"] == "mu0"][0]
    Tc, kappa = dr0["Tc_over_v"], dr0["kappa_dm2_dT2_at_Tc"]
    Tc_range = tr["adopted"]["Tc_over_v_range_resummed_running"]
    lp = BENCH["lp"]
    epsG = tr["adopted"]["ginzburg"]["t_G"]
    xi0 = 1.0 / (Tc * np.sqrt(2 * kappa))         # symmetric-phase mean-field amplitude
    tau0 = xi0                                    # convention: microscopic time = light crossing of xi0

    # m_phi(T) near Tc from the DR potential
    V = Veff(BENCH, MU0, "DR", JBSeries())
    mrows = []
    for e in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 1e-5, 1e-6]:
        Ts, Tb = Tc * (1 + e), Tc * (1 - e)
        m_sym = np.sqrt(2 * V.dW(0.0, Ts)[0])
        g = V.dW(np.geomspace(1e-14, 1, 20000), Tb)
        idx = np.where((g[:-1] < 0) & (g[1:] >= 0))[0]
        s_grid = np.geomspace(1e-14, 1, 20000)
        sb = brentq(lambda s: V.dW(s, Tb)[0], s_grid[idx[0]], s_grid[idx[0] + 1], xtol=1e-20, rtol=1e-14)
        w2 = (V.dW(sb * (1 + 1e-5), Tb)[0] - V.dW(sb * (1 - 1e-5), Tb)[0]) / (2e-5 * sb)
        m_brk = np.sqrt(4 * sb * w2)
        mrows.append({"eps": e, "m_sym_over_Tc": m_sym / Tc, "mf_sqrt(2 kappa eps)": np.sqrt(2 * kappa * e),
                      "m_broken_over_Tc": m_brk / Tc, "mf_sqrt(4 kappa eps)": np.sqrt(4 * kappa * e),
                      "phi_broken_over_Tc": np.sqrt(sb) / Tc,
                      "eps_phi_sym_3lpT_over_4pi_m": 3 * lp * Tc * (1 + e) / (4 * PI * m_sym),
                      "corr_length_sym_GeV^-1": 1.0 / (m_sym * V_GEV)})

    rows = []
    for xi in XI_GRID:
        Tsm = Tc / xi
        H = hubble(Tsm, G_STAR)
        tauQ = 1.0 / H
        hub = 1.0 / H
        row = {"xi": xi, "T_SM_GeV": Tsm * V_GEV, "H_GeV": H * V_GEV, "tau_Q_GeV^-1": tauQ / V_GEV,
               "H_correction_factor_with_hidden_energy": np.sqrt(1 + G_HID * xi**4 / G_STAR),
               "Hubble_length_GeV^-1": hub / V_GEV, "Hubble_length_cm": hub / V_GEV * HBARC_GEV_CM}
        for z in (1, 2):
            lK, eK = kz_naive(xi0, tau0, tauQ, NU_ISING, z)
            lM, eM = kz_matched(xi0, tau0, tauQ, NU_ISING, z, epsG)
            lF, eF = kz_naive(xi0, tau0, tauQ, 0.5, z)
            for tag, (l_, e_) in {"ising_naive": (lK, eK), "ising_ginzburg_matched": (lM, eM), "meanfield_nu_half": (lF, eF)}.items():
                l_c = min(l_, hub)
                row["z%d_%s" % (z, tag)] = {"eps_freeze": e_, "xi_KZ_over_Tc^-1": l_ * Tc, "xi_KZ_GeV^-1": l_ / V_GEV,
                                            "causal_cap_active": bool(l_ > hub), "N_H": (hub / l_c) ** 3,
                                            "log10_N_H": 3 * np.log10(hub / l_c)}
            # sensitivity to tau0 by factors 1e-2, 1e2 (Ising, matched)
            row["z%d_log10_N_H_tau0_x1e-2_x1e2" % z] = [
                3 * np.log10(hub / min(kz_matched(xi0, tau0 * f, tauQ, NU_ISING, z, epsG)[0], hub)) for f in (1e-2, 1e2)]
        lG = xi0 * epsG ** -0.5
        row["kibble_ginzburg"] = {"xi_G_GeV^-1": lG / V_GEV, "log10_N_H": 3 * np.log10(hub / lG)}
        rows.append(row)
    decades = [r for r in rows if abs(np.log10(r["xi"]) - round(np.log10(r["xi"]))) < 1e-9]
    kz = {"inputs": {"Tc_over_v_DR_mu0": Tc, "Tc_range_over_v": Tc_range, "kappa": kappa, "xi0_over_Tc^-1": xi0 * Tc,
                     "tau0_convention": "tau0 = xi0 (c = 1)", "nu_ising": NU_ISING, "nu_source": "1/(3 - Delta_eps), Delta_eps = 1.412625 (arXiv:1603.04436 abstract)",
                     "z_values": [1, 2], "eps_G": epsG, "g_star": G_STAR, "M_P_GeV": MP_RED * V_GEV,
                     "assumptions": ["T_h proportional to 1/a in both sectors (separately conserved entropy, constant g's) so tau_Q = 1/H",
                                     "visible sector radiation dominated at T_SM = Tc/xi (requires reheating above Tc/xi: preparation data outside the action)",
                                     "z = 1: underdamped relativistic order-parameter dynamics; z = 2: relaxational (model A-like) dynamics",
                                     "tau0 = xi0 and xi0 = mean-field symmetric-phase amplitude; Ising exponent applied naive or matched at the Ginzburg point",
                                     "the bath is internally thermal (see thermalization.json); otherwise KZ is not applicable"]},
          "m_phi_near_Tc_DR": mrows, "rows": rows, "decades": decades}
    lo = {k: min(r[k]["log10_N_H"] for r in rows) for k in ["z1_ising_ginzburg_matched", "z2_ising_ginzburg_matched", "z1_ising_naive", "z2_ising_naive"]}
    hi = {k: max(r[k]["log10_N_H"] for r in rows) for k in lo}
    kz["log10_N_H_range_over_xi_grid"] = {k: [lo[k], hi[k]] for k in lo}
    kz["causal_bound_ever_active"] = any(r[k]["causal_cap_active"] for r in rows for k in r if k.startswith("z") and isinstance(r[k], dict))
    (OUT / "kibble_zurek.json").write_text(json.dumps(to_builtin(kz), indent=1))

    # ---------------------------------------------------------------- thermalization
    n1 = zeta(3) / PI**2                      # n per real bosonic dof / T^3
    Eav = PI**4 / (30 * zeta(3))               # <E>/T
    s_av = 2 * Eav**2                          # <s>/T^2 (isotropic average of 2 E1 E2 (1 - cos))
    pref = n1 / (16 * PI * s_av)               # Gamma/T per unit sum |M|^2 S_f
    p = BENCH
    lpp, lx, ld, c, d = p["lp"], p["lx"], p["ld"], p["c"], p["d"]
    sums = {
        "phi_elastic_total": 18 * lpp**2 + 12 * c**2 + 6 * d**2,
        "X_elastic_total": 24 * lx**2 + 6 * c**2,
        "D_exchange_with_phi": 6 * d**2,
        "D_self": 18 * ld**2,
    }
    rates = {k: pref * v_ for k, v_ in sums.items()}            # Gamma / T_h
    # number changing 2 -> 4 (phi phi -> 4 phi) estimate
    Nd, Sf = 10.0, 24.0
    Phi4 = 1.0 / (24576 * PI**5)                                # Phi_4 / s^2
    def g24(kp, nch):
        return n1 * nch * Nd**2 * (6 * lpp) ** 4 * kp**2 * Phi4 / (2 * Sf * s_av)    # Gamma/T
    rates["number_changing_2to4_central"] = g24(4.0, 1.0)
    rates["number_changing_2to4_low"] = g24(1.0, 1.0)
    rates["number_changing_2to4_high"] = g24(16.0, 4.0)
    H_over_xi2 = hubble(Tc, G_STAR)                             # H at T_SM = Tc/xi is H_over_xi2 / xi^2
    xi_min = {k: np.sqrt(H_over_xi2 / (g * Tc)) for k, g in rates.items()}
    th = {"formulas": {
              "rate": "Gamma_a = sum_b n_b sum_f |M_ab->f|^2 S_f / (16 pi <s>), n_b = zeta(3) T^3/pi^2, <s> = 2 <E>^2, <E> = pi^4 T/(30 zeta(3))",
              "amplitudes": "phi phi->phi phi: 6 lambda_phi; phi phi->q_i q_i, phi q_i->phi q_i: 2 g_X^2; q q->q q: 6 lambda_X (same), 2 lambda_X (q1q1->q2q2, q1q2->q1q2); "
                            "phi phi->DD, phi D->phi D: 2 h_D^2; DD->DD: 6 lambda_D; alpha = 0 at mu0",
              "sum_|M|^2_S_f": sums,
              "2to4": "sigma = |M4|^2 Phi_4/(2 s S_f), Phi_4 = s^2/(24576 pi^5), |M4| = N_d (6 lambda_phi)^2 kappa_p/s, N_d = 10, S_f = 4!; "
                      "central kappa_p = 4, N_ch = 1; band kappa_p in [1,16], N_ch in [1,4]",
              "criterion": "Gamma(T_h = Tc) > H(T_SM = Tc/xi); Gamma/H scales as xi^2 / T_h for contact interactions, so a bath thermal at Tc stays thermal below it",
              "estimate_quality": "order of magnitude; collinear logs, Bose enhancement and thermal masses not resolved (factor ~3-10 on 2->2 rates)"},
          "Gamma_over_T": rates, "Gamma_over_H_at_xi=1": {k: g * Tc / H_over_xi2 for k, g in rates.items()},
          "xi_min_for_Gamma>H": xi_min,
          "number_chemical_potential_sensitivity": {
              "statement": "kinetic equilibrium alone fixes T_h but allows a number chemical potential mu_N <= 0 if 2->4 is slow; "
                           "the leading Debye coefficient scales as Li_2(z)/zeta(2), z = exp(mu_N/T_h), so Tc (in kinetic temperature) "
                           "scales as sqrt(zeta(2)/Li_2(z))",
              "example_z": {str(zz): float(np.sqrt((PI**2 / 6) / spence(1 - zz))) for zz in (1.0, 0.9, 0.5, 0.1)}},
          "grid": [{"xi": xi, "phi_elastic_Gamma_over_H": rates["phi_elastic_total"] * Tc * xi**2 / H_over_xi2,
                    "X_elastic_Gamma_over_H": rates["X_elastic_total"] * Tc * xi**2 / H_over_xi2,
                    "D_exchange_Gamma_over_H": rates["D_exchange_with_phi"] * Tc * xi**2 / H_over_xi2,
                    "two_to_four_central_Gamma_over_H": rates["number_changing_2to4_central"] * Tc * xi**2 / H_over_xi2}
                   for xi in XI_GRID]}
    (OUT / "thermalization.json").write_text(json.dumps(to_builtin(th), indent=1))
    print("xi_min:", json.dumps(to_builtin(xi_min), indent=1))
    print("log10 N_H ranges:", json.dumps(to_builtin(kz["log10_N_H_range_over_xi_grid"])), "causal cap ever:", kz["causal_bound_ever_active"])
    for r in decades:
        print("xi=%.0e T_SM=%.3e GeV H=%.3e GeV tauQ*Tc=%.3e  log10N_H z1 matched %.2f naive %.2f | z2 matched %.2f naive %.2f | KG %.2f; eps1 %.2e" % (
            r["xi"], r["T_SM_GeV"], r["H_GeV"], r["tau_Q_GeV^-1"] * Tc * V_GEV,
            r["z1_ising_ginzburg_matched"]["log10_N_H"], r["z1_ising_naive"]["log10_N_H"],
            r["z2_ising_ginzburg_matched"]["log10_N_H"], r["z2_ising_naive"]["log10_N_H"], r["kibble_ginzburg"]["log10_N_H"],
            r["z1_ising_ginzburg_matched"]["eps_freeze"]))
    for m in mrows:
        print(m)
    print("time %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
