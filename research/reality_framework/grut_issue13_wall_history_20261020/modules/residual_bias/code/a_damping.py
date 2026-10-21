"""Route A: thermal restoration of a previously selected uniform +v vacuum.  Damping of the soft phi zero mode in
the restored phase, the integrated damping over the restored epoch, and the residual systematic bias m0 at the exit.

1. On-shell width of a phi quantum at rest (k = 0) in a mu_N = 0 hidden bath (leading order, 2 <-> 2):
       Gamma = sum_ch mult |M|^2 S /(64 pi^3 m_a) Int dE_c dE_d F_stat Theta_kin,
       F_stat = f_b (1+f_c)(1+f_d) - (1+f_b) f_c f_d = f_b (1+f_c)(1+f_d)(1 - exp(-m_a/T)),
   with E_b = E_c + E_d - m_a, p_b = p_c + p_d and Theta_kin = [|K| <= p_c p_d],
   K = E_c E_d - m_a (E_c + E_d) + (m_a^2 - m_b^2 + m_c^2 + m_d^2)/2   (derivation in NOTES section A.2).
   Variants: 'full' (Bose), 'classical' (f = T/E; reproduces Parwani's lambda^2 T^2/(1536 pi m) exactly for one
   field), 'nsv' (Boltzmann collision rate of the k = 0 quantum: f_b only, no final-state factors, no inverse).
   1 <-> 2 processes: absent at linear order (no cubic vertex at phi = 0, exact Z2); with a coherent background
   phi_c the induced cubic 2 g_X^2 phi_c phi q q gives phi -> q q only if m_phi > 2 m_q, checked below (closed).
2. Amplitude damping exponent D = Int (Gamma/2) dt over the restored epoch (overdamped root m^2/Gamma used where
   Gamma > 2m), cut at the Ginzburg reduced temperature t_G (no damping credited inside it).
3. Coherent amplitude at the decision time (KZ freeze-out, z = 1 Ginzburg-matched xi_KZ from thermal_potential),
   linear WKB a^-3/2 m^-1/2 from production (phi_0 = v, sudden production on the cold +v vacuum [P]) times exp(-D);
   Hubble-only case with the nonlinear (quartic) cap; soft thermal rms on the correlation scale
   sigma^2 = c_w T_c m_hat, c_w = (sqrt(pi)/2 - (pi/2) e erfc(1))/(2 pi^2) = 0.01087 (Gaussian window radius 1/m_hat).
   m0 = erf(r/sqrt2), r = phi_coh/sigma.
4. Dilute (non-thermal, free-streaming) restoring gas with equal sharing: Gamma_d = C J/m_phi,
   C = (18 l^2 + 12 g^4 + 6 h^4)/(128 pi), Pi = c_eff J  =>  D = (C/(2 c_eff)) pi A/(2 sqrt(l) v), A = Pi/H.
5. Comparisons: production_map WKB memory estimate 1.7 v sqrt(l) v / A vs 0.11 sqrt(l) v [Q].
"""
import json
import math
import os
import sys
import time
import warnings

import numpy as np
from scipy.integrate import quad
from scipy.special import erfc

from common import (CHANNELS, EPS_G, GSTAR, LPHI, GX, HD, MP, PI_PHI, PI_Q, TC, TC_RANGE, V, dump, hubble,
                    log10_m0, masses)

warnings.filterwarnings("ignore", category=RuntimeWarning)
TP_RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "thermal_potential", "results")


# ------------------------------------------------------------------------------------------------ kinematics
def intervals(Ec, ma, mb, mc, md):
    """Allowed E_d intervals for fixed E_c: K(E_d)^2 <= p_c^2 p_d^2 (quadratic in E_d), E_d >= m_d, E_b >= m_b."""
    pc2 = max(Ec * Ec - mc * mc, 0.0)
    c0 = 0.5 * (ma * ma - mb * mb + mc * mc + md * md)
    A = ma * ma + mc * mc - 2 * ma * Ec
    B = 2 * (Ec - ma) * (c0 - ma * Ec)
    C = (c0 - ma * Ec) ** 2 + pc2 * md * md
    lo = max(md, mb + ma - Ec)
    out = []
    if abs(A) < 1e-300:
        return out
    disc = B * B - 4 * A * C
    if disc < 0:
        if A < 0:
            out.append((lo, np.inf))
        return out
    s = math.sqrt(disc)
    r0, r1 = sorted([(-B - s) / (2 * A), (-B + s) / (2 * A)])
    if A > 0:
        a_, b_ = max(lo, r0), r1
        if b_ > a_:
            out.append((a_, b_))
    else:
        if r0 > lo:
            out.append((lo, r0))
        out.append((max(lo, r1), np.inf))
    return out


def bose(E):
    return 1.0 / np.expm1(E)


def channel_integral(ma, mb, mc, md, mode, rtol=1e-7):
    """Int dE_c dE_d F Theta in units T = 1 (returns the dimensionless double integral)."""
    def F(Ec, Ed):
        Eb = Ec + Ed - ma
        if mode == "full":
            return bose(Eb) * (1 + bose(Ec)) * (1 + bose(Ed)) * (-math.expm1(-ma))
        if mode == "classical":
            return ma / (Eb * Ec * Ed)
        if mode == "nsv":
            return bose(Eb)
        raise ValueError(mode)

    def inner(Ec):
        tot = 0.0
        for a_, b_ in intervals(Ec, ma, mb, mc, md):
            tot += quad(lambda Ed: F(Ec, Ed), a_, b_, limit=200, epsabs=0, epsrel=rtol * 0.1)[0]
        return tot

    return quad(inner, mc, np.inf, limit=400, epsabs=0, epsrel=rtol)[0]


def gamma_over_T(T_h, mode="full", per_channel=False, Tc=TC):
    """On-shell zero-momentum width Gamma/T_h (energy/occupation relaxation rate; amplitude decays at Gamma/2)."""
    m = masses(T_h, Tc)
    x = {k: val / T_h for k, val in m.items()}
    rows = {}
    tot = 0.0
    for name, b, c, d, M2, S, mult in CHANNELS:
        I = channel_integral(x["phi"], x[b], x[c], x[d], mode)
        g = mult * M2 * S / (64 * math.pi ** 3 * x["phi"]) * I
        rows[name] = g
        tot += g
    return (tot, rows, x) if per_channel else tot


# ------------------------------------------------------------------------------------------------ checks
def parwani_check():
    """One field, all masses m: classical integral must equal pi^2/(6 m) (Parwani 1992 eq. 3.29:
    gamma = g^2 m/(64 pi) with m^2 = g^2 T^2/24, i.e. gamma = g^4 T^2/(1536 pi m); Gamma = 2 gamma)."""
    vals = {}
    for m in (0.01, 0.1, 1.0):
        I = channel_integral(m, m, m, m, "classical") / m      # = Int 1/(EbEcEd)
        vals[str(m)] = dict(I_times_m=I * m, pi2_over_6=math.pi ** 2 / 6, rel_dev=I * m / (math.pi ** 2 / 6) - 1)
    # full Bose for the pure lambda phi^4 theory at its own leading thermal mass: g^2 = 6 lambda_phi
    g2 = 6 * LPHI
    mT = math.sqrt(g2 / 24.0)
    Ifull = channel_integral(mT, mT, mT, mT, "full")
    G_full = g2 ** 2 * 0.5 / (64 * math.pi ** 3 * mT) * Ifull
    G_parwani = 2 * g2 ** 2 / (1536 * math.pi * mT)
    vals["pure_phi4_g2=6lambda_phi"] = dict(m_over_T=mT, Gamma_full_over_T=G_full, Gamma_Parwani_over_T=G_parwani,
                                           ratio_full_over_LO=G_full / G_parwani)
    return vals


def one_to_two_check():
    """phi -> q q with the background-induced cubic vertex needs m_phi > 2 m_q; phi -> D D needs m_phi > 2 m_D."""
    rows = []
    for t in (1.0005, 1.01, 1.1, 2.0, 10.0, 1e3):
        m = masses(t * TC)
        rows.append(dict(T_over_Tc=t, m_phi=m["phi"], two_m_q=2 * m["q"], two_m_D=2 * m["D"],
                         phi_to_qq_open=m["phi"] > 2 * m["q"], phi_to_DD_open=m["phi"] > 2 * m["D"]))
    asym = math.sqrt(PI_PHI) / (2 * math.sqrt(PI_Q))
    # nonlinear (background-vertex) decay phi_c phi_c -> D D where open: rate per condensate quantum
    # Gamma_nl ~ (h_D^2 phi_c)^2/(8 pi m_phi) sqrt(1 - 4 m_D^2/m_phi^2) (1 + 2 f_D(m_phi/2)), at phi_c = v
    Th = 32.1 * TC
    m = masses(Th)
    fD = 1.0 / math.expm1(0.5 * m["phi"] / Th)
    g_nl = (HD ** 2 * V) ** 2 / (8 * math.pi * m["phi"]) * math.sqrt(max(1 - 4 * m["D"] ** 2 / m["phi"] ** 2, 0)) \
        * (1 + 2 * fD)
    lin = 3.84e-5 * Th
    return dict(rows=rows, m_phi_over_2m_q_at_T_to_inf=asym,
                nonlinear_phi_c_to_DD_at_Hstar_production=dict(T_h_GeV=Th, Gamma_nl_GeV_at_phi_c_eq_v=g_nl,
                                                               Gamma_linear_2to2_GeV=lin, ratio=g_nl / lin),
                note="leading-order masses. phi -> q q: m_phi/(2 m_q) <= 0.946, closed at every T. phi -> D D is "
                     "kinematically open for T >~ 3 T_c (m_D grows only with Pi_D = 2.6e-4), but at phi = 0 there is "
                     "no cubic vertex (exact Z2): the decay exists only through the coherent background itself "
                     "(h_D^2 phi_c phi D D), i.e. at second order in phi_c. It adds damping while phi_c is large and "
                     "is negligible (ratio below) against the linear 2<->2 width; it is not credited.")


# ------------------------------------------------------------------------------------------------ tables
def gamma_table():
    eps = np.concatenate([[EPS_G, 2e-3, 5e-3], np.geomspace(1e-2, 1e3, 34)])
    rows = []
    t0 = time.time()
    for e in eps:
        Th = TC * (1 + e)
        gf, per, x = gamma_over_T(Th, "full", per_channel=True)
        gc = gamma_over_T(Th, "classical")
        gn = gamma_over_T(Th, "nsv")
        mphi = x["phi"] * Th
        rows.append(dict(eps=float(e), T_h_GeV=Th, m_phi_over_T=x["phi"], m_q_over_T=x["q"], m_D_over_T=x["D"],
                         Gamma_full_over_T=gf, Gamma_classical_over_T=gc, Gamma_nsv_over_T=gn,
                         Gamma_full_over_m_phi=gf * Th / mphi, per_channel_full_over_T=per))
    return rows, time.time() - t0


class GammaInterp:
    def __init__(self, rows, key):
        self.le = np.log(np.array([r["eps"] for r in rows]))
        self.lg = np.log(np.array([r[key] for r in rows]))

    def __call__(self, eps):
        le = np.log(np.clip(eps, math.exp(self.le[0]), math.exp(self.le[-1])))
        g = np.exp(np.interp(le, self.le, self.lg))
        # beyond the table at high T: Gamma/T approaches a constant (checked: slope -> 0)
        return g


def amp_rate(Gamma, m):
    """Decay rate of the coherent amplitude: Gamma/2 (underdamped) or the slow overdamped root."""
    if Gamma < 2 * m:
        return Gamma / 2
    return Gamma / 2 - math.sqrt(Gamma ** 2 / 4 - m ** 2)


def damping_exponent(xi, T_prod_SM, gfun, eps_cut=EPS_G, n=4000):
    """D = Int_{t_prod}^{t_exit,cut} r dt, dt = dT_SM/(H T_SM) (RD, constant g*)."""
    T_end = TC * (1 + eps_cut) / xi
    if T_prod_SM <= T_end:
        return 0.0
    lT = np.linspace(math.log(T_end), math.log(T_prod_SM), n)
    vals = np.empty(n)
    for i, l in enumerate(lT):
        T = math.exp(l)
        Th = xi * T
        eps = Th / TC - 1
        G = gfun(eps) * Th
        m = masses(Th)["phi"]
        vals[i] = amp_rate(G, m) / hubble(T)
    return float(np.trapezoid(vals, lT))


# ------------------------------------------------------------------------------------------------ KZ input
def kz_scale(xi):
    with open(os.path.join(TP_RES, "kibble_zurek.json")) as fh:
        kz = json.load(fh)
    xs = np.array([r["xi"] for r in kz["rows"]])
    L = np.array([r["z1_ising_ginzburg_matched"]["xi_KZ_GeV^-1"] for r in kz["rows"]])
    E = np.array([r["z1_ising_ginzburg_matched"]["eps_freeze"] for r in kz["rows"]])
    lx = math.log(xi)
    return (float(math.exp(np.interp(lx, np.log(xs), np.log(L)))),
            float(math.exp(np.interp(lx, np.log(xs), np.log(E)))))


C_W = (math.sqrt(math.pi) / 2 - (math.pi / 2) * math.e * erfc(1.0)) / (2 * math.pi ** 2)


def coherent_ratio(xi, T_prod_SM, D, phi0=V):
    """log10 r_damped, r_hubble_only and auxiliary numbers at the KZ decision point."""
    xiKZ, eps_hat = kz_scale(xi)
    m_hat = 1.0 / xiKZ
    Th0 = xi * T_prod_SM
    m0 = masses(Th0)["phi"]
    T_dec = TC * (1 + eps_hat)
    a_ratio = T_dec / Th0                        # a_0/a_dec = T_dec/T_0 (both baths scale as 1/a)
    A_lin = phi0 * a_ratio ** 1.5 * math.sqrt(m0 / m_hat)
    sigma = math.sqrt(C_W * T_dec * m_hat)
    # Hubble-only: cap where lambda A^2 = m^2 (quartic oscillation thereafter, amplitude ~ 1/a)
    Cc = phi0 ** 2 * (1.0 / Th0) ** 3 * m0        # A^2 = Cc T^3 / m(T)
    def gap(lT):
        T = math.exp(lT)
        mm = masses(T)["phi"]
        return LPHI * Cc * T ** 3 / mm - mm ** 2
    lo, hi = math.log(T_dec), math.log(Th0)
    if gap(hi) > 0:          # already nonlinear at production
        Tstar = Th0
    elif gap(lo) < 0:        # never nonlinear before the decision
        Tstar = None
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if gap(mid) > 0:
                lo = mid
            else:
                hi = mid
        Tstar = math.exp(0.5 * (lo + hi))
    if Tstar is None:
        A_hub = A_lin
    else:
        Astar = math.sqrt(Cc * Tstar ** 3 / masses(Tstar)["phi"])
        A_hub = Astar * T_dec / Tstar
    lr_d = math.log10(A_lin / sigma) - D / math.log(10)
    return dict(xi_KZ_GeV_inv=xiKZ, eps_hat=eps_hat, m_hat_GeV=m_hat, sigma_GeV=sigma, A_lin_undamped_GeV=A_lin,
                A_hubble_only_GeV=A_hub, T_nonlinear_GeV=Tstar, log10_r_damped=lr_d,
                log10_m0_damped=log10_m0(lr_d), r_hubble_only=A_hub / sigma,
                m0_hubble_only=math.erf(A_hub / sigma / math.sqrt(2)))


# ------------------------------------------------------------------------------------------------ dilute gas
C_DIL = (18 * LPHI ** 2 + 12 * GX ** 4 + 6 * HD ** 4) / (128 * math.pi)
C_EFF = (3 * LPHI + 2 * GX ** 2 + HD ** 2) / 4


def dilute_D(A_GeV):
    """Amplitude exponent for a free-streaming restoring gas, Pi = A H in RD, zero mode responding from early on:
    D = (C/(2 c_eff)) Int Pi/m dt = (C/(2 c_eff)) (A/(2 sqrt(l) v)) * 2 arctan(sqrt(u0 - 1)) -> pi/2 limit."""
    return C_DIL / (2 * C_EFF) * math.pi * A_GeV / (2 * math.sqrt(LPHI) * V)


def main():
    out = {"module": "residual_bias", "part": "A thermal restoration of a cold +v vacuum",
           "labels": "[D] derived, [C] computed, [Q] quoted, [E] estimate, [P] preparation (outside the action)"}
    out["parwani_check [C]"] = parwani_check()
    out["one_to_two_kinematics [C]"] = one_to_two_check()
    rows, dt = gamma_table()
    out["Gamma_table [C]"] = rows
    out["Gamma_table_runtime_s"] = dt
    hiT = rows[-1]
    out["Gamma_high_T [C]"] = dict(full=hiT["Gamma_full_over_T"], classical=hiT["Gamma_classical_over_T"],
                                   nsv=hiT["Gamma_nsv_over_T"],
                                   thermal_potential_hard_phi_rate_Q=4.981371739469606e-07,
                                   ratio_full_over_hard_rate=hiT["Gamma_full_over_T"] / 4.981371739469606e-07)
    gf = GammaInterp(rows, "Gamma_full_over_T")
    gc = GammaInterp(rows, "Gamma_classical_over_T")
    gn = GammaInterp(rows, "Gamma_nsv_over_T")

    # ---- H-star history [P]: xi = 5e-4, bath from T_SM = 1e8 GeV
    xi_s, Tp_s = 5e-4, 1e8
    Dfull = damping_exponent(xi_s, Tp_s, gf)
    Dcl = damping_exponent(xi_s, Tp_s, gc)
    Dns = damping_exponent(xi_s, Tp_s, gn)
    Dhard = damping_exponent(xi_s, Tp_s, lambda e: 4.981371739469606e-07 + 0 * e)
    T_exit = TC / xi_s
    star = dict(xi=xi_s, T_SM_prod_GeV=Tp_s, T_SM_exit_GeV=T_exit, H_exit_GeV=hubble(T_exit),
                Gamma_over_H_at_exit_full=gf(EPS_G * 10) * TC / hubble(T_exit),
                D_full=Dfull, D_classical=Dcl, D_nsv=Dns, D_with_hard_particle_rate_floor=Dhard,
                log10_suppression_full=-Dfull / math.log(10))
    for lab, D in (("full", Dfull), ("classical", Dcl), ("nsv", Dns), ("hard_rate_floor", Dhard),
                   ("full_x0.5", 0.5 * Dfull), ("full_x2", 2 * Dfull)):
        star["decision_" + lab] = coherent_ratio(xi_s, Tp_s, D)
    # Hubble-friction-only (no plasma damping)
    star["hubble_only"] = coherent_ratio(xi_s, Tp_s, 0.0)
    # production_map WKB estimate [Q], evaluated at H-star
    A_pi = LPHI * V ** 2 / hubble(T_exit)
    amp_pm = 1.7 * V * math.sqrt(LPHI) * V / A_pi
    star["production_map_WKB [Q formula]"] = dict(A_GeV=A_pi, coherent_amp_GeV=amp_pm,
                                                  soft_level_GeV=0.11 * math.sqrt(LPHI) * V,
                                                  ratio=amp_pm / (0.11 * math.sqrt(LPHI) * V),
                                                  m0_if_used=math.sqrt(2 / math.pi) * amp_pm / (0.11 * math.sqrt(LPHI) * V),
                                                  applicability="assumes the bath existed when H >> A (zero mode "
                                                  "frozen until H ~ A); for production at T_SM = 1e8 GeV (H << A) the "
                                                  "zero mode starts oscillating at production with amplitude ~ v")
    star["dilute_gas_formula_D (same A)"] = dilute_D(A_pi)
    out["H_star [P inputs]"] = star

    # ---- scan in xi and production temperature
    scan = []
    for xi in (1e-6, 3e-6, 1e-5, 2e-5, 3e-5, 4.5e-5, 6.6e-5, 8e-5, 1e-4, 1.5e-4, 2e-4, 3e-4, 5e-4, 7.6e-4, 1e-3,
               3e-3, 1e-2):
        for f in (1.001, 1.01, 1.1, 2.0, 10.0, 32.1, 1e3):
            Tp = TC * f / xi
            D = damping_exponent(xi, Tp, gf)
            cr = coherent_ratio(xi, Tp, D)
            scan.append(dict(xi=xi, T_h_prod_over_Tc=f, T_SM_prod_GeV=Tp, T_SM_exit_GeV=TC / xi, D_full=D,
                             log10_m0=cr["log10_m0_damped"], m0_hubble_only=cr["m0_hubble_only"],
                             r_hubble_only=cr["r_hubble_only"], kinetic_bath=xi >= 6.6e-5))
    out["scan_xi_Tprod [C]"] = scan

    # critical xi for thresholds (production long before exit, f = 32.1)
    def lm0(xi, f=32.1):
        Tp = TC * f / xi
        return coherent_ratio(xi, Tp, damping_exponent(xi, Tp, gf))["log10_m0_damped"]
    def xi_crit(target, f):
        lo, hi = 1e-6, 1e-1
        if lm0(hi, f) > target or lm0(lo, f) < target:
            return None
        for _ in range(40):
            mid = math.sqrt(lo * hi)
            if lm0(mid, f) > target:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)
    crit = {str(t): xi_crit(t, 32.1) for t in (-14, -17, -20, -25, -30)}
    out["xi_min_for_log10_m0_below [C]"] = crit
    out["xi_min_vs_production_factor [C]"] = {
        "f=T_h,prod/T_c=%g" % f: {str(t): xi_crit(t, f) for t in (-17, -20, -30)} for f in (1.01, 1.1, 2.0, 10.0)}

    # ---- dilute (non-thermal) gas: D(A) and the A, T_SM,end needed
    dil = dict(C_dil=C_DIL, c_eff=C_EFF, D_over_A_per_GeV=dilute_D(1.0),
               formula="D = (C/(2 c_eff)) * (pi/2) * A/(sqrt(lambda_phi) v), C = (18 l^2 + 12 g^4 + 6 h^4)/(128 pi)")
    need = {}
    for target in (-17, -20, -30):
        # need D >= ln(r_H sqrt(2/pi)) + |target| ln 10, with r_H ~ 1 (conservative lower) .. 1e2
        for rH in (1.0, 1e2):
            Dn = math.log(rH * math.sqrt(2 / math.pi)) + abs(target) * math.log(10)
            A = Dn / dilute_D(1.0)
            H_exit = LPHI * V ** 2 / A
            T_end = math.sqrt(H_exit * MP / math.sqrt(math.pi ** 2 * GSTAR / 90))
            need["log10m0=%d,r_H=%g" % (target, rH)] = dict(D_needed=Dn, A_min_GeV=A, T_SM_end_max_GeV=T_end)
    dil["requirements"] = need
    dil["production_map_B_min_edge_A~200GeV_D"] = dilute_D(200.0)
    out["dilute_gas [D, E]"] = dil
    dump("a_damping.json", out)
    return out


if __name__ == "__main__":
    t0 = time.time()
    o = main()
    s = o["H_star [P inputs]"]
    print("H-star D_full %.1f D_cl %.1f D_nsv %.1f D_hard %.1f" % (s["D_full"], s["D_classical"], s["D_nsv"],
                                                                   s["D_with_hard_particle_rate_floor"]))
    print("decision full", s["decision_full"])
    print("hubble only", s["hubble_only"])
    print("prod map", s["production_map_WKB [Q formula]"])
    print("crit", o["xi_min_for_log10_m0_below [C]"])
    print("dilute", o["dilute_gas [D, E]"])
    print("runtime %.1f s" % (time.time() - t0))
