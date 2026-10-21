#!/usr/bin/env python3
"""v_production independent check, parts (1)-(3): kinetic T_c(z), MB exit coefficient, T_SM,min.

Written from scratch by the verifier (no import of the module's code).
Hartree curvature at phi = 0 with a common fugacity z = exp(mu_N/T_h) for the four real hidden dof
(phi, q1, q2, D):  Pi = 3 lam_phi <phi^2> + g_X^2 (<q1^2> + <q2^2>) + h_D^2 <D^2>,
<chi^2> = int d^3k/(2pi)^3 f/E, f = 1/(z^-1 e^{E/T} - 1).  Thermal parts only (no CW), as in the module.
"""
import json, math, os
import numpy as np
from scipy import integrate, optimize, special

os.environ.setdefault("OMP_NUM_THREADS", "2")

v = 1000.0; gX = 0.1; lphi = 0.01; lX = 0.001; hD = 0.01; lD = 0.001; mD0 = 100.0
gstar = 106.75
Ymax = 4.373e-12
zeta2 = math.pi**2 / 6; zeta3 = special.zeta(3); zeta4 = math.pi**4 / 90


def polylog(s, z):
    """Li_s(z) for 0 <= z <= 1 via series (fast for z<=0.5) or integral representation."""
    if z == 0:
        return 0.0
    if z <= 0.5:
        k = np.arange(1, 400)
        return float(np.sum(z**k / k**s))
    # Bose integral: Li_s(z) = 1/Gamma(s) int_0^inf x^{s-1}/(z^-1 e^x - 1) dx
    f = lambda x: x**(s - 1) / (np.exp(x) / z - 1.0)
    val, _ = integrate.quad(f, 0, 60, limit=400, epsabs=0, epsrel=1e-13)
    val2, _ = integrate.quad(f, 60, np.inf, limit=200)
    return (val + val2) / math.gamma(s)


def I_th(m, T, z):
    """<chi^2> thermal for one real dof, Bose, fugacity z (mu = T ln z <= 0 <= m)."""
    a = m / T
    def f(x):
        E = math.sqrt(x * x + a * a)
        if E == 0.0:
            return 0.0
        return x * x / (E * (math.exp(E) / z - 1.0))
    v1, _ = integrate.quad(f, 0, 50, limit=400, epsabs=0, epsrel=1e-12)
    v2, _ = integrate.quad(f, 50, np.inf, limit=200)
    return T * T * (v1 + v2) / (2 * math.pi**2)


def solve_gap(T, z):
    """Fixed-point solve of m_q^2, m_D^2 at phi=0 with m_phi = 0 (boundary). Returns Pi, mq, mD."""
    Iphi = I_th(0.0, T, z)
    mq2 = gX**2 * Iphi; mD2 = mD0**2 + hD**2 * Iphi
    for _ in range(200):
        Iq = I_th(math.sqrt(mq2), T, z)
        ID = I_th(math.sqrt(mD2), T, z)
        mq2n = gX**2 * Iphi + 4 * lX * Iq
        mD2n = mD0**2 + hD**2 * Iphi + 3 * lD * ID
        if abs(mq2n - mq2) < 1e-13 * mq2n and abs(mD2n - mD2) < 1e-13 * mD2n:
            mq2, mD2 = mq2n, mD2n
            break
        mq2, mD2 = mq2n, mD2n
    Iq = I_th(math.sqrt(mq2), T, z); ID = I_th(math.sqrt(mD2), T, z)
    Pi = 3 * lphi * Iphi + 2 * gX**2 * Iq + hD**2 * ID
    return Pi, math.sqrt(mq2), math.sqrt(mD2)


def Tc_hartree(z):
    lead = v * math.sqrt(lphi / ((3 * lphi + 2 * gX**2 + hD**2) * polylog(2, z) / (2 * math.pi**2)))
    g = lambda T: solve_gap(T, z)[0] - lphi * v * v
    T = optimize.brentq(g, 0.9 * lead, 1.2 * lead, xtol=1e-12, rtol=1e-13)
    Pi, mq, mD = solve_gap(T, z)
    return T, lead, mq / T, mD / T


out = {}
# ---- (1) leading-order and Hartree T_c(z)
c_sum = 3 * lphi + 2 * gX**2 + hD**2
Tc1_lead = v * math.sqrt(lphi / (c_sum / 12))
out["leading_Pi_over_T2_at_z1"] = c_sum / 12
out["Tc1_leading_over_v"] = Tc1_lead / v
rows = []
Tc1_H = None
for z in [1.0, 0.9, 0.5, 0.1, 0.01, 1e-6]:
    T, lead, xq, xD = Tc_hartree(z)
    if z == 1.0:
        Tc1_H = T
    rows.append(dict(z=z, Li2=polylog(2, z), Tc_hartree_over_v=T / v, Tc_leading_over_v=lead / v,
                     factor_leading=math.sqrt(zeta2 / polylog(2, z)),
                     ratio_hartree_to_hartree_z1=T / Tc1_H,
                     hartree_over_leading=T / lead, mq_over_T=xq, mD_over_T=xD))
out["Tc_of_z"] = rows

# ---- (2) MB-regime exit coefficient: Pi = c_eff n_h/(2 T_h), rho_h = 3 n_h T_h
c_eff = c_sum / 4
s_over_T3 = 2 * math.pi**2 * gstar / 45
rho_over_T4 = math.pi**2 * gstar / 30
K = 1.5 * c_eff * s_over_T3**2 / rho_over_T4
coef = v * math.sqrt(lphi / K)
out["MB"] = dict(c_eff=c_eff, s_over_T3=s_over_T3, rho_over_T4=rho_over_T4, K=K,
                 Tend_coefficient_GeV=coef,
                 formula="T_SM,end = sqrt(lam_phi v^2 eps/(K Y^2)) = coef*sqrt(eps)/Y_N")

# exact-Bose (massless, leading) vs MB ratio of T_end at fixed (Y, eps): T_end ∝ Pi^-1/2,
# Pi_exact/Pi_MB = Li2 Li4 / Li3^2
mbrows = []
for lnz in [0.0, -0.5, -1.0, -2.0, -3.0, -5.0]:
    z = math.exp(lnz)
    r = polylog(2, z) * polylog(4, z) / polylog(3, z)**2
    mbrows.append(dict(mu_over_T=lnz, Pi_exact_over_MB=r, Tend_exact_over_MB=r**-0.5))
out["MB"]["exact_vs_MB_leading"] = mbrows


# ---- (3) mu_N = 0 exit and T_SM,min; G(z) = T_end Y^{1/3}, F(z) = eps/Y^{4/3} using Hartree T_c(z)
def nh_rho_over_Th(z):
    n = 4 * polylog(3, z) / math.pi**2      # n_h / T_h^3 (massless, 4 dof)
    r = 12 * polylog(4, z) / math.pi**2     # rho_h / T_h^4
    return n, r


def GF(z, Tc):
    n, r = nh_rho_over_Th(z)
    # n Th^3 = Y s_over_T3 T^3  -> xi = (Y s/n)^(1/3);  T_end = Tc/xi  => G = Tc (n/s)^(1/3)
    G = Tc * (n / s_over_T3)**(1 / 3)
    # eps = r xi^4/rho_over_T4 = r (Y s/n)^(4/3)/rho  => F = r (s/n)^(4/3)/rho
    F = r * (s_over_T3 / n)**(4 / 3) / rho_over_T4
    return G, F


G1, F1 = GF(1.0, Tc1_H)
out["mu0"] = dict(G1_GeV=G1, F1=F1, Tc1_hartree_over_v=Tc1_H / v,
                  xi_exit_at_Ymax=(Ymax * s_over_T3 / nh_rho_over_Th(1.0)[0])**(1 / 3),
                  TSMmin_hartree=G1 * Ymax**(-1 / 3),
                  TSMmin_resummed_range=[1.5571 * v * (nh_rho_over_Th(1.0)[0] / s_over_T3)**(1 / 3) * Ymax**(-1 / 3),
                                         1.5641 * v * (nh_rho_over_Th(1.0)[0] / s_over_T3)**(1 / 3) * Ymax**(-1 / 3)],
                  TSMmin_leading=Tc1_lead * (nh_rho_over_Th(1.0)[0] / s_over_T3)**(1 / 3) * Ymax**(-1 / 3))
# monotonicity of G(z): G must increase as z decreases for mu=0 to be the minimum
Gs = []
for row in rows:
    G, F = GF(row["z"], row["Tc_hartree_over_v"] * v)
    Gs.append(dict(z=row["z"], G=G, F=F, MB_closed_form_Tend_over_Y13=coef * math.sqrt(F)))
out["G_F_of_z"] = Gs
out["G_monotone_increasing_as_z_decreases"] = all(Gs[i + 1]["G"] > Gs[i]["G"] for i in range(len(Gs) - 1))

# check the exit-table rows that are quoted: (Y, eps/eps_th) -> T_end, using exact Bose with Hartree Tc(z)
def exit_T(Y, eps_ratio):
    eps = eps_ratio * F1 * Y**(4 / 3)
    target = eps / Y**(4 / 3)
    # find z with F(z) = target (F with leading Tc is z-independent of Tc)
    h = lambda lnz: GF(math.exp(lnz), 1.0)[1] - target
    if eps_ratio == 1:
        lnz = 0.0
    else:
        lnz = optimize.brentq(h, -80, -1e-12)
    z = math.exp(lnz)
    Tc = Tc_hartree(z)[0]
    return dict(Y=Y, eps_over_eps_th=eps_ratio, mu_over_T=lnz, Tend=GF(z, Tc)[0] * Y**(-1 / 3),
                Tend_MB_closed_form=coef * math.sqrt(eps) / Y)

out["exit_rows"] = [exit_T(4.373e-12, 1), exit_T(4.373e-12, 2), exit_T(4.373e-12, 10), exit_T(4.373e-12, 1e4),
                    exit_T(4.373e-12, 1e6), exit_T(1e-12, 1), exit_T(1e-13, 1), exit_T(1e-14, 1)]

json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vp_thermal.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
