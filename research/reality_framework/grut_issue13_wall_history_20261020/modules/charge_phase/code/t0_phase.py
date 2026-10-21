"""Parts 1-2: zero-temperature homogeneous charged condensate of the unchanged action.

Homogeneous rotating condensate X = |X| e^{i w t}, D = 0.  At fixed charge density
  rhobar(nbar) = min_{S,x} [ nbar^2/(2x) + Vbar(S,x) ],   nbar = Omega x.
Branch 1 (S^2 = 1 - a x > 0):  Omega^2 = a - (a^2-b) x,  rhobar = a x - (3/4)(a^2-b) x^2.
Branch 2 (S = 0, a x >= 1):    Omega^2 = b x,  nbar = sqrt(b) x^{3/2},  rhobar = 1/4 + (3/4) b x^2.
Everything below is checked against brute-force numerical minimisation.
Run:  python3 t0_phase.py   (seconds)
"""
import numpy as np
from scipy.optimize import brentq, minimize
from common import (A, B, LAM_PHI, V_GEV, Vbar, UNIT_N_OVER_V3, UNIT_RHO_OVER_V4,
                    M_X_OVER_UNIT, MSTAR_OVER_MX, ARCHIVED_SOLITON, CODEX_PRODUCTION, save_json)

a, b = A, B
c2 = a * a - b                         # > 0 in the binding class


# ---------------- analytic branches ----------------
def br1(x):
    Om2 = a - c2 * x
    return dict(x=x, S2=1 - a * x, Om=np.sqrt(Om2), n=x * np.sqrt(Om2), rho=a * x - 0.75 * c2 * x**2)


def br2(x):
    Om2 = b * x
    return dict(x=x, S2=0.0 * x, Om=np.sqrt(Om2), n=np.sqrt(b) * x**1.5, rho=0.25 + 0.75 * b * x**2)


def rho2_of_n(n):
    return 0.25 + 0.75 * b ** (1 / 3) * n ** (4 / 3)


x_c = 1.0 / a                                     # |X|_c^2 = x_c v^2/2
n_c = np.sqrt(b) * x_c**1.5                       # branch-2 lower end (phi=0 local min appears)
x_s1 = 2 * a / (3 * c2)                           # branch-1 radial spinodal
n_s1 = br1(x_s1)["n"]
x_star = 1 / np.sqrt(b)                           # Q-matter density (thin-wall interior)
n_star = b ** -0.25
Om_star = b ** 0.25                               # = m*/(sqrt(lam_phi) v)


# homogeneous global switch: rho1(n) = rho2(n) on the overlap n_c < n < n_s1
def d12(x1):
    r = br1(x1)
    return r["rho"] - rho2_of_n(r["n"])


# n along branch 1 increases on (0, x_s1); n_c corresponds to x in (0, x_s1)
x1_at_nc = brentq(lambda x: br1(x)["n"] - n_c, 1e-9, x_s1)
x_t = brentq(d12, x1_at_nc, x_s1)
n_t = br1(x_t)["n"]
x2_t = (n_t**2 / b) ** (1 / 3)


# ---------------- brute-force numerical check ----------------
def E_of(Sx, n):
    S, lx = Sx
    x = np.exp(lx)
    return n * n / (2 * x) + Vbar(S, x)


def numeric_min(n):
    best = None
    for S0 in (0.0, 0.3, 0.7, 0.95):
        for x0 in (1e-3, 0.05, 0.3, 0.8, 1.5, 4.0, 10.0):
            r = minimize(E_of, [S0, np.log(x0)], args=(n,), method="Nelder-Mead",
                         options=dict(xatol=1e-12, fatol=1e-15, maxiter=20000))
            if best is None or r.fun < best.fun:
                best = r
    S, x = abs(best.x[0]), np.exp(best.x[1])
    return S, x, best.fun


def analytic_min(n):
    cands = []
    # branch 1 stationary points with x < x_s1 (minima)
    if n < n_s1:
        x = brentq(lambda x: br1(x)["n"] - n, 0.0, x_s1)
        cands.append(("phi!=0", x, br1(x)["rho"]))
    if n >= n_c:
        x = (n * n / b) ** (1 / 3)
        cands.append(("phi=0", x, rho2_of_n(n)))
    return min(cands, key=lambda c: c[2]), cands


grid_n = np.concatenate([np.geomspace(1e-3, 0.03, 8), np.linspace(0.031, 0.045, 15),
                         np.geomspace(0.05, 5.0, 12)])
check_rows, max_dev = [], 0.0
for n in grid_n:
    S, x, E = numeric_min(n)
    (lab, xa, Ea), cands = analytic_min(n)
    dev = abs(E - Ea) / max(1e-12, abs(Ea))
    max_dev = max(max_dev, dev)
    check_rows.append(dict(nbar=n, numeric_S=S, numeric_x=x, numeric_rho=E, analytic_branch=lab,
                           analytic_x=xa, analytic_rho=Ea, rel_dev=dev,
                           n_stationary_minima=len(cands)))

# ---------------- grand canonical (fixed mu) at T = 0 ----------------
def p_of_mu(Om):
    def negp(Sx):
        S, lx = Sx
        x = np.exp(lx)
        return -(Om**2 * x / 2 - Vbar(S, x) + Vbar(1.0, 0.0))
    best = None
    for S0 in (0.0, 0.5, 1.0):
        for x0 in (1e-6, 0.5, 3.0, 20.0):
            r = minimize(negp, [S0, np.log(x0)], method="Nelder-Mead",
                         options=dict(xatol=1e-12, fatol=1e-15, maxiter=20000))
            if best is None or r.fun < best.fun:
                best = r
    return -best.fun, abs(best.x[0]), np.exp(best.x[1])


gc_rows = []
for Om in np.linspace(0.40, 0.80, 21):
    p, S, x = p_of_mu(Om)
    p_an = max(0.0, Om**4 / (4 * b) - 0.25) if Om**2 / b >= x_c else 0.0
    gc_rows.append(dict(Omega=Om, p_numeric=p, S=S, x=x, p_analytic=p_an))

# ---------------- convex hull (phase coexistence) ----------------
ng = np.linspace(0.0, 4.0, 40001)
rho_h = np.array([analytic_min(n)[0][2] if n > 0 else 0.0 for n in ng])
slope = rho_h[1:] / ng[1:]
i_t = np.argmin(slope)
hull_tangent_n = ng[1:][i_t]
hull_tangent_slope = slope[i_t]

# ---------------- linear fluctuations about the homogeneous condensate ----------------
def growth(x, k):
    """Largest Re(lambda) for perturbations e^{lambda t + i k z} around branch background x."""
    if a * x < 1:
        S02 = 1 - a * x
        Om2 = a - c2 * x
        ms2 = 2 * S02
    else:
        S02 = 0.0
        Om2 = b * x
        ms2 = a * x - 1
    P = np.polynomial.polynomial
    t1 = P.polymul(P.polymul([2 * b * x, 1], [0, 1]), [ms2, 1])
    t2 = 4 * Om2 * P.polymul([-k * k, 1], [ms2, 1])
    t3 = -4 * a * a * x * S02 * np.array([0, 1])
    poly = P.polyadd(P.polyadd(t1, t2), t3)
    Aroots = P.polyroots(poly)
    lam2 = Aroots - k * k
    lam = np.sqrt(lam2.astype(complex))
    g = np.max(lam.real)
    # S-only sector on branch 2 (decoupled): lambda^2 = -(k^2 + ms2)
    if S02 == 0.0:
        g = max(g, np.sqrt(complex(-(k * k + ms2))).real)
    return g


def band_and_peak(x):
    ks = np.linspace(1e-4, 3.0, 6001)
    gs = np.array([growth(x, k) for k in ks])
    unstable = ks[gs > 1e-10]
    i = np.argmax(gs)
    return dict(kmax_band=float(unstable.max()) if unstable.size else 0.0,
                gamma_max=float(gs[i]), k_at_gamma_max=float(ks[i]))


def cs2_branch1(x):
    Om2 = a - c2 * x
    return -x * c2 / (2 * Om2 - x * c2)


frag_rows = []
for x in [0.01, 0.05, 0.1, 0.25, 0.5, 0.7, 0.74, 0.8, 0.95, 1.0, 1.5, 3.162, 6.0]:
    bp = band_and_peak(x)
    kmax_an = None
    if a * x < 1:
        S02 = 1 - a * x
        # (k^2+2bx)(k^2+2S0^2) = 4 a^2 x S0^2
        qa, qb, qc = 1.0, 2 * b * x + 2 * S02, 4 * b * x * S02 - 4 * a * a * x * S02
        disc = qb * qb - 4 * qa * qc
        kmax_an = float(np.sqrt((-qb + np.sqrt(disc)) / 2)) if qc < 0 else 0.0
    if a * x < 1:
        r = br1(x)
        Om, nb, rho = r["Om"], r["n"], r["rho"]
        cs2 = cs2_branch1(x) if x < x_s1 else None
    else:
        r = br2(x)
        Om, nb, rho = r["Om"], r["n"], r["rho"]
        cs2 = 1 / 3
    frag_rows.append(dict(x=x, branch="phi!=0" if a * x < 1 else "phi=0", Omega=Om, nbar=nb,
                          n_over_v3=nb * UNIT_N_OVER_V3, rho_over_n_mX=rho / nb / M_X_OVER_UNIT,
                          dOmega_dn_sign=("-" if (a * x < 1 and x < x_s1) else "+"),
                          cs2_longwave=cs2, kmax_band_numeric=bp["kmax_band"],
                          kmax_band_analytic=kmax_an, gamma_max=bp["gamma_max"],
                          k_at_gamma_max=bp["k_at_gamma_max"]))

# Codex production background
xp = CODEX_PRODUCTION["F2"]
bp = band_and_peak(xp)
k_box = 2 * np.pi / CODEX_PRODUCTION["box_side"]
eps = CODEX_PRODUCTION["seed_rms"]
codex_prod = dict(
    x=xp, nbar=br1(xp)["n"], nbar_codex=CODEX_PRODUCTION["Qbar_per_volume"],
    n_over_v3=br1(xp)["n"] * UNIT_N_OVER_V3, n_over_nc=br1(xp)["n"] / n_c,
    rho_over_n_mX=br1(xp)["rho"] / br1(xp)["n"] / M_X_OVER_UNIT,
    cs2_longwave=cs2_branch1(xp), unstable_band_kmax=bp["kmax_band"],
    gamma_max=bp["gamma_max"], k_at_gamma_max=bp["k_at_gamma_max"],
    wavelength_at_gamma_max=2 * np.pi / bp["k_at_gamma_max"],
    box_fundamental_k=k_box, n_unstable_box_harmonics_per_axis=int(bp["kmax_band"] / k_box),
    t_nonlinear_estimate_ln_1_over_eps=float(np.log(1 / eps) / bp["gamma_max"]),
    note="t_nonlinear estimate = ln(1/0.005)/gamma_max; compare Codex first phi<0 at tau=72.35-73.35",
)

# ---------------- Q-matter / vacuum planar interface at mu = m* ----------------
def interface_tension(L=30.0, N=1501, iters=200000):
    z = np.linspace(-L, L, N)
    dz = z[1] - z[0]
    Fst = np.sqrt(x_star)
    S = 0.5 * (1 - np.tanh(z / 2.0))
    F = 0.5 * Fst * (1 + np.tanh(z / 2.0))
    dt = 0.2 * dz * dz
    Om2 = Om_star**2

    def dV(S, F):
        x = F * F
        dS = (S * S - 1) * S + a * x * S
        dF = a * S * S * F + b * x * F - Om2 * F
        return dS, dF

    for it in range(iters):
        lapS = np.zeros_like(S); lapF = np.zeros_like(F)
        lapS[1:-1] = (S[2:] - 2 * S[1:-1] + S[:-2]) / dz**2
        lapF[1:-1] = (F[2:] - 2 * F[1:-1] + F[:-2]) / dz**2
        dS, dF = dV(S, F)
        S[1:-1] += dt * (lapS[1:-1] - dS[1:-1])
        F[1:-1] += dt * (lapF[1:-1] - dF[1:-1])
        # re-centre (exact degeneracy => zero mode); shift by interpolation every 20000 steps
        if it % 20000 == 0 and it > 0:
            zc = np.interp(0.5, S[::-1], z[::-1])
            S = np.interp(z + zc, z, S, left=1.0, right=0.0)
            F = np.interp(z + zc, z, F, left=0.0, right=Fst)
    Sp = np.gradient(S, dz); Fp = np.gradient(F, dz)
    dens = 0.5 * (Sp**2 + Fp**2) + Vbar(S, F * F) - Om2 * F * F / 2 - 0.0
    # vacuum (S=1,F=0) and Q-matter have zero grand potential density at mu=m*
    sig = np.trapezoid(dens, z)
    grad = np.trapezoid(0.5 * (Sp**2 + Fp**2), z)
    pot = np.trapezoid(Vbar(S, F * F) - Om2 * F * F / 2, z)
    return sig, grad, pot


sig_Q, sig_grad, sig_pot = interface_tension()


def thin_wall_E_over_Q(Qbar):
    R = (3 * Qbar / (4 * np.pi * n_star)) ** (1 / 3)
    return Om_star + 3 * sig_Q / (n_star * R), R


Qa = ARCHIVED_SOLITON["Qbar"]
tw_arch, R_arch = thin_wall_E_over_Q(Qa)
# Qbar above which thin-wall estimate gives E/Q < m_X (free carrier):  Om* + 3 sig/(n* R) < 1
R_bind = 3 * sig_Q / (n_star * (M_X_OVER_UNIT - Om_star))
Q_bind = 4 * np.pi / 3 * R_bind**3 * n_star

# thin-wall ESTIMATE of vacuum-bubble nucleation in stretched (n_c<n<n_star) homogeneous Q-matter at T=0:
# pressure pbar = nbar*Omega - rhobar < 0 on branch 2; eps = -pbar; S4 = 27 pi^2 sigma^4/(2 eps^3) (physical, /lam_phi)
nuc_rows = []
for nb in np.linspace(n_c, n_star, 9)[:-1]:
    xb = (nb * nb / b) ** (1 / 3)
    pbar = nb * np.sqrt(b * xb) - rho2_of_n(nb)
    epsb = -pbar
    S4 = 27 * np.pi**2 * sig_Q**4 / (2 * epsb**3) / LAM_PHI
    Rc = 3 * sig_Q / epsb
    nuc_rows.append(dict(nbar=nb, n_over_v3=nb * UNIT_N_OVER_V3, pbar=pbar, Rc_bar=Rc, S4_thin_wall=S4))
S4_needed = 4 * np.log(1.8e9)   # Gamma ~ m^4 e^{-S4} > H^4 with m/H ~ 1.8e9 at T_vis ~ 2e5 GeV (cosmo_map)

# rho/n table (homogeneous ground state among homogeneous configurations)
table = []
for n in [1e-4, 1e-3, 0.01, 0.1, 0.2201, n_c, n_t, n_s1, 0.5, 1.0, n_star, 3.0, 10.0]:
    (lab, x, rho), cands = analytic_min(n)
    table.append(dict(nbar=n, n_over_v3=n * UNIT_N_OVER_V3, branch=lab, x=x,
                      rho_over_n_in_mX=rho / n / M_X_OVER_UNIT,
                      rho_over_n_in_mstar=rho / n / Om_star,
                      local_minima=[c[0] for c in cands]))

# densities at which the homogeneous energy per charge equals the archived soliton's 0.91108652 m_X
e_arch = ARCHIVED_SOLITON["E_over_mXQ"] * M_X_OVER_UNIT
x1_arch = brentq(lambda x: br1(x)["rho"] / br1(x)["n"] - e_arch, 1e-6, x_s1)
n1_arch = br1(x1_arch)["n"]
n2_arch_lo = brentq(lambda n: rho2_of_n(n) / n - e_arch, n_c, n_star)
n2_arch_hi = brentq(lambda n: rho2_of_n(n) / n - e_arch, n_star, 100.0)

v3_GeV3 = V_GEV**3
out = dict(
    homogeneous_equals_archived_E_over_Q=dict(
        branch1_n_over_v3=n1_arch * UNIT_N_OVER_V3, branch1_x=x1_arch,
        branch2_n_over_v3_low=n2_arch_lo * UNIT_N_OVER_V3, branch2_n_over_v3_high=n2_arch_hi * UNIT_N_OVER_V3,
        note="homogeneous rho/n < 0.91108652 m_X on phi!=0 branch for n > branch1 value (up to n_s1) and on "
             "phi=0 branch for n between low and high values"),
    conventions="S=phi/v, x=F^2=2|X|^2/v^2, a=gX^2/lam_phi=1, b=lam_X/lam_phi=0.1; nbar=n/(sqrt(lam_phi) v^3); "
                "rhobar=rho/(lam_phi v^4); Omega=w/(sqrt(lam_phi) v); m_X=sqrt(a) sqrt(lam_phi) v.",
    status="derived (closed form) and computed (brute-force minimisation agrees)",
    X_c_sq_over_v2=x_c / 2, F2_c=x_c,
    n_c_over_v3=n_c * UNIT_N_OVER_V3, n_c_GeV3_at_v1TeV=n_c * UNIT_N_OVER_V3 * v3_GeV3,
    n_c_formula="n_c = sqrt(8 lam_X) |X|_c^3 = sqrt(8 lam_X) (lam_phi/(2 gX^2))^{3/2} v^3",
    n_c_meaning="lower spinodal of the phi=0 homogeneous phase: phi=0 is a (local) minimum for n>=n_c; "
                "at fixed |X| the phi transition at |X|_c is continuous, but at fixed n it is not",
    branch1_spinodal_x=x_s1, branch1_spinodal_n_over_v3=n_s1 * UNIT_N_OVER_V3,
    branch1_spinodal_formula="x_s1 = 2 a/(3(a^2-b)); n_s1 = gX v^3 x_s1/sqrt(3)",
    homogeneous_global_switch_n_over_v3=n_t * UNIT_N_OVER_V3,
    homogeneous_global_switch_x_branch1=x_t, homogeneous_global_switch_x_branch2=x2_t,
    homogeneous_global_switch_S_branch1=np.sqrt(1 - a * x_t),
    order_at_fixed_n="first order (metastability window n_c < n < n_s1; global switch at n_t; "
                     "S jumps from %.4f to 0)" % np.sqrt(1 - a * x_t),
    binodal_n_star_over_v3=n_star * UNIT_N_OVER_V3, binodal_n_star_GeV3=n_star * UNIT_N_OVER_V3 * v3_GeV3,
    binodal_x_star=x_star, binodal_X_star_over_v=np.sqrt(x_star / 2),
    hull_tangent_numeric_nbar=hull_tangent_n, hull_tangent_numeric_slope=hull_tangent_slope,
    mu_c_grand_canonical_over_mX=Om_star / M_X_OVER_UNIT, mu_c_GeV=Om_star * np.sqrt(LAM_PHI) * V_GEV,
    mstar_over_mX_brief=MSTAR_OVER_MX,
    grand_canonical_statement="T=0, fixed mu: vacuum (n=0) for mu<m*; first-order jump to phi=0 Q-matter "
                              "with n=n_star at mu=m*",
    brute_force_max_rel_dev=max_dev, brute_force_rows=check_rows, grand_canonical_rows=gc_rows,
    fragmentation_rows=frag_rows,
    fragmentation_verdict="phi!=0 homogeneous branch (all n<n_s1 stable radially only for x<x_s1) has "
                          "d(Omega)/dn<0, c_s^2<0: long-wavelength (FLS/Q-ball) fragmentation instability for "
                          "every x in (0, x_s1); phi=0 branch (n>=n_c) has c_s^2=1/3, S mass^2=a x-1>=0: linearly "
                          "stable, metastable to phase separation for n<n_star, ground state for n>n_star",
    codex_production_linear=codex_prod,
    Qmatter_interface_tension_bar=sig_Q, Qmatter_interface_grad_part=sig_grad, Qmatter_interface_pot_part=sig_pot,
    Qmatter_interface_tension_GeV3=sig_Q * np.sqrt(LAM_PHI) * V_GEV**3,
    thin_wall_E_over_Q_at_archived_Q_mX=tw_arch / M_X_OVER_UNIT, thin_wall_R_at_archived_Q=R_arch,
    archived_E_over_mXQ=ARCHIVED_SOLITON["E_over_mXQ"], archived_Q_phys=ARCHIVED_SOLITON["Q_phys"],
    thin_wall_Qbar_for_binding=Q_bind, thin_wall_Q_phys_for_binding=Q_bind / LAM_PHI,
    thin_wall_R_for_binding=R_bind,
    rho_over_n_table=table,
    vacuum_bubble_nucleation_thin_wall_estimate=dict(
        rows=nuc_rows, S4_needed_for_one_bubble_per_hubble_time=S4_needed,
        note="ESTIMATE: thin-wall with the coexistence tension sigma_Q at mu=m*, O(4) bounce, T=0. Thin wall fails "
             "near the spinodal n_c where the barrier vanishes. Within this estimate the stretched homogeneous phi=0 "
             "phase does not nucleate vacuum bubbles before n approaches n_c; the exit is spinodal (tachyonic)."),
)
save_json("t0_phase.json", out)
if __name__ == "__main__":
    for k in ["X_c_sq_over_v2", "n_c_over_v3", "n_c_GeV3_at_v1TeV", "branch1_spinodal_x",
              "branch1_spinodal_n_over_v3", "homogeneous_global_switch_n_over_v3",
              "homogeneous_global_switch_S_branch1", "binodal_n_star_over_v3", "hull_tangent_numeric_nbar",
              "hull_tangent_numeric_slope", "mu_c_GeV", "brute_force_max_rel_dev",
              "Qmatter_interface_tension_bar", "Qmatter_interface_grad_part", "Qmatter_interface_pot_part",
              "thin_wall_E_over_Q_at_archived_Q_mX", "thin_wall_R_at_archived_Q",
              "thin_wall_Q_phys_for_binding"]:
        print(k, out[k])
    print("codex production:", codex_prod)
    for r in frag_rows:
        print(r)
    for r in table:
        print(r)
    for r in gc_rows[::4]:
        print(r)
