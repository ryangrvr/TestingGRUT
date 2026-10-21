"""Item 3: hidden heating and vacuum randomisation from inflation-era gravity with zero direct hidden branching.
All inflaton/inflation quantities (H_inf, H_e, m_I, duration N_obs) are OUTSIDE-ACTION preparation data.

(i)  Inflaton-condensate gravitational scattering phi_I phi_I -> h -> (hidden pairs) for a quadratic minimum:
     energy rate per real minimally coupled scalar dof = (1-6 xi)^2 m_I rho_I^2/(128 pi M_P^4)
     [derived here via the conformal-frame effective mass -(1/6 - xi) R with R oscillating at 2 m_I; agrees with
     Clery et al. 2112.15214 eqs.(22),(55),(62); Mambrini-Olive 2102.06214 eq.(7) as fetched is a factor 8 lower].
     Numerical mode-function check of the golden-rule pair rate A^2/(64 pi).
(ii) Stochastic (Starobinsky-Yokoyama astro-ph/9407016) equilibrium of the light hidden spectators in de Sitter:
     P_eq ∝ exp(-8 pi^2 U/(3 H^4)) for (phi, q1, q2, D) with U + curvature terms xi_a(mu=H) R chi_a^2/2, R=12 H^2.
     Marginal F(phi), barrier exponent, sign-relaxation rate Lambda_1 (Fokker-Planck eigenvalue and Kramers
     integral, checked against SY eq.(58)), equilibrium energy <U> = c H^4 (formal T_h,max) and vacuum
     randomisation thresholds.
"""
import numpy as np
from scipy.special import erfcx
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from common import (M_P, V_BENCH, LAM_PHI, LAM_X, LAM_D, C_PORT, D_PORT, M_D0, G_H, C_AVG, xi_running, tc,
                    T_of_rho, dump, TC_BAND)


# ----------------------------------------------------------------- (i) condensate scattering
def mode_function_check(A0=0.001, omega=2.0, tf=200.0, ramp=30.0, nk=1601, dt=0.005):
    """Solve sigma_k'' + (k^2 + A(t) cos(omega t)) sigma_k = 0 for a k-grid around omega/2 (RK4, vectorised);
    compare n(t_f) = int d^3k/(2pi)^3 |beta_k|^2 with the golden rule A0^2 T_eff/(32 pi) (particles of one real
    scalar; pairs = half)."""
    k = np.linspace(0.6 * omega / 2, 1.4 * omega / 2, nk)
    s = np.exp(-1j * 0) / np.sqrt(2 * k) * np.ones_like(k, dtype=complex)
    ds = -1j * k * s

    def Aof(t):
        x = np.clip(t / ramp, 0, 1)
        y = np.clip((tf - t) / ramp, 0, 1)
        return A0 * np.sin(0.5 * np.pi * x)**2 * np.sin(0.5 * np.pi * y)**2

    def acc(t, s):
        return -(k**2 + Aof(t) * np.cos(omega * t)) * s

    t = 0.0
    nsteps = int(tf / dt)
    for _ in range(nsteps):
        k1s, k1d = ds, acc(t, s)
        k2s, k2d = ds + 0.5 * dt * k1d, acc(t + 0.5 * dt, s + 0.5 * dt * k1s)
        k3s, k3d = ds + 0.5 * dt * k2d, acc(t + 0.5 * dt, s + 0.5 * dt * k2s)
        k4s, k4d = ds + dt * k3d, acc(t + dt, s + dt * k3s)
        s = s + dt * (k1s + 2 * k2s + 2 * k3s + k4s) / 6
        ds = ds + dt * (k1d + 2 * k2d + 2 * k3d + k4d) / 6
        t += dt
    beta2 = (np.abs(ds)**2 + k**2 * np.abs(s)**2) / (2 * k) - 0.5
    n_num = np.trapezoid(k**2 * beta2, k) / (2 * np.pi**2)
    tt = np.linspace(0, tf, 200001)
    Teff = np.trapezoid((Aof(tt) / A0)**2, tt)
    n_gr = A0**2 * Teff / (32 * np.pi)
    return {"n_numeric": n_num, "n_golden_rule": n_gr, "ratio": n_num / n_gr, "A0": A0, "omega": omega, "t_f": tf}


def condensate_peak(H_e, m_I, N_h=4, xi=0.0):
    """Peak hidden energy (MD, no decay): rho_h(x) = 18 N_h m H^3/(128 pi) (1-6xi)^2 x^4 (1 - sqrt x), x = a_e/a."""
    pref = 18 * N_h * m_I * H_e**3 * (1 - 6 * xi)**2 / (128 * np.pi)
    x = (8.0 / 9.0)**2
    return pref * x**4 * (1 - np.sqrt(x)), 1 / x


def condensate_Pi_threshold(m_over_H, lv2, era="MD"):
    """chi2(a) = (9 N_h 0.4 /(128 pi)) H_e^3/m (a_e/a)^2 (sum over 4 dof), Pi = C_avg chi2.
    Hubble criterion: exists a with Pi - lambda v^2 > H^2; MD: H^2 = H_e^2 x^-3; RD-after-a_e (upper case):
    H^2 = H_e^2 x^-4."""
    P0c = C_AVG * 9 * 4 * 0.4 / (128 * np.pi)

    def margin(H_e):
        m = m_over_H * H_e
        x = np.logspace(0, 30, 6000)
        Pi = P0c * H_e**3 / m * x**-2.0
        H2 = H_e**2 * (x**-3.0 if era == "MD" else x**-4.0)
        return np.max((Pi - lv2) / H2)
    return 10**brentq(lambda l: np.log(max(margin(10**l), 1e-300)), 0.0, 19.0)


def part_condensate():
    mc = mode_function_check()
    mc["note"] = ("A0 = 0.01 gives ratio 1.064 (t_f=200) growing with t_f: parametric (Bose) enhancement, Floquet "
                  "mu = A/(4k); at A0 = 0.001 the residual (mu t)^2/3 ~ 6e-4 matches the golden rule")
    out = {"mode_function_check": mc}
    v = V_BENCH
    lv2 = LAM_PHI * v**2
    rows = []
    for mH in [1.0, 3.0, 10.0]:
        row = {"m_I/H_e": mH}
        for lab, f in [("nominal", 1.0), ("minus30", 1 - TC_BAND), ("plus30", 1 + TC_BAND)]:
            T_c = tc(v, f)
            g = lambda l: np.log(T_of_rho(condensate_peak(10**l, mH * 10**l)[0], G_H) / T_c)
            row[f"H_e_max_formal_thermal_{lab}"] = 10**brentq(g, -2, 19)
        row["H_e_max_nonthermal_Pi_MD"] = condensate_Pi_threshold(mH, lv2, "MD")
        row["H_e_max_nonthermal_Pi_RD"] = condensate_Pi_threshold(mH, lv2, "RD")
        rows.append(row)
    out["bounds"] = rows
    rho, ap = condensate_peak(1.0, 1.0)
    out["peak_rho_h_coefficient"] = rho
    out["peak_rho_h_formula"] = "rho_h,peak = 7.76e-3 (m_I/H_e) H_e^4 at a = 1.27 a_e (N_h = 4, xi = 0)"
    out["peak_scale_factor_over_a_e"] = ap
    out["formal_Th_peak_coeff"] = "T_h,peak = 0.277 (m_I H_e^3)^(1/4)"
    # Clery eq.(61) cross-check: same process into the 4 Higgs dof thermalised with g* = 106.75
    rho_end = 1e64
    He = np.sqrt(rho_end / 3) / M_P
    rho_peak, _ = condensate_peak(He, 3e13)
    out["clery_eq61_check"] = {"T_max_Higgs_ours": T_of_rho(rho_peak, 106.75), "T_max_Higgs_clery": 3.0e12,
                               "note": "equal energy goes to the 4 hidden dof (both minimally coupled scalar quartets)"}
    out["hidden_equals_higgs_note"] = ("With xi_H = xi_hidden = 0 at mu_0 the condensate deposits the same energy in "
                                       "the 4 hidden real scalars as in the 4 Higgs components (gauge bosons and "
                                       "massless fermions are conformal at tree level).")
    return out


# ----------------------------------------------------------------- (ii) stochastic spectators
def log_erfcx(z):
    z = np.asarray(z, float)
    out = np.empty_like(z)
    m = z > -5
    out[m] = np.log(erfcx(z[m]))
    zz = z[~m]
    out[~m] = zz**2 + np.log(2 - erfcx(-zz) * np.exp(-zz**2))
    return out


def U_phi(phi, H, xi_phi):
    return LAM_PHI / 4 * (phi**2 - V_BENCH**2)**2 + 0.5 * xi_phi * 12 * H**2 * phi**2


def logZq(phi, a, H, xi_X):
    """ln of int d^2q exp(-a[(1/2) M r^2 + (lambda_X/4) r^4]), M = c phi^2 + xi_X R (alpha = 0)."""
    M = C_PORT * phi**2 + xi_X * 12 * H**2
    A = a * M / 2.0
    B = a * LAM_X / 4.0
    return np.log(np.pi) + np.log(0.5 * np.sqrt(np.pi / B)) + log_erfcx(A / (2 * np.sqrt(B)))


def logZD(phi, a, H, xi_D, nD=601):
    phi = np.atleast_1d(phi)
    M = M_D0**2 + D_PORT * phi**2 + xi_D * 12 * H**2
    Mmin = M.min()
    scale = max(np.sqrt(1.0 / (a * abs(Mmin))) if Mmin > 0 else 0.0, (1.0 / (a * LAM_D))**0.25)
    if Mmin < 0:
        scale = max(scale, np.sqrt(-Mmin / LAM_D))
    D = np.linspace(-12 * scale, 12 * scale, nD)
    e = -a * (0.5 * M[:, None] * D[None, :]**2 + LAM_D / 4 * D[None, :]**4)
    emax = e.max(axis=1)
    return emax + np.log(np.trapezoid(np.exp(e - emax[:, None]), D, axis=1))


def marginal_F(H, phis, include_XD=True, include_xi=True):
    a = 8 * np.pi**2 / (3 * H**4)
    xi = xi_running(H) if include_xi else np.zeros(3)
    F = U_phi(phis, H, xi[0])
    if include_XD:
        F = F - (logZq(phis, a, H, xi[1]) + logZD(phis, a, H, xi[2])) / a
    return F - F.min()


def fp_lambda1(x, Fx_over_3H4):
    """Lowest nonzero FP eigenvalue per e-fold for dx/dN = -u'(x) + (1/2pi) noise, u = F/(3H^4), x = phi/H."""
    Dd = 1.0 / (8 * np.pi**2)
    u = Fx_over_3H4
    h = x[1] - x[0]
    up = np.gradient(u, h)
    upp = np.gradient(up, h)
    W = (up / (2 * Dd))**2 - upp / (2 * Dd)
    diag = Dd * (2.0 / h**2 + W)
    off = -Dd / h**2 * np.ones(len(x) - 1)
    w = eigh_tridiagonal(diag, off, select="i", select_range=(0, 1), eigvals_only=True)
    return w[1] - w[0], w[0]


def kramers_lambda1(x, Fx, H):
    """SY eq.(57) integral form: Lambda_1/H = 1/(4 pi^2 N_x I_x), N_x = int e^{-aF} dx, I_x = int_0^{x0} e^{aF} dx."""
    a = 8 * np.pi**2 / (3 * H**4)
    e1 = -a * Fx
    m1 = e1.max()
    lnN = m1 + np.log(np.trapezoid(np.exp(e1 - m1), x))
    i0 = np.argmin(np.abs(x))
    imin = i0 + np.argmin(Fx[i0:])
    xs, Fs = x[i0:imin + 1], Fx[i0:imin + 1]
    e2 = a * Fs
    m2 = e2.max()
    lnI = m2 + np.log(np.trapezoid(np.exp(e2 - m2), xs)) if len(xs) > 2 else -np.inf
    return np.exp(-np.log(4 * np.pi**2) - lnN - lnI)


def stochastic_point(H, include_XD=True, include_xi=True, nx=4001):
    v = V_BENCH
    a = 8 * np.pi**2 / (3 * H**4)
    # field range: cover both wells and the equilibrium tail
    xi = xi_running(H) if include_xi else np.zeros(3)
    mu2 = LAM_PHI * v**2 - xi[0] * 12 * H**2          # effective tachyonic coefficient
    phimin_est = np.sqrt(max(mu2, 0) / LAM_PHI)
    tail = (40.0 / (a * LAM_PHI))**0.25
    L = max(phimin_est + 6 * tail, 3 * v, 8 * tail)
    phis = np.linspace(-L, L, nx)
    F = marginal_F(H, phis, include_XD, include_xi)
    i0 = nx // 2
    imin = i0 + np.argmin(F[i0:])
    dF = F[i0] - F[imin]
    S = a * dF
    x = phis / H
    lam_fp, lam0 = fp_lambda1(x, F / (3 * H**4))
    lam_kr = kramers_lambda1(x, F, H)
    # SY eq.(58) analytic (pure double well, no X/D, no xi)
    m2 = LAM_PHI * v**2
    lam_sy = np.sqrt(2) * m2 / (3 * np.pi * H**2) * np.exp(-2 * np.pi**2 * m2**2 / (3 * LAM_PHI * H**4))
    # equilibrium energy <U_true> (curvature terms removed after inflation) via 3D quadrature on (phi, r, D)
    U_true_mean, U_eq_mean = equilibrium_energy(H, a, xi, np.linspace(-L, L, 801))
    return {"H": H, "H_over_v": H / v, "phi_min_over_v": phis[imin] / v, "barrier_exponent_S": S,
            "Lambda1_per_efold_FP": lam_fp, "FP_lambda0_residual": lam0, "Lambda1_per_efold_Kramers": lam_kr,
            "Lambda1_SY_eq58_pure_double_well": lam_sy, "U_true_mean": U_true_mean, "U_true_over_H4": U_true_mean / H**4,
            "U_eq_over_H4": U_eq_mean / H**4, "T_h_formal_from_U": T_of_rho(U_true_mean, G_H),
            "P_phi_negative": 0.5}


def equilibrium_energy(H, a, xi, phis, nr=600, nD=301):
    """<U_true> and <U_eq> under P_eq ∝ exp(-a U_eq) (alpha = 0: q and D integrals factorise at fixed phi).
    U_true excludes the curvature terms (energy left after inflation; gradients negligible)."""
    v = V_BENCH
    R = 12 * H**2
    Uphi_t = LAM_PHI / 4 * (phis**2 - v**2)**2
    Uphi_e = Uphi_t + 0.5 * R * xi[0] * phis**2
    # q sector on (phi, r)
    rmax = max((60 / (a * LAM_X))**0.25, 3 * np.sqrt(max(-xi[1] * R, 0) / LAM_X) + 1e-12)
    r = np.linspace(0, rmax, nr)
    P2, R2 = np.meshgrid(phis, r, indexing="ij")
    Uq_t = C_PORT * P2**2 * R2**2 / 2 + LAM_X / 4 * R2**4
    Uq_e = Uq_t + 0.5 * R * xi[1] * R2**2
    eq = -a * Uq_e
    eq -= eq.max(axis=1, keepdims=True)
    wq = np.exp(eq) * R2
    lnZq = np.log(wq.sum(axis=1)) + (-a * Uq_e).max(axis=1)
    mq_t = (wq * Uq_t).sum(axis=1) / wq.sum(axis=1)
    mq_e = (wq * Uq_e).sum(axis=1) / wq.sum(axis=1)
    # D sector on (phi, D)
    Dmax = max((60 / (a * LAM_D))**0.25, np.sqrt(60 / (a * M_D0**2)))
    D = np.linspace(-Dmax, Dmax, nD)
    P3, D3 = np.meshgrid(phis, D, indexing="ij")
    UD_t = 0.5 * (M_D0**2 + D_PORT * P3**2) * D3**2 + LAM_D / 4 * D3**4
    UD_e = UD_t + 0.5 * R * xi[2] * D3**2
    ed = -a * UD_e
    edm = ed.max(axis=1, keepdims=True)
    wd = np.exp(ed - edm)
    lnZD = np.log(wd.sum(axis=1)) + edm[:, 0]
    mD_t = (wd * UD_t).sum(axis=1) / wd.sum(axis=1)
    mD_e = (wd * UD_e).sum(axis=1) / wd.sum(axis=1)
    lw = -a * Uphi_e + lnZq + lnZD
    w = np.exp(lw - lw.max())
    Z = w.sum()
    Ut = float((w * (Uphi_t + mq_t + mD_t)).sum() / Z)
    Ue = float((w * (Uphi_e + mq_e + mD_e)).sum() / Z)
    return Ut, Ue


def part_stochastic(N_obs=60.0):
    out = {"N_obs_prep_data": N_obs}
    Hs = np.concatenate([np.linspace(0.08, 0.6, 27), [0.7, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]]) * V_BENCH
    rows = [stochastic_point(H) for H in Hs]
    rows_pure = [stochastic_point(H, include_XD=False, include_xi=False) for H in Hs[:27:3]]
    out["rows"] = rows
    out["rows_pure_double_well_check"] = rows_pure
    # checks: Kramers vs FP where both are reliable; SY eq58 vs Kramers for the pure double well
    chk = []
    for r in rows_pure:
        if 1e-8 < r["Lambda1_per_efold_FP"] and r["barrier_exponent_S"] > 3:
            chk.append({"H_over_v": r["H_over_v"], "S": r["barrier_exponent_S"], "FP": r["Lambda1_per_efold_FP"],
                        "Kramers": r["Lambda1_per_efold_Kramers"], "SY58": r["Lambda1_SY_eq58_pure_double_well"]})
    out["lambda1_method_checks"] = chk
    # thresholds (use Kramers when S>8, FP otherwise)
    def lam(r):
        return r["Lambda1_per_efold_Kramers"] if r["barrier_exponent_S"] > 8 else r["Lambda1_per_efold_FP"]
    lHs = np.log10(Hs)
    lam_arr = np.array([lam(r) for r in rows])
    llam = np.log10(lam_arr)

    def H_at(target):
        # first crossing as H increases
        for i in range(len(Hs) - 1):
            if (llam[i] - target) * (llam[i + 1] - target) <= 0:
                t = (target - llam[i]) / (llam[i + 1] - llam[i])
                return 10**(lHs[i] + t * (lHs[i + 1] - lHs[i]))
        return None
    thr = {}
    thr["H_Lambda1_Nobs_eq_1_raw"] = H_at(np.log10(1.0 / N_obs))
    thr["H_Lambda1_Nobs_eq_0.1"] = H_at(np.log10(0.1 / N_obs))
    for Dl in [5.0, 10.0, 20.0]:
        thr[f"H_Ndom_Delta{int(Dl)}_eq_1"] = H_at(np.log10(2 * 3 / (np.exp(3 * Dl) - 1)))
    out["thresholds"] = thr
    out["threshold_definitions_note"] = ("Lambda1-based thresholds apply only in the bimodal regime H < H_restoration; "
                                          "above it the end-of-inflation phi is a unimodal, nearly scale-invariant random "
                                          "field whose zero set is a wall network on all scales")
    out["threshold_definitions"] = {
        "Lambda1": "sign relaxation rate per e-fold (odd FP mode); sign correlation ~ exp(-Lambda1 dN)",
        "N_dom(Delta)": "(Lambda1/2)(e^{3 Delta}-1)/3 flip events in our patch creating domains that exited within "
                        "Delta e-folds of the observable-universe scale (heuristic, O(1))"}
    # formal thermal bound from <U_true>
    U = np.array([r["U_true_mean"] for r in rows])
    Th = T_of_rho(U, G_H)
    form = {}
    for lab, f in [("nominal", 1.0), ("minus30", 1 - TC_BAND), ("plus30", 1 + TC_BAND)]:
        T_c = tc(V_BENCH, f)
        lt = np.log10(Th) - np.log10(T_c)
        hit = None
        for i in range(len(Hs) - 1):
            if lt[i] * lt[i + 1] <= 0:
                t = -lt[i] / (lt[i + 1] - lt[i])
                hit = 10**(lHs[i] + t * (lHs[i + 1] - lHs[i]))
                break
        form[lab] = hit
    out["H_max_formal_thermal_from_dS_energy"] = form
    out["c_of_H4_large_H"] = rows[-1]["U_true_over_H4"]
    out["c_of_H4_range"] = [min(r["U_true_over_H4"] for r in rows), max(r["U_true_over_H4"] for r in rows)]
    # sign-change probability between points separated by dN e-folds (labelled approximations)
    for row in rows:
        lam1 = row["Lambda1_per_efold_FP"] if row["barrier_exponent_S"] <= 8 else row["Lambda1_per_efold_Kramers"]
        if row["barrier_exponent_S"] > 0:
            row["regime"] = "bimodal"
            row["P_sign_differs_dN1"] = 0.5 * (1 - np.exp(-lam1))
            row["P_sign_differs_dN60"] = 0.5 * (1 - np.exp(-lam1 * N_obs))
        else:
            row["regime"] = "unimodal (Z2 restored in dS)"
            row["P_sign_differs_dN1"] = np.arccos(np.exp(-lam1)) / np.pi     # Gaussian level-set estimate
            row["P_sign_differs_dN60"] = np.arccos(np.exp(-lam1 * N_obs)) / np.pi
    # equilibrium de Sitter restoration of Z2 by the light X fluctuations: marginal F has its minimum at phi=0
    def curv0(H):
        ph = np.array([0.0, 0.02 * V_BENCH])
        F = marginal_F(H, ph)
        return F[1] - F[0]
    Hr = brentq(curv0, 0.25 * V_BENCH, 1.0 * V_BENCH, xtol=1e-6)
    out["H_origin_locally_stable"] = Hr
    out["H_origin_locally_stable_over_v"] = Hr / V_BENCH

    def glob(H):
        ph = np.linspace(0, 2.0 * V_BENCH, 4001)
        F = marginal_F(H, ph)
        return F[0] - F[1:].min() - 1e-12 * H**4
    Hg = brentq(glob, Hr, 1.0 * V_BENCH, xtol=1e-6)
    out["H_restoration_equilibrium"] = Hg
    out["H_restoration_equilibrium_over_v"] = Hg / V_BENCH
    out["restoration_note"] = ("between H_origin_locally_stable and H_restoration_equilibrium the marginal F(phi) has "
                               "minima at 0 and at +-phi_b (trimodal P_eq); above, phi = 0 is the global minimum")
    Hr = Hg
    raw = thr["H_Lambda1_Nobs_eq_1_raw"]
    thr["H_vacuum_randomised_min(H_restoration, Lambda1_Nobs=1)"] = min(Hr, raw) if raw else Hr
    # analytic leading estimate: Pi_phi = g_X^2 <q1^2+q2^2>, <r^2>_eq = H^2 sqrt(3/(2 pi^3 lambda_X)) (quartic X)
    r2eq = np.sqrt(3.0 / (2 * np.pi**3 * LAM_X))
    out["X_variance_eq_over_H2"] = r2eq
    out["H_restoration_leading_estimate_over_v"] = np.sqrt(LAM_PHI / (C_PORT * r2eq))
    # finite duration: free diffusion from X=0 for N e-folds, <q1^2+q2^2> = N H^2/(2 pi^2) (until equilibrium)
    fd = {}
    for N in [30.0, 60.0, 100.0]:
        r2 = min(N / (2 * np.pi**2), r2eq)
        fd[f"N={int(N)}"] = np.sqrt(LAM_PHI / (C_PORT * r2))
    out["H_restoration_over_v_after_N_efolds_free_diffusion"] = fd
    out["X_equilibration_efolds_estimate"] = 1.0 / (1.36859 * np.sqrt(LAM_X / (24 * np.pi**2)))
    return out


def run():
    out = {"condensate": part_condensate(), "stochastic": part_stochastic()}
    out["labels"] = "H_inf, H_e, m_I, N_obs, initial field: outside-action preparation data"
    dump("c4_inflation_sources.json", out)
    return out


if __name__ == "__main__":
    import json
    r = run()
    print(json.dumps(r["condensate"], indent=1, default=float))
    st = r["stochastic"]
    for row in st["rows"]:
        print({k: (round(v, 5) if isinstance(v, float) and abs(v) > 1e-3 else v) for k, v in row.items()})
    print(st["lambda1_method_checks"])
    print(st["thresholds"])
    print(st["H_max_formal_thermal_from_dS_energy"], st["c_of_H4_large_H"])
