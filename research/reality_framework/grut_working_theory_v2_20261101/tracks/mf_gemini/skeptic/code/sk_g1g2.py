"""mf_gemini SKEPTIC: independent re-derivation of the G1/G2 numbers, the loop criteria, the Higgs routes,
the holder algebra and two loopholes.  Own code; nothing imported from mf_gemini/code.

Methods differ from the track where possible:
  * loop criterion: RENORMALISATION CONDITIONS (V(0) = V0, V(v) = 0, V'(v) = 0 [+ V''(v) = m^2, V'(0) = 0 when
    odd counterterms are allowed]) instead of a least-squares counterterm fit; the measure is the induced change of
    the tachyonic curvature at chi = 0 (which sets the waterfall rate) and of m_chi^2, plus max|dV|/V0 on [0, v];
  * Higgs routes: both portal signs, threshold and Landau pole, and the Higgs loop on chi with the track's own 1%
    least-squares criterion;
  * species-bound check of the large-N_f loophole;
  * the -v domain of the holder on (chi - v) with M1r's heavy-sector mass formula (light Y mass^2 = M_Y^2 - sqrt2 |mu| h);
  * K_crit with the Y deposits also rescaled (K^{1/2}, inf_potential skeptic scenario B).
Units: GeV, reduced M_P.  Run: python3 -B sk_g1g2.py -> ../results/sk_g1g2.json
"""
import json
import os
import numpy as np
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'sk_g1g2.json')
MP = 2.435e18


def proposal(H):
    """m = sqrt(2 lam) v = 3H, lam v^4/4 = 3 MP^2 H^2  ->  solve numerically (not by the closed form)."""
    # unknowns (v, lam): 2 lam v^2 = 9 H^2 ; lam v^4 = 12 MP^2 H^2 -> divide: v^2/2 = 12 MP^2/9
    f = lambda lv: [2 * np.exp(lv[1]) * np.exp(2 * lv[0]) - 9 * H**2, np.exp(lv[1]) * np.exp(4 * lv[0]) - 12 * MP**2 * H**2]
    from scipy.optimize import fsolve
    sol = fsolve(f, [np.log(MP), np.log(1e-14)], xtol=1e-14)
    v, lam = np.exp(sol[0]), np.exp(sol[1])
    return v, lam


def cw_conditions(kind, y, lam, nodd=False):
    """Fermion (Majorana, 2 dof) CW on chi in units x = chi/v, V0 = lam v^4/4 = 1:
    dV = -(y^4/(8 pi^2 lam)) f(x),  f = x^4 [ln(y^2 x^2) - 3/2]   (kind 'Z2', m_N = y chi)
                                    f = (x-1)^4 [ln(y^2 (x-1)^2) - 3/2] (kind 'shift', m_N = y (chi - v))
    Counterterms: Z2 {1, x^2, x^4} fixed by V(0) = 1, V(1) = 0, V'(1) = 0.
                  shift {1, x, x^2, x^3, x^4} fixed by V(0) = 1, V'(0) = 0, V(1) = 0, V'(1) = 0, V''(1) = 8.
    Returns relative change of V''(0) (tree -4), of V''(1) (tree 8), and max |V - V_tree| on [0, 1]."""
    x = np.linspace(1e-9, 1.0, 20001)
    pref = y**4 / (8 * np.pi**2 * lam)

    def f(xx, d=0):
        # numerical derivatives by complex-step-free central differences on a fine analytic form
        if kind == 'Z2':
            g = lambda t: -pref * (t**4 * (np.log(y**2 * t**2) - 1.5))
        else:
            g = lambda t: -pref * np.where(np.abs(t - 1) > 0, (t - 1)**4 * (np.log(y**2 * (t - 1)**2 + 1e-300) - 1.5), 0.0)
        h = 1e-5
        if d == 0:
            return g(xx)
        if d == 1:
            return (g(xx + h) - g(xx - h)) / (2 * h)
        return (g(xx + h) - 2 * g(xx) + g(xx - h)) / h**2
    tree = lambda t: (t**2 - 1)**2
    if kind == 'Z2':
        # V = tree + f + a + b x^2 + c x^4 ; conditions at x->0 (f(0)=0), x = 1
        A = np.array([[1, 0, 0], [1, 1, 1], [0, 2, 4]], float)
        rhs = np.array([0.0, -f(1.0), -f(1.0, 1)])
        a, b, c = np.linalg.solve(A, rhs)
        ct = lambda t: a + b * t**2 + c * t**4
        ct2 = lambda t: 2 * b + 12 * c * t**2
    else:
        # unknowns a0..a4 ; conditions V(0)=1 (tree 1), V'(0)=0, V(1)=0, V'(1)=0, V''(1)=8 (tree values)
        A = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [1, 1, 1, 1, 1], [0, 1, 2, 3, 4], [0, 0, 2, 6, 12]], float)
        rhs = -np.array([f(0.0), f(0.0, 1), 0.0, 0.0, 0.0])   # f and derivatives up to 3 vanish at x = 1
        co = np.linalg.solve(A, rhs)
        ct = lambda t: sum(co[i] * t**i for i in range(5))
        ct2 = lambda t: sum(i * (i - 1) * co[i] * t**(i - 2) for i in range(2, 5))
    dV = f(x) + ct(x)
    d2_0 = (f(1e-4, 2) if kind == 'Z2' else f(0.0, 2)) + ct2(0.0)
    d2_1 = f(1.0 - 1e-4, 2) + ct2(1.0) if kind == 'shift' else f(1.0, 2) + ct2(1.0)
    return dict(dV2_at_0_over_tree=float(d2_0 / -4.0), dV2_at_v_over_tree=float(d2_1 / 8.0),
                max_absdV_over_V0=float(np.max(np.abs(dV))))


def higgs_cw_lsq(lchiH, lam):
    """Higgs loop on chi for the OPPOSITE sign, m_H^2 = |l| (v^2 - chi^2) >= 0 on [0, v), 4 real dof:
    dV/V0 = (l^2/(4 pi^2 lam)) (1-x^2)^2 [ln(1-x^2) - 3/2] (mu = sqrt(|l|) v); track's criterion: least-squares
    fit of {1, x^2, x^4} on the grid, max residual.  Grid [0, 0.999] (the IR log at x -> 1 is excluded)."""
    x = np.linspace(0, 0.999, 20001)
    u = 1 - x**2
    dV = (lchiH**2 / (4 * np.pi**2 * lam)) * u**2 * (np.log(u) - 1.5)
    A = np.vstack([x**0, x**2, x**4]).T
    c, *_ = np.linalg.lstsq(A, dV, rcond=None)
    return float(np.max(np.abs(dV - A @ c)))


def main():
    out = {}
    # ---------------- G1 ----------------
    g1 = {}
    for b in [2, 3, 4]:
        # oscillator at rest with all energy potential: m^2 Phi^2/2 = 3 MP^2 H^2
        Phi = brentq(lambda P: 0.5 * b**2 * P**2 - 3.0, 1e-3, 10.0)
        g1[f'm={b}H'] = float(Phi)
    # Majorana width y^2 m beta^3 / 16 pi (beta = 1 unblocked); Gamma = 0.5 H at m = 3H
    y = brentq(lambda yy: yy**2 * 3.0 / (16 * np.pi) - 0.5, 0.1, 10)
    out['G1'] = dict(Phi_over_MP=g1, y_Majorana_Gamma_half_H=y,
                     blocking_yPhi_over_mhalf=float(y * g1['m=3H'] * MP / (1.5 * 1.3e11)),
                     CW_over_tree=float(y**4 * (g1['m=3H'] * MP)**4 / (32 * np.pi**2) / (3 * MP**2 * (1.3e11)**2)))
    # ---------------- G2 ----------------
    g2 = {}
    for H in [1.25e11, 1.3e11, 1.55e11]:
        v, lam = proposal(H)
        V0 = lam * v**4 / 4
        d = dict(v_over_MP=v / MP, lam=lam, V0_quarter_GeV=V0**0.25, H_check=np.sqrt(V0 / 3) / MP,
                 tach_mass2_at_0_over_H2=-lam * v**2 / H**2, m_over_H=np.sqrt(2 * lam) * v / H,
                 MN_over_mhalf_y2=2 * v / (1.5 * H))
        for yy in [1.5, 2.0]:
            d[f'y={yy}'] = dict(Z2=cw_conditions('Z2', yy, lam), shift=cw_conditions('shift', yy, lam))
        # caps: |relative change of V''(0)| = 1% (waterfall rate) ; and max|dV| = 1% V0
        cap = {}
        for kind in ['Z2', 'shift']:
            for key in ['dV2_at_0_over_tree', 'max_absdV_over_V0']:
                fn = (lambda ly, k=kind, kk=key: np.log(abs(cw_conditions(k, np.exp(ly), lam)[kk]) + 1e-300) - np.log(0.01))
                try:
                    cap[f'{kind}:{key}'] = float(np.exp(brentq(fn, np.log(1e-6), np.log(0.5))))
                except ValueError:
                    cap[f'{kind}:{key}'] = None
        d['y_cap_1pct'] = cap
        ymax = max(c for c in cap.values() if c)
        d['Gamma_over_H_at_largest_cap_unblocked'] = ymax**2 * 3 / (16 * np.pi)
        g2[f'H={H:.3g}'] = {k: (float(v_) if isinstance(v_, (float, np.floating)) else v_) for k, v_ in d.items()}
    out['G2'] = g2
    # ---------------- Higgs portal (both signs) at H = 1.3e11 ----------------
    H = 1.3e11
    v, lam = proposal(H)
    m = 3 * H
    hp = {}
    for G in [0.5, 1.0]:
        A = np.sqrt(8 * np.pi * m * G * H)          # 4 real dof: Gamma = A^2/(8 pi m)
        l = A / (2 * v)
        thr = l**2 / lam                             # same for both signs: lam_H(UV) = lam_SM + l^2/lam
        tLP = 16 * np.pi**2 / (24 * thr)             # one loop, 24 lam^2/16pi^2
        hp[f'Gamma={G}H'] = dict(lam_chiH=l, threshold=thr, LandauPole_GeV=m * np.exp(tLP),
                                 LP_over_V0quarter=m * np.exp(tLP) / (lam * v**4 / 4)**0.25,
                                 Higgs_CW_on_chi_lsq_residual_over_V0_opposite_sign=higgs_cw_lsq(l, lam),
                                 lam_chi_running_per_efold_over_lam=l**2 / (2 * np.pi**2) / lam,
                                 valley_H_at_overshoot_GeV=np.sqrt(lam * ((1.147 * v)**2 - v**2) / (2 * l)),
                                 note='valley_H: |H| on the flat direction chi^2 - v^2 = 2|l||H|^2/lam at the overshoot')
    out['Higgs'] = {k: ({kk: float(vv) if not isinstance(vv, str) else vv for kk, vv in d.items()}) for k, d in hp.items()}
    # ---------------- large-N_f loophole (perturbative Yukawa, bare mass tuned) ----------------
    ycap = g2['H=1.3e+11']['y_cap_1pct']['shift:max_absdV_over_V0'] or 5.5e-4
    # Gamma_tot = N_f y^2 m/16pi = 0.5 H ; loop ~ N_f y^4 <= ycap^4
    Nf = (8.0 * np.pi / 3.0 * 2) ** 2 / ycap**4 / 4      # placeholder overwritten below (keep exact algebra)
    Ny2 = 0.5 * 16 * np.pi / 3.0                          # N_f y^2 = 8.38
    Nf = Ny2**2 / ycap**4
    yNf = np.sqrt(Ny2 / Nf)
    out['largeNf'] = dict(Nf_y2=Ny2, y_cap_used=ycap, Nf_min=Nf, y=yNf,
                          blocking_at_overshoot=yNf * 0.147 * v / (m / 2),
                          species_cutoff_GeV=MP / np.sqrt(Nf), species_cutoff_over_H=MP / np.sqrt(Nf) / H)
    # ---------------- holder on (chi - v): the -v domain ----------------
    MY, h, mu = 30.0, 40.0, 7.159456
    u_inf = 10j
    Gv = mu - u_inf
    mu_minus = u_inf - Gv
    mYL2 = lambda mm: MY**2 - np.sqrt(2) * abs(mm) * h
    out['holder_minus_v'] = dict(mu_plus=[mu, 0.0], mu_minus=[mu_minus.real, mu_minus.imag], abs_mu_minus=abs(mu_minus),
                                 mYL2_plus_H2=mYL2(mu), mYL_plus=np.sqrt(mYL2(mu)), mYL2_minus_H2=mYL2(mu_minus),
                                 mYL2_inflation_H2=mYL2(u_inf), mu_threshold_tachyonic=MY**2 / (np.sqrt(2) * h),
                                 x_where_minus_path_turns_tachyonic=float(brentq(
                                     lambda xx: abs(u_inf - Gv * xx) - MY**2 / (np.sqrt(2) * h), 0.0, 1.0)),
                                 note='light Y mass^2 = M_Y^2 - sqrt2 |mu| h (reproduces omega_YL = 22.25 H at mu = 7.16); '
                                      'path mu(chi) = 10i + G chi along chi in [0, -v]: x = |chi|/v')
    # ---------------- K_crit with the Y deposits rescaled ----------------
    W = ((3.4716 + 0.00428) / 0.190255)**(1 / 0.75013) + 0.97494
    R = W / 15.6036 - 1
    Kc_A = 1.98581 * np.sqrt((R - 0.379 - 0.0186 - 4e-4) / 3.608)
    Kc_B = 1.98581 * brentq(lambda xx: 3.608 * xx**2 + 0.379 * np.sqrt(xx) + 0.0186 + 4e-4 - R, 0.1, 2.0)
    out['K_crit'] = dict(W=W, R_total=R, K_crit_Rvec_only=Kc_A, K_crit_Rvec_and_Ydeposit=Kc_B)
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
