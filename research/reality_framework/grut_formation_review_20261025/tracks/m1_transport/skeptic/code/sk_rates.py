"""Skeptic check of m1_transport T4-T7, T9-T11: own background, rates, free streaming, scalings.

Own code (not derived from the reviewer's transport_check.py nor Codex's free_stream.py).
Background: photons + e+- (exact Fermi-Dirac integrals, entropy conservation of the gamma-e plasma
after instantaneous nu decoupling), N_eff = 3.044 neutrinos, Omega_m h^2 = 0.143, flat Lambda.
Above 2 MeV a piecewise g_* table (10.75 / 14.25 / 17.25 / 61.75 / 75.75 / 86.25) is used only for
the formation-clock scaling test.  Codex's checkpoint distribution (late_kin_h0025_96.npz, 6dc9556)
is read for a labelled reproduction; everything else is Maxwellian/analytic.
"""
import json, math, os
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

RES = Path(__file__).resolve().parents[1] / 'results'
SP = Path(os.environ.get('SP', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')) / 'review' / 'm1_transport_skeptic'

kB = 8.617333262e-5 * 1e-9          # GeV/K
T0 = 2.7255 * kB                      # GeV
Mpl = 2.435e18
hbar = 6.582119569e-25                # GeV s
c_kms = 299792.458
Mpc_km = 3.0856775814913673e19
pc_per_Mpc = 1e6
h = 0.674; Omh2 = 0.143; Om = Omh2 / h**2; Neff = 3.044
H0 = 100 * h / Mpc_km * hbar          # GeV
rho_crit = 3 * H0**2 * Mpl**2
me = 0.51099895e-3
mX = 0.09240896373753621
YQ = 5.079414232725842e-9
Qref = 2e5
Nplus = 13189.13146074248
f_charge = Nplus / Qref
sigma = (4e-3 - 0.04 / 3)**2 / (128 * math.pi)       # confirmed in sk_xsec.py
vrms = 0.06837528835889453
fgas = 0.058767411685782034


def g_e(y):
    """(g_rho, g_s) of e+- at m_e/T = y (4 dof)."""
    E = lambda x: math.sqrt(x * x + y * y)
    ir = quad(lambda x: x * x * E(x) / (math.exp(E(x)) + 1), 0, 80, epsrel=1e-11, epsabs=0)[0]
    ip = quad(lambda x: x**4 / (3 * E(x)) / (math.exp(E(x)) + 1), 0, 80, epsrel=1e-11, epsabs=0)[0]
    grho = 60 / math.pi**4 * ir
    gs = 45 / (2 * math.pi**2) * (2 / math.pi**2) * (ir + ip)
    return grho, gs


Ts = np.geomspace(2e-6, 0.02, 500)
GE = np.array([g_e(me / T) for T in Ts])
gs_pl = 2 + GE[:, 1]; gr_pl = 2 + GE[:, 0]
a_of_Ts = T0 / Ts * (2 / gs_pl)**(1 / 3)                 # gamma-e plasma entropy conservation
lnT_of_lna = PchipInterpolator(np.log(a_of_Ts[::-1]), np.log(Ts[::-1]))
gr_of_lnT = PchipInterpolator(np.log(Ts), gr_pl)
A_MIN, A_MAX = a_of_Ts[-1], a_of_Ts[0]


def T_of_a(a):
    a = np.asarray(a, float)
    return np.where(a > A_MAX, T0 / a, np.exp(lnT_of_lna(np.log(np.clip(a, A_MIN, A_MAX)))))


def rho_r(a):
    T = T_of_a(a)
    gr = np.where(np.asarray(a) > A_MAX, 2.0, gr_of_lnT(np.log(np.clip(T, Ts[0], Ts[-1]))))
    Tnu = (4 / 11)**(1 / 3) * T0 / np.asarray(a)
    return math.pi**2 / 30 * (gr * T**4 + 7 / 8 * 2 * Neff * Tnu**4)


def H_GeV(a):
    a = np.asarray(a, float)
    rho = rho_r(a) + Om * rho_crit / a**3 + (1 - Om) * rho_crit
    return np.sqrt(rho / (3 * Mpl**2))


def a_of_T(T):
    return float(T0 / T * (2 / (2 + g_e(me / T)[1]))**(1 / 3))


# ---- above 2 MeV: crude g_* table for the clock test only
def gstar_hi(T):
    if T <= 0.1: return 10.75
    if T <= 0.15: return 14.25
    if T <= 0.2: return 17.25
    if T <= 1.0: return 61.75
    if T <= 4.0: return 75.75
    return 86.25


def vrel_iso(v, w):
    vi, vj = v[:, None], v[None, :]
    big = np.maximum(vi, vj); small = np.minimum(vi, vj)
    rel = big + small**2 / (3 * np.where(big > 0, big, 1.0))
    W = w[:, None] * w[None, :]
    return float((W * rel).sum() / W.sum())


def main():
    out = {}
    T_i = 2e-3
    a_i = a_of_T(T_i)
    H_i = float(H_GeV(a_i))
    a_eq = brentq(lambda a: rho_r(a) - Om * rho_crit / a**3, 1e-5, 1e-2)
    out['background'] = {'a_2MeV': a_i, 'H_2MeV_GeV': H_i, 'H_over_mX': H_i / mX, 'a_eq': a_eq, 'z_eq': 1 / a_eq - 1,
                         'g_rho_e_2MeV': g_e(me / T_i)[0], 'g_s_e_2MeV': g_e(me / T_i)[1],
                         'H_2MeV_pure_g_10.75_formula': math.sqrt(math.pi**2 * 10.75 / 90) * T_i**2 / Mpl,
                         'codex_H_initial': 1.785634392341076e-24, 'reviewer_H_initial': 1.7859349018710206e-24,
                         'codex_z_eq': 3420.5, 'reviewer_z_eq': 3417.0}
    # ---- T4: collision rate
    s_i = 2 * math.pi**2 / 45 * 10.75 * T_i**3
    nX = f_charge * YQ * s_i / mX**3
    mean_v = math.sqrt(8 / (3 * math.pi)) * vrms
    GH_max = nX * sigma * math.sqrt(2) * mean_v * mX / H_i
    out['T4_collision_rate'] = {'n_X_over_mX3': nX, 'codex_n_X_from_ncore': 1.2141071422639381e-18 * Nplus,
                                'Gamma_over_H_maxwell': GH_max, 'codex': 15.72, 'codex_NR': 15.876924286306776, 'reviewer_maxwell': 16.00}
    npz = SP / 'codex' / 'late_kin_h0025_96.npz'
    w = None
    if npz.exists():
        z = np.load(npz); E = z['E']; n = len(E); Y = z['Y'][:, -1]
        w = np.maximum(Y[:n] + Y[n:2 * n], 0); p = np.sqrt(E * E - 1); v = p / E
        vr = vrel_iso(v, w)
        out['T4_collision_rate']['reproduction_codex_shape'] = {
            'N_plus': float(w.sum()), 'N_minus_raw': float(Y[n:2 * n].sum()), 'v_rms': float(math.sqrt(w @ (v * v) / w.sum())),
            'mean_v': float(w @ v / w.sum()), 'mean_vrel_iso_NR': vr,
            'Gamma_over_H': nX * sigma * vr * mX / H_i,
            'fraction_floor_bin': float(w[0] / w.sum()), 'p_rms_mX': float(math.sqrt(w @ (p * p) / w.sum()))}
    # ---- T5/T6: integrated collisions, viscosity crossing (frozen Maxwellian shape)
    x = np.geomspace(1, 10, 20001); a = a_i * x
    GH = GH_max * x**-4 * H_i / H_GeV(a)
    Ncoll = np.trapezoid(GH, np.log(x))
    xv = brentq(lambda xx: np.interp(xx, x, 2 / 3 * GH) - 1, 1, 10)
    xt = brentq(lambda xx: np.interp(xx, x, GH) - 1, 1, 10)
    out['T5_T6'] = {'collisions_x1_10': float(Ncoll), 'codex': 7.91173, 'reviewer': 7.901,
                    'viscosity_crossing_x': xv, 'viscosity_crossing_T_MeV': float(T_of_a(a_i * xv)) * 1e3, 'codex_T': 0.62156, 'reviewer_T': 0.6227,
                    'total_rate_crossing_x': xt, 'T_at_x10_MeV': float(T_of_a(10 * a_i)) * 1e3, 'codex_T_x10': 0.22086,
                    'pure_RD_collisions': GH_max * 0.5 * (1 - 0.01)}
    # ---- T7: free streaming
    def I_unit(a0, a1, N=30000):
        aa = np.geomspace(a0, a1, N)
        return float(np.trapezoid(a_i / (aa**2 * H_GeV(aa) / hbar), np.log(aa)) * c_kms / Mpc_km)   # Mpc per unit v at a_i
    I_eq = I_unit(a_i, a_eq); I_0 = I_unit(a_i, 1.0, 60000)
    RH = c_kms / (H_i / hbar) / Mpc_km / a_i
    fs = {'I_eq_pc_per_v': I_eq * 1e6, 'I_today_pc_per_v': I_0 * 1e6, 'comoving_Hubble_radius_2MeV_pc': RH * 1e6,
          'ln_aeq_ai': math.log(a_eq / a_i),
          'maxwell_mean_pc_eq': mean_v * I_eq * 1e6, 'maxwell_rms_pc_eq': vrms * I_eq * 1e6, 'maxwell_mean_pc_today': mean_v * I_0 * 1e6,
          'codex': {'mean': 35.51, 'rms': 39.86, 'today': 39.61}, 'reviewer_maxwell': {'mean': 36.58, 'rms': 39.70, 'today': 40.81}}
    if w is not None:
        aa = np.geomspace(a_i, a_eq, 30000)
        q = p[:, None] * a_i / aa[None, :]; vv = q / np.sqrt(1 + q * q)
        lam = np.trapezoid(vv / (aa * H_GeV(aa) / hbar)[None, :], np.log(aa), axis=1) * c_kms / Mpc_km
        fs['reproduction_codex_distribution'] = {'mean_pc': float(w @ lam / w.sum() * 1e6), 'rms_pc': float(math.sqrt(w @ lam**2 / w.sum()) * 1e6)}
    fs['from_x_1000_pc_per_v'] = I_unit(a_i * 1000, a_eq) * 1e6
    fs['maxwell_mean_from_x1000_pc'] = mean_v * I_unit(a_i * 1000, a_eq) * 1e6
    out['T7_free_streaming'] = fs
    # ---- T11 clock scaling: p_f fixed by the pulse, formation at T_f; lambda_fs to equality
    # a(T) above 2 MeV: entropy with g_s = g_* table, continuous at 2 MeV
    def a_hi(T):
        return a_i * (T_i / T) * (10.75 / gstar_hi(T))**(1 / 3)
    def H_hi(T):
        return math.sqrt(math.pi**2 * gstar_hi(T) / 90) * T**2 / Mpl
    p20 = vrms * 10.0                                    # p_rms at the 20 MeV clock, Codex's own redshift
    clock = {}
    for Tf in (2e-3, 2e-2, 2e-1, 2.0):
        af = a_hi(Tf) if Tf > T_i else a_i
        pf = p20                                          # same momentum at formation, set by the pulse energy
        # part 1: a_f -> a_i (RD, table), integrate in ln a with H(T(a)); T a = const within a g-plateau (approx)
        lam1 = 0.0
        if af < a_i:
            Ta = np.geomspace(Tf, T_i, 4000)[::-1] if False else np.geomspace(T_i, Tf, 4000)   # from T_i up to Tf
            aa = np.array([a_hi(T) for T in Ta]); HH = np.array([H_hi(T) for T in Ta])
            q = pf * af / aa; vv = q / np.sqrt(1 + q * q)
            # integrate over ln a (aa decreasing with T); sort ascending
            order = np.argsort(aa)
            lam1 = float(np.trapezoid((vv / (aa * HH / hbar))[order], np.log(aa[order])) * c_kms / Mpc_km)
        aa2 = np.geomspace(a_i, a_eq, 30000)
        q2 = pf * af / aa2; vv2 = q2 / np.sqrt(1 + q2 * q2)
        lam2 = float(np.trapezoid(vv2 / (aa2 * H_GeV(aa2) / hbar), np.log(aa2)) * c_kms / Mpc_km)
        v_at_2MeV = pf * af / a_i / math.sqrt(1 + (pf * af / a_i)**2)
        clock[f'T_f_{Tf*1e3:g}_MeV'] = {'a_f': af, 'p_rms_at_2MeV_mX': pf * af / a_i, 'v_rms_at_2MeV': v_at_2MeV,
                                       'lambda_rms_pc_to_eq': (lam1 + lam2) * 1e6, 'part_before_2MeV_pc': lam1 * 1e6,
                                       'reviewer_lambda_fs_mean_pc': {2e-3: 3.66, 2e-2: 36.6, 2e-1: 366., 2.0: 3658.}[Tf],
                                       'reviewer_v_rms_at_2MeV': {2e-3: 0.00684, 2e-2: 0.0684, 2e-1: 0.684, 2.0: 6.84}[Tf]}
    out['T11_clock_scaling'] = {'rows': clock,
                                'statement': 'p(2 MeV) = p_f a_f/a(2 MeV) ~ p_f (2 MeV/T_f): an EARLIER clock (larger T_f) gives a COLDER gas at 2 MeV and a SHORTER free-streaming length, lambda_fs ~ 1/T_f (times a log); the reviewer wrote v_rms(2 MeV) ~ T_f/20 MeV and lambda_fs ~ T_f, which is inverted.'}
    # pulse energy -> momentum (c481530 numbers)
    KE, N1000 = 2526.289, 11629.389
    Ebar = 1 + KE / N1000
    out['T11_velocity_origin'] = {'mean_E_minus_1_at_20MeV': Ebar - 1, 'p_if_monoenergetic': math.sqrt(Ebar**2 - 1), 'v_if_monoenergetic': math.sqrt(Ebar**2 - 1) / Ebar,
                                  'reviewer_NR_v_rms_0.66_is_p_over_m_not_v': True, 'p_over_10_vs_checkpoint_p_rms': math.sqrt(Ebar**2 - 1) / 10}
    # ---- T9: halo masses and the half-mode conventions
    rho_m = Omh2 * 2.77536627e11
    hm = {f'{k:g}': {'M_half_mode_Msun': 4 * math.pi / 3 * rho_m * (math.pi / k)**3, 'lambda_pc': 2 * math.pi / k * 1e6} for k in (1e4, 3e4, 5e4, 1e5, 3e5, 1e6, 40.)}
    def alpha_viel(m_keV):
        return 0.049 * m_keV**-1.11 * (Om / 0.25)**0.11 * (h / 0.7)**1.22 / h
    nu = 1.12
    k_T_half = lambda m: (2**(nu / 5) - 1)**(1 / (2 * nu)) / alpha_viel(m)      # T = 1/2 (Schneider 2012 half-mode)
    k_P_half = lambda m: (2**(nu / 10) - 1)**(1 / (2 * nu)) / alpha_viel(m)     # T^2 = 1/2
    Tnu0 = (4 / 11)**(1 / 3) * T0 * 1e9
    v0 = vrms * a_i
    m_eq = brentq(lambda m: 3.597 * Tnu0 * (0.12 * 94 / m)**(1 / 3) / m - v0, 1e3, 1e9)
    out['T9_scales'] = {'halo_masses': hm, 'thermal_equivalent_mass_eV': m_eq,
                        'k_half_mode_T_eq_half_Mpc^-1': k_T_half(m_eq / 1e3), 'k_half_mode_T2_eq_half_Mpc^-1': k_P_half(m_eq / 1e3),
                        'reviewer_k_half_27627': 27627.76, 'reviewer_label': "'T^2 = 1/2' in transport_check.json, but the code uses 2**(1.12/5), i.e. T = 1/2",
                        'k_half_5.3keV_T_half': k_T_half(5.3), 'k_half_5.3keV_T2_half': k_P_half(5.3),
                        'BGK_half_suppression_target': (1 + (1 - fgas)**2) / 2,
                        'warm_T_at_BGK_half_suppression_ignoring_backreaction': (math.sqrt((1 + (1 - fgas)**2) / 2) - (1 - fgas)) / fgas}
    # ---- T10: Lyman-alpha toy bound
    lam2 = (vrms * I_eq)**2
    out['T10_lya_toy'] = {f'{kh:g} h/Mpc': fgas * (kh * h)**2 * lam2 / 3 for kh in (10., 30.)}
    RES.mkdir(parents=True, exist_ok=True)
    (RES / 'sk_rates.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
