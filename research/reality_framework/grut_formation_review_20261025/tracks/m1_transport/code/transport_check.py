"""Track m1_transport, items 2-3: collisions, decoupling, free streaming, halo masses.

Independent code: own standard background (photons, e+- with exact Fermi-Dirac integrals and
entropy conservation, N_eff = 3.044 neutrinos with T_nu = (4/11)^(1/3) T0/a, Omega_m h^2 = 0.143,
flat Lambda), own collision integrals and free-streaming quadrature.  Nothing from Codex's code is
imported.  Two inputs are taken from Codex's packets and labelled as such when used:
  * the checkpoint gas distribution late_kin_h0025_96.npz (6c504de) -- used for a reproduction of
    Codex's non-Maxwellian numbers next to the Maxwellian independent numbers;
  * the model constants Y_Q, m_X and the gas charge fraction (sigma is re-derived in xsec.py).

Usage: python3 -I transport_check.py  (reads the npz from $SP/review/m1_transport/codex/ if present).
"""
import json, math, os
from pathlib import Path
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.stats import chi

ROOT = Path(__file__).resolve().parents[1]
SP = Path(os.environ.get('SP', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')) / 'review' / 'm1_transport'

# ---------------------------------------------------------------- constants
kB_eV_per_K = 8.617333262e-5
T0_K = 2.7255
T0_GeV = T0_K * kB_eV_per_K * 1e-9
me_GeV = 0.51099895e-3
Mpl_red = 2.435e18                      # GeV, as in Codex's clock
GeV_to_s_inv = 1 / 6.582119569e-25      # 1 GeV = 1.519e24 s^-1
c_km_s = 299792.458
Mpc_km = 3.0856775814913673e19
Neff = 3.044
Omh2 = 0.143
h = 0.674
Om = Omh2 / h**2
H0_GeV = 100 * h / Mpc_km * 6.582119569e-25   # km/s/Mpc -> s^-1 -> GeV
rho_crit = 3 * H0_GeV**2 * Mpl_red**2          # GeV^4
rho_gamma0 = math.pi**2 / 15 * T0_GeV**4
Og = rho_gamma0 / rho_crit
Onu = Neff * 7 / 8 * (4 / 11)**(4 / 3) * Og
Or = Og + Onu
OL = 1 - Om - Or
a_eq = Or / Om

# model constants (Codex CANDIDATE_M1.json / 6c504de; labelled [Q])
mX_GeV = 0.09240896373753621
YQ = 5.079414232725842e-9
Qref = 200000.0
N_plus = 13189.13146074248
f_charge = N_plus / Qref                   # 6.5946 % [Q 6c504de]
sigma_mX2 = (4e-3 - 0.04 / 3)**2 / (128 * math.pi)   # re-derived in xsec.py
T_i = 2e-3                                 # GeV, checkpoint clock
v_rms = 0.0684
fgas = 0.058767411685782034                # Codex test-model gas fraction of all matter [Q]


# ---------------------------------------------------------------- e+- thermodynamics
def fd_integrals(y):
    """(g_rho_e, g_s_e) for e+e- (4 dof) at m_e/T = y; massless limit (3.5, 3.5)."""
    rho = quad(lambda x: x * x * math.sqrt(x * x + y * y) / (math.exp(math.sqrt(x * x + y * y)) + 1), 0, 60, epsabs=0, epsrel=1e-11)[0]
    P = quad(lambda x: x**4 / (3 * math.sqrt(x * x + y * y)) / (math.exp(math.sqrt(x * x + y * y)) + 1), 0, 60, epsabs=0, epsrel=1e-11)[0]
    rho *= 2 / math.pi**2
    P *= 2 / math.pi**2
    return 30 / math.pi**2 * rho, 45 / (2 * math.pi**2) * (rho + P)


class Background:
    """Photon temperature from entropy conservation of the gamma e+- plasma (neutrinos decoupled)."""

    def __init__(self):
        Ts = np.geomspace(1e-6, 0.05, 700)     # GeV
        gr = np.empty_like(Ts); gs = np.empty_like(Ts)
        for i, T in enumerate(Ts):
            gr[i], gs[i] = fd_integrals(me_GeV / T)
        self.gs_tot = 2 + gs; self.gr_tot = 2 + gr
        a_of_T = T0_GeV / Ts * (2 / self.gs_tot)**(1 / 3)
        self.lnT_of_lna = PchipInterpolator(np.log(a_of_T[::-1]), np.log(Ts[::-1]))
        self.gr_of_lnT = PchipInterpolator(np.log(Ts), self.gr_tot)
        self.gs_of_lnT = PchipInterpolator(np.log(Ts), self.gs_tot)
        self.a_min = a_of_T[-1]; self.a_max = a_of_T[0]
        self.a_i = float(self.a_of_T(T_i))

    def a_of_T(self, T):
        gs = 2 + fd_integrals(me_GeV / T)[1]
        return T0_GeV / T * (2 / gs)**(1 / 3)

    def T(self, a):
        a = np.asarray(a, float)
        late = T0_GeV / np.maximum(a, 1e-300)
        early = np.exp(self.lnT_of_lna(np.log(np.clip(a, self.a_min, self.a_max))))
        return np.where(a > self.a_max, late, early)

    def E2(self, a):
        """H^2/H0^2."""
        a = np.asarray(a, float); T = self.T(a)
        gr = np.where(a <= self.a_max, self.gr_of_lnT(np.log(np.clip(T, 1e-6, 0.05))), 2.0)
        rho_ge = math.pi**2 / 30 * gr * T**4
        rho_nu = Onu * rho_crit / a**4
        return (rho_ge + rho_nu) / rho_crit + Om / a**3 + OL

    def H_GeV(self, a):
        return H0_GeV * np.sqrt(self.E2(a))

    def H_s(self, a):
        return self.H_GeV(a) * GeV_to_s_inv


def vrel_isotropic(v, w):
    """<|v1 - v2|> for isotropic speeds v with weights w (NR): v> + v<^2/(3 v>)."""
    vi = v[:, None]; vj = v[None, :]
    big = np.maximum(vi, vj); small = np.minimum(vi, vj)
    rel = big + small**2 / (3 * np.where(big > 0, big, 1))
    W = w[:, None] * w[None, :]
    return float((W * rel).sum() / W.sum())


def main():
    out = {'conventions': 'independent code; see module docstring'}
    bg = Background()
    a_i = bg.a_i
    H_i = float(bg.H_GeV(a_i))
    out['background'] = {'T0_GeV': T0_GeV, 'Omega_gamma_h2': Og * h**2, 'Omega_r_h2': Or * h**2, 'Omega_m_h2': Omh2,
                         'h_assumed': h, 'a_equality': a_eq, 'z_equality': 1 / a_eq - 1,
                         'a_at_2MeV': a_i, 'H_at_2MeV_GeV': H_i, 'H_over_mX_at_2MeV': H_i / mX_GeV,
                         'g_s_gamma_e_at_2MeV': float(bg.gs_of_lnT(math.log(T_i))),
                         'codex_a_initial': 8.394521281378524e-11, 'codex_H_initial_GeV': 1.785634392341076e-24,
                         'codex_z_eq': 3420.5}
    # ------------------------------------------------------------ item 2: collision rate at 2 MeV
    s_i = 2 * math.pi**2 / 45 * 10.75 * T_i**3           # GeV^3, standard g_*s = 10.75
    n_X = f_charge * YQ * s_i                            # GeV^3
    n_X_mX3 = n_X / mX_GeV**3
    mean_v = math.sqrt(8 / (3 * math.pi)) * v_rms
    mean_vrel = math.sqrt(2) * mean_v
    Gamma_maxwell = n_X_mX3 * sigma_mX2 * mean_vrel * mX_GeV   # GeV
    out['collisions_2MeV'] = {
        'entropy_density_GeV3': s_i, 'n_X_GeV3': n_X, 'n_X_over_mX3': n_X_mX3,
        'codex_n_core_over_mX3_at_20MeV': 1.2141071422639381e-15,
        'implied_n_X_over_mX3_codex': 1.2141071422639381e-15 / 1000 * N_plus,
        'sigma_mX2': sigma_mX2, 'mean_v_maxwell': mean_v, 'mean_vrel_maxwell': mean_vrel,
        'Gamma_per_particle_maxwell_GeV': Gamma_maxwell, 'Gamma_over_H_maxwell': Gamma_maxwell / H_i,
        'codex_Gamma_over_H': {'96bin_relativistic_kernel': 15.72, 'NR_threshold_elastic_py': 15.876924286306776}}
    npz = SP / 'codex' / 'late_kin_h0025_96.npz'
    w = None
    if npz.exists():
        z = np.load(npz); E = z['E']; n = len(E); Y = z['Y'][:, -1]
        w = np.maximum(Y[:n] + Y[n:2 * n], 0)
        p = np.sqrt(E * E - 1); v = p / E
        vr_chk = math.sqrt(float(w @ (v * v) / w.sum()))
        vrel_chk = vrel_isotropic(v, w)
        Gamma_chk = n_X_mX3 * sigma_mX2 * vrel_chk * mX_GeV
        out['collisions_2MeV']['checkpoint_reproduction'] = {
            'source': 'late_kin_h0025_96.npz (Codex 6c504de), weights = plus+minus clipped at 0; own integrals',
            'N_plus': float(w.sum()), 'v_rms': vr_chk, 'mean_v': float(w @ v / w.sum()), 'mean_vrel_isotropic_NR': vrel_chk,
            'kurtosis_p4_over_p2sq': float((w @ p**4) * w.sum() / (w @ p**2)**2),
            'Gamma_over_H': Gamma_chk / H_i,
            'ratio_vrel_checkpoint_over_maxwell_same_vrms': vrel_chk / (math.sqrt(2) * math.sqrt(8 / (3 * math.pi)) * vr_chk)}
    # ------------------------------------------------------------ integrated collisions, viscosity crossing
    x = np.geomspace(1, 10, 4001); a = a_i * x
    Hx = bg.H_GeV(a)
    shapes = [('maxwell_shape', Gamma_maxwell)]
    if w is not None:
        shapes.append(('checkpoint_shape_frozen', Gamma_chk))
    for label, G0 in shapes:
        GH = G0 * x**-4 / Hx
        Ncoll = cumulative_trapezoid(GH, np.log(x), initial=0)
        visc = 2 / 3 * GH
        xdec = float(brentq(lambda xx: np.interp(xx, x, visc) - 1, 1, 10)) if visc[-1] < 1 < visc[0] else None
        xcoll1 = float(brentq(lambda xx: np.interp(xx, x, GH) - 1, 1, 10)) if GH[-1] < 1 < GH[0] else None
        out['integrated_' + label] = {
            'Gamma_over_H_initial': float(GH[0]), 'Gamma_over_H_at_x10': float(GH[-1]),
            'collisions_per_particle_x1_to_10': float(Ncoll[-1]), 'collisions_x1_to_3.262': float(np.interp(3.2619, x, Ncoll)),
            'viscosity_proxy_crossing_x': xdec, 'viscosity_proxy_crossing_T_MeV': float(bg.T(a_i * xdec)) * 1e3 if xdec else None,
            'total_rate_crossing_x': xcoll1, 'total_rate_crossing_T_MeV': float(bg.T(a_i * xcoll1)) * 1e3 if xcoll1 else None,
            'T_at_x10_MeV': float(bg.T(a_i * 10)) * 1e3,
            'pure_RD_estimate_collisions': float(GH[0]) * 0.5 * (1 - 1 / 100),
            'codex': {'collisions_x1_to_10': 7.91173, 'viscosity_crossing_x': 3.26194, 'viscosity_crossing_T_MeV': 0.62156, 'T_at_x10_MeV': 0.22086}}
    # ------------------------------------------------------------ item 3: free streaming
    def I_unit(a_start, a_stop, N=20000):
        """comoving distance (Mpc) per unit velocity at a_i for v = v_i a_i/a:  c a_i int dlna /(a^2 H)."""
        aa = np.geomspace(a_start, a_stop, N)
        integrand = a_i / (aa**2 * bg.H_s(aa))
        return float(np.trapezoid(integrand, np.log(aa)) * c_km_s / Mpc_km)

    def lambda_exact(p_i, a_start, a_stop, N=20000):
        """exact relativistic v(a) = q/sqrt(1+q^2), q = p_i a_i/a (p_i in m_X units); comoving Mpc."""
        aa = np.geomspace(a_start, a_stop, N)
        q = p_i[:, None] * a_i / aa[None, :]
        vv = q / np.sqrt(1 + q * q)
        integ = vv / (aa * bg.H_s(aa))[None, :]
        return np.trapezoid(integ, np.log(aa), axis=1) * c_km_s / Mpc_km

    I_eq = I_unit(a_i, a_eq)
    I_today = I_unit(a_i, 1.0, 40000)
    I_from_dec = I_unit(a_i * 3.2619, a_eq)
    I_from_x10 = I_unit(a_i * 10, a_eq)
    I_20_to_2MeV = I_unit(a_i / 10, a_i)
    x_dec_core = 36.0
    I_from_core_dec = I_unit(a_i * x_dec_core, a_eq)
    scale = v_rms / math.sqrt(3)
    fs = {'I_eq_pc_per_unit_v': I_eq * 1e6, 'I_today_pc_per_unit_v': I_today * 1e6,
          'I_from_0.62MeV_pc_per_unit_v': I_from_dec * 1e6, 'I_from_0.2MeV_pc_per_unit_v': I_from_x10 * 1e6,
          'I_from_x36_pc_per_unit_v': I_from_core_dec * 1e6,
          'I_20MeV_to_2MeV_pc_per_unit_v_at_2MeV': I_20_to_2MeV * 1e6,
          'comoving_Hubble_radius_at_2MeV_pc': c_km_s / bg.H_s(a_i) / Mpc_km / a_i * 1e6,
          'ln_aeq_over_ai': math.log(a_eq / a_i),
          'codex_I_eq_from_rms_pc': 39.86 / 0.06837528835889453}
    fs['maxwellian_v_rms_0.0684'] = {
        'mean_pc_to_eq': mean_v * I_eq * 1e6, 'rms_pc_to_eq': v_rms * I_eq * 1e6,
        'p10_pc_to_eq': float(chi.ppf(.1, 3) * scale * I_eq * 1e6), 'p50_pc_to_eq': float(chi.ppf(.5, 3) * scale * I_eq * 1e6),
        'p90_pc_to_eq': float(chi.ppf(.9, 3) * scale * I_eq * 1e6),
        'mean_pc_to_today': mean_v * I_today * 1e6,
        'mean_pc_from_0.62MeV_to_eq': mean_v * I_from_dec * 1e6, 'rms_pc_from_0.62MeV_to_eq': v_rms * I_from_dec * 1e6,
        'mean_pc_from_0.2MeV_to_eq': mean_v * I_from_x10 * 1e6,
        'mean_pc_if_dragged_by_cores_until_x36': mean_v * I_from_core_dec * 1e6,
        'codex': {'mean_pc': 35.51, 'rms_pc': 39.86, 'median_pc': 34.06, 'p10_pc': 12.12, 'p90_pc': 61.36, 'mean_today_pc': 39.61,
                  'delayed_x2': 33.69, 'delayed_x4': 31.89, 'delayed_x10': 29.56}}
    if w is not None:
        lam = lambda_exact(p, a_i, a_eq)
        order = np.argsort(lam); cw = (np.cumsum(w[order]) - .5 * w[order]) / w.sum()
        fs['checkpoint_reproduction_exact_relativistic'] = {
            'mean_pc': float(w @ lam / w.sum() * 1e6), 'rms_pc': float(math.sqrt(w @ lam**2 / w.sum()) * 1e6),
            'p10_pc': float(np.interp(.1, cw, lam[order]) * 1e6), 'p50_pc': float(np.interp(.5, cw, lam[order]) * 1e6),
            'p90_pc': float(np.interp(.9, cw, lam[order]) * 1e6),
            'mean_pc_NR_formula': float(w @ v / w.sum() * I_eq * 1e6),
            'mean_pc_to_today': float(w @ lambda_exact(p, a_i, 1.0, 40000) / w.sum() * 1e6)}
    out['free_streaming'] = fs
    lam_mean = mean_v * I_eq
    aH_i = a_i * bg.H_s(a_i) * Mpc_km / c_km_s
    out['free_streaming_scales'] = {
        'lambda_fs_mean_Mpc': lam_mean, 'k_1_over_lambda_Mpc^-1': 1 / lam_mean, 'k_2pi_over_lambda_Mpc^-1': 2 * math.pi / lam_mean,
        'comoving_aH_at_2MeV_Mpc^-1': aH_i, 'k_fs_at_2MeV_aH_over_vrms_Mpc^-1': aH_i / v_rms,
        'codex_BGK_P_mix_over_P_cold': {'1e4': 0.997879, '3e4': 0.979835, '1e5': 0.893210, '3e5': 0.849806, '1e6': 0.764169},
        'implied_gas_transfer_T_gas': {}}
    for k, r in out['free_streaming_scales']['codex_BGK_P_mix_over_P_cold'].items():
        out['free_streaming_scales']['implied_gas_transfer_T_gas'][k] = (math.sqrt(r) - (1 - fgas)) / fgas
    out['free_streaming_scales']['note_T_gas'] = 'T_gas = (sqrt(P_mix/P_cold)-(1-f))/f ignores the slowed growth of the cold part; it goes negative where that effect dominates (k >= 3e5).'
    # ------------------------------------------------------------ halo masses
    rho_m = Omh2 * 2.77536627e11       # Msun / Mpc^3 comoving (h-free)
    hm = {}
    for k in (1e4, 3e4, 1e5, 3e5, 1e6, 1e7, 20.0, 27.0, 40.0):
        hm[f'{k:g}'] = {'lambda_pc': 2 * math.pi / k * 1e6,
                        'M_half_wavelength_sphere_Msun': 4 * math.pi / 3 * rho_m * (math.pi / k)**3,
                        'M_full_wavelength_sphere_Msun': 4 * math.pi / 3 * rho_m * (2 * math.pi / k)**3}
    out['halo_masses'] = {'rho_m_Msun_Mpc3': rho_m, 'table': hm,
                          'definition': 'M(k) = (4 pi/3) rho_m (pi/k)^3 (half-mode convention, Schneider et al 2012 arXiv:1112.0330); full-wavelength sphere is 8x larger'}
    # ------------------------------------------------------------ equivalent thermal WDM
    Tnu0_eV = (4 / 11)**(1 / 3) * T0_GeV * 1e9
    v0 = v_rms * a_i

    def vrms_thermal(m_eV, Omh2w=0.12):
        Tw = Tnu0_eV * (Omh2w * 94 / m_eV)**(1 / 3)
        return 3.597 * Tw / m_eV
    m_eq = brentq(lambda m: vrms_thermal(m) - v0, 1e2, 1e9)

    def alpha_viel(m_keV):
        return 0.049 * m_keV**-1.11 * (Om / 0.25)**0.11 * (h / 0.7)**1.22 / h   # Mpc

    def k_half(m_keV):
        return (2**(1.12 / 5) - 1)**(1 / 2.24) / alpha_viel(m_keV)
    out['equivalent_thermal_WDM'] = {
        'gas_v_rms_today_c': v0, 'thermal_relic_mass_with_same_v_rms_today_eV': m_eq,
        'k_half_mode_of_that_relic_Mpc^-1': k_half(m_eq / 1e3),
        'k_half_mode_Lyman_alpha_bound_5.3keV_Mpc^-1': k_half(5.3), 'k_half_mode_3.5keV_Mpc^-1': k_half(3.5),
        'gas_half_mode_from_codex_BGK_Mpc^-1_approx': 5e4,
        'thermal_mass_with_k_half_5e4_keV': brentq(lambda m: k_half(m) - 5e4, 1, 1e6),
        'formula': 'Viel et al 2005 (astro-ph/0501562) alpha, nu = 1.12, T^2 = 1/2; v_rms = 3.597 T_w/m, T_w from Omega h^2 = (m/94 eV)(T_w/T_nu)^3'}
    # ------------------------------------------------------------ Lyman-alpha weighted bound
    for kh in (10., 20., 30.):
        k = kh * h
        out.setdefault('lyman_alpha_toy', {})[f'{kh:g} h/Mpc'] = {
            'k_Mpc^-1': k, 'no_gravity_bound_abs_dP_over_P': fgas * k**2 * (v_rms * I_eq)**2 / 3,
            'two_f_times_k_lambda_squared': 2 * fgas * (k * lam_mean)**2}
    # ------------------------------------------------------------ gas-core elastic scattering (E)
    n_core_2MeV = 1.2141071422639381e-15 / 1000
    out['gas_core_elastic_estimate'] = {}
    for sig in (100., 300., 1000., 2560.):
        GH = n_core_2MeV * sig * v_rms * mX_GeV / H_i
        out['gas_core_elastic_estimate'][f'sigma_{sig:g}_mX2'] = {
            'Gamma_over_H_2MeV': GH, 'x_decouple_sqrt': math.sqrt(GH), 'T_decouple_MeV': float(bg.T(a_i * math.sqrt(GH))) * 1e3 if GH > 1 else None,
            'fractional_energy_transfer_to_cores_per_Hubble': GH * mX_GeV / 17324.83235}
    out['gas_core_elastic_estimate']['note'] = ('s-wave unitarity bound 4 pi/k^2 at k = m_X v_rms = 0.0684 m_X is 2686/m_X^2; geometric 4 pi R^2 with R ~ 10/m_X is ~1260/m_X^2. '
                                                'Codex bf8ab3b found 99% elastic reflection of l=0 Xbar packets at k=0.4-1; X has the same diagonal potential S^2+2aF^2. Not computed here: the actual X-core phase shifts at k ~ 0.07.')
    (ROOT / 'results' / 'transport_check.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
