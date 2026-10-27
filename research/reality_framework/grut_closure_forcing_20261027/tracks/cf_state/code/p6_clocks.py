"""P6 (no external time): relational clocks, realistic-clock decoherence, and comparisons (own code).

Part A [C, D]: Page-Wootters conditional state of a qubit read against a free-particle clock (hbar = 1 units).
  For a non-interacting clock + system, the conditional state at clock reading x_T is
      rho_S(T) = int dt w(t) U_S(t) rho_S U_S(t)^dag / int dt w(t),   w(t) = |psi_C(x_T, t)|^2,
  so |rho_01(T)|/|rho_01(0)| = |int dt w(t) e^{-i omega t}| / int dt w(t).
  (i) ideal clock (linear dispersion, no spreading): coherence independent of T, -> 1 as the packet width -> 0;
  (ii) free-particle clock: packet spreading gives coherence loss that grows with T; at the optimal initial width
      s0^2 = T/(2M) the timing variance is dT^2 = T/(M v^2) (Salecker-Wigner form).
  The integral is evaluated (a) with the analytic free Gaussian packet and (b) with an FFT split-step evolution of the
  packet on a grid, as an independent check.
Part B [C]: predicted decoherence exponents D (|rho_12(T)/rho_12(0)| = exp(-D)) for the realistic-clock models and
  for standard time-dilation effects, for an optical clock, a 25 kDa molecule interferometer, a proposed 1e-14 kg
  superposition, a half-metre atom interferometer and the GRUT release epoch; required energy splittings.
Usage: python3 -B p6_clocks.py ../results/p6_clocks.json
"""
import sys
import numpy as np
from scipy.integrate import quad
from common import hbar, c, G, kB, u, g_earth, t_P, GeV, QUOTED, dump


# ---------------- Part A: Page-Wootters conditional coherence ----------------
def coherence_analytic(omega, M, s0, v, T, ideal=False):
    x = v * T

    def w(t):
        s2 = s0**2 if ideal else s0**2 + (t / (2 * M * s0))**2
        return np.exp(-(x - v * t)**2 / (2 * s2)) / np.sqrt(2 * np.pi * s2)
    # integrate over a window of +-12 widths around T
    s2T = s0**2 if ideal else s0**2 + (T / (2 * M * s0))**2
    half = 12 * np.sqrt(s2T) / v + 1e-12
    a, b = T - half, T + half
    norm = quad(w, a, b, limit=400, epsabs=0, epsrel=1e-12)[0]
    re = quad(lambda t: w(t) * np.cos(omega * t), a, b, limit=800, epsabs=0, epsrel=1e-11)[0]
    im = quad(lambda t: -w(t) * np.sin(omega * t), a, b, limit=800, epsabs=0, epsrel=1e-11)[0]
    return abs(re + 1j * im) / norm


def coherence_grid(omega, M, s0, v, T, N=2**16):
    """Independent check: evolve the clock packet with an exact FFT free propagator on a grid and read
    w(t) = |psi(x_T, t)|^2 at x_T = v T, for t on a grid around T; then the same ratio of integrals."""
    p0 = M * v
    sT = np.sqrt(s0**2 + (T / (2 * M * s0))**2)
    L = v * T + 40 * sT + 40 * s0
    xs = np.linspace(-20 * s0 - 20 * sT, L, N, endpoint=False)
    dx = xs[1] - xs[0]
    k = 2 * np.pi * np.fft.fftfreq(N, dx)
    psi0 = (2 * np.pi * s0**2)**-0.25 * np.exp(-xs**2 / (4 * s0**2) + 1j * p0 * xs)
    phi0 = np.fft.fft(psi0)
    half = 10 * sT / v
    ts = np.linspace(T - half, T + half, 1201)
    xT = v * T
    ix = int(round((xT - xs[0]) / dx))
    # evaluate psi(x_T, t) by direct Fourier sum at the single point x_T (exact free evolution)
    phase_x = np.exp(1j * k * (xs[ix] - xs[0]))
    wv = np.empty_like(ts)
    for j, t in enumerate(ts):
        val = np.sum(phi0 * np.exp(-1j * k**2 * t / (2 * M)) * phase_x) / N
        wv[j] = abs(val)**2
    # shift: the grid point x_ix differs from x_T by < dx/2; use t-grid trapezoid
    norm = np.trapezoid(wv, ts)
    z = np.trapezoid(wv * np.exp(-1j * omega * ts), ts)
    return abs(z) / norm, float(xs[ix] - xT)


def part_A():
    out = {'units': 'hbar = 1; clock mass M, speed v, initial width s0; system qubit splitting omega', 'rows': []}
    M, v = 50.0, 1.0
    for omega in (0.5, 1.0, 2.0):
        for T in (10.0, 100.0, 1000.0):
            s0opt = np.sqrt(T / (2 * M))
            dT2_opt = T / (M * v**2)            # Salecker-Wigner form at optimal s0
            free = coherence_analytic(omega, M, s0opt, v, T)
            ideal = coherence_analytic(omega, M, s0opt, v, T, ideal=True)
            ideal_narrow = coherence_analytic(omega, M, 1e-3, v, T, ideal=True)
            gauss_pred = np.exp(-omega**2 * dT2_opt / 2)
            row = dict(omega=omega, T=T, s0_opt=s0opt, dT2_SW=dT2_opt,
                       free_clock=free, gaussian_prediction=gauss_pred,
                       rel_diff=(free - gauss_pred) / gauss_pred if gauss_pred > 1e-6 else None,
                       ideal_clock_same_width=ideal, ideal_clock_width_1em3=ideal_narrow)
            out['rows'].append(row)
    # grid check on two points
    checks = []
    for omega, T in ((1.0, 100.0), (2.0, 100.0), (0.5, 1000.0)):
        s0opt = np.sqrt(T / (2 * M))
        g, off = coherence_grid(omega, M, s0opt, v, T)
        a = coherence_analytic(omega, M, s0opt, v, T)
        checks.append(dict(omega=omega, T=T, grid=g, analytic=a, abs_diff=abs(g - a), xgrid_minus_xT=off))
    out['grid_checks'] = checks
    # T-independence of the ideal clock and T-growth of the free clock at fixed s0
    s0 = 1.0
    out['fixed_s0_scan'] = [dict(T=T, free=coherence_analytic(1.0, M, s0, v, T),
                                 ideal=coherence_analytic(1.0, M, s0, v, T, ideal=True))
                            for T in (1.0, 10.0, 100.0, 1000.0)]
    return out


# ---------------- Part B: model predictions for real systems ----------------
def D_GPP_opt(omega, T):          # hep-th/0406260 eq. (6): ln|rho12(T)/rho12(0)| = -(3/2) t_P^{4/3} T^{2/3} w^2
    return 1.5 * t_P**(4 / 3) * T**(2 / 3) * omega**2


def D_GPP_NJP(omega, T):          # gr-qc/0402118 eq. (55): (w1 - w2)^2 t_P t_life (sigma ~ t_P)
    return omega**2 * t_P * T


def D_Milburn(omega, T, tau=t_P):  # Milburn PRA 44, 5401 (1991), second-order form: rate w^2 tau / 2  [R]
    return 0.5 * omega**2 * tau * T


def D_SW(omega, T, M):            # free-particle (light-speed) clock of mass M: dT^2 = hbar T/(M c^2), D = w^2 dT^2/2
    return 0.5 * omega**2 * hbar * T / (M * c**2)


def tau_Pikovski(N, Tint, dx):    # arXiv:1311.1095 eq. (3): sqrt(2/N) hbar c^2/(kB T g dx), Gaussian decay exp(-(t/tau)^2)
    return np.sqrt(2.0 / N) * hbar * c**2 / (kB * Tint * g_earth * dx)


def t_CastroRuiz(N, dE, x):       # arXiv:1507.01955 eq. (8): 2 hbar c^4 x/(sqrt(N) G dE^2)
    return 2 * hbar * c**4 * x / (np.sqrt(N) * G * dE**2)


def required_omega(model, T, eps):
    if model == 'GPP_opt':
        return np.sqrt(eps / (1.5 * t_P**(4 / 3) * T**(2 / 3)))
    if model == 'GPP_NJP':
        return np.sqrt(eps / (t_P * T))
    if model == 'Milburn_tP':
        return np.sqrt(2 * eps / (t_P * T))


def part_B():
    res = {}
    eV = 1.602176634e-19
    # S1: Sr optical clock (Bothwell et al. 2022)
    f = QUOTED['Sr_clock_freq'][0]
    om = 2 * np.pi * f
    T = QUOTED['Bothwell_coherence'][0]
    D_bound = T / QUOTED['Bothwell_extra_decoherence_time'][0]
    m_Sr = 86.909 * u
    s1 = dict(system='87Sr optical clock superposition (Bothwell et al. 2022)', omega=om, T=T,
              D_bound_all_extra_decoherence=D_bound,
              D_GPP_opt=D_GPP_opt(om, T), D_GPP_NJP=D_GPP_NJP(om, T), D_Milburn_tP=D_Milburn(om, T),
              D_SW_clock_1kg=D_SW(om, T, 1.0),
              t_CastroRuiz_two_atoms_1mm=t_CastroRuiz(2, hbar * om, 1e-3),
              t_CastroRuiz_1e4_atoms_1mm=t_CastroRuiz(1e4, hbar * om, 1e-3))
    for k in ('D_GPP_opt', 'D_GPP_NJP', 'D_Milburn_tP', 'D_SW_clock_1kg'):
        s1['gap_' + k] = D_bound / s1[k]
    # bound on a Milburn-type tau from the Sr data: D <= D_bound
    s1['Milburn_tau_upper_bound_s'] = 2 * D_bound / (om**2 * T)
    s1['Milburn_tau_upper_bound_over_tP'] = s1['Milburn_tau_upper_bound_s'] / t_P
    # clock-comparison timing floor (Ng-van Dam / GPP optimal clock)
    for Tav in (1e4, 3.3e5):
        dT_NvD = np.sqrt(1.5) * t_P**(2 / 3) * Tav**(1 / 3)
        s1[f'NvD_fractional_timing_floor_T{Tav:.0e}'] = dT_NvD / Tav
    s1['measured_fractional_uncertainty'] = QUOTED['Bothwell_frac_uncertainty'][0]
    s1['ratio_measured_to_NvD_floor_T1e4'] = s1['measured_fractional_uncertainty'] / s1['NvD_fractional_timing_floor_T1e+04']
    res['S1_optical_clock'] = s1

    # S2: Fein et al. 2019 molecule interferometer
    m = QUOTED['Fein_mass_u'][0] * u
    v = 2 * np.pi * hbar / (m * QUOTED['Fein_lambda_dB'][0])
    T = QUOTED['Fein_length'][0] / v
    d = QUOTED['Fein_grating_period'][0]
    dE = m * g_earth * d            # worst case: separation taken vertical [E]
    om = dE / hbar
    D_bound = -np.log(QUOTED['Fein_visibility_fraction'][0])
    N_modes = 3 * 2000
    s2 = dict(system='25 kDa oligoporphyrin, Talbot-Lau (Fein et al. 2019)', mass_kg=m, v=v, T=T, dE_J=dE, omega=om,
              D_bound_from_visibility=D_bound,
              D_GPP_opt=D_GPP_opt(om, T), D_GPP_NJP=D_GPP_NJP(om, T), D_Milburn_tP=D_Milburn(om, T),
              tau_Pikovski_500K=tau_Pikovski(N_modes, 500.0, d),
              D_Pikovski_500K=(T / tau_Pikovski(N_modes, 500.0, d))**2)
    for k in ('D_GPP_opt', 'D_GPP_NJP', 'D_Milburn_tP', 'D_Pikovski_500K'):
        s2['gap_' + k] = D_bound / s2[k]
    res['S2_molecule'] = s2

    # S3: proposed macroscopic superposition (Bose et al. 2017 parameters), hypothetical
    m = QUOTED['Bose_mass'][0]
    dx = QUOTED['Bose_dx'][0]
    T = QUOTED['Bose_tau'][0]
    dE = m * g_earth * dx           # worst case vertical [E]; Bose et al. separate horizontally
    om = dE / hbar
    N_modes = 3 * m / (20 * u)
    s3 = dict(system='proposed 1e-14 kg, 250 um, 2.5 s (Bose et al. 2017), not yet realized', dE_J=dE, omega=om, T=T,
              D_GPP_opt=D_GPP_opt(om, T), D_GPP_NJP=D_GPP_NJP(om, T), D_Milburn_tP=D_Milburn(om, T),
              tau_Pikovski_1K=tau_Pikovski(N_modes, 1.0, dx), tau_Pikovski_300K=tau_Pikovski(N_modes, 300.0, dx),
              D_Pikovski_1K=(T / tau_Pikovski(N_modes, 1.0, dx))**2,
              D_Pikovski_300K=(T / tau_Pikovski(N_modes, 300.0, dx))**2)
    res['S3_macroscopic_proposed'] = s3

    # S4: half-metre atom interferometer (Kovachy et al. Nature 528, 530 (2015)) [R]: 87Rb, 0.54 m, ~2 s
    m = 86.909 * u
    dE = m * g_earth * 0.54
    om = dE / hbar
    T = 2.0
    res['S4_atom_interferometer_R'] = dict(system='87Rb, 0.54 m, ~2 s (Kovachy et al. 2015) [R]', omega=om, T=T,
                                          D_GPP_opt=D_GPP_opt(om, T), D_GPP_NJP=D_GPP_NJP(om, T),
                                          D_Milburn_tP=D_Milburn(om, T))

    # S5: GRUT release epoch, heavy spot point H_i = 1.3335e11 GeV: omega ~ H/hbar, T ~ 1/H
    Hi = 1.3335e11 * GeV / hbar
    res['S5_GRUT_release'] = dict(system='GRUT heavy point release: omega ~ H_i, T ~ 1/H_i', omega=Hi, T=1 / Hi,
                                  D_GPP_opt=D_GPP_opt(Hi, 1 / Hi), D_GPP_NJP=D_GPP_NJP(Hi, 1 / Hi),
                                  omega_over_MP=Hi * t_P)

    # required energy splittings for D = 1e-3 in T = 1 s (and 10 s)
    req = {}
    for model in ('GPP_opt', 'GPP_NJP', 'Milburn_tP'):
        for T in (1.0, 10.0):
            omr = required_omega(model, T, 1e-3)
            req[f'{model}_T{T:g}s'] = dict(omega=omr, dE_eV=hbar * omr / eV,
                                           mass_for_dx_10um_kg=hbar * omr / (g_earth * 1e-5))
    res['required_for_D_1e-3'] = req
    return res


def main(out):
    res = dict(part_A=part_A(), part_B=part_B(), constants=dict(t_P=t_P, hbar=hbar))
    dump(res, out)
    A = res['part_A']
    print('Part A (free clock vs Gaussian SW prediction; ideal clock):')
    for r in A['rows']:
        print(f"  w={r['omega']:.1f} T={r['T']:7.1f}: free {r['free_clock']:.6f} gauss {r['gaussian_prediction']:.6f} "
              f"rel {r['rel_diff']} | ideal(same width) {r['ideal_clock_same_width']:.6f} ideal(1e-3) {r['ideal_clock_width_1em3']:.9f}")
    print('  grid checks:', [(c_['analytic'], c_['grid'], c_['abs_diff']) for c_ in A['grid_checks']])
    print('  fixed s0 scan:', [(r['T'], round(r['free'], 6), round(r['ideal'], 6)) for r in A['fixed_s0_scan']])
    B = res['part_B']
    for k, v in B.items():
        print(k)
        for kk, vv in v.items():
            print('   ', kk, vv)


if __name__ == '__main__':
    main(sys.argv[1])
