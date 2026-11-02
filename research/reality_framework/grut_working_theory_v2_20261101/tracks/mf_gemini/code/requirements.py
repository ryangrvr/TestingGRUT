"""mf_gemini: assemble the audit (requirements matrix, holder algebra, S/Q threshold in K, minimal changes).

Reads results/{g1g2, waterfall_dn, waterfall_stoch, waterfall_bg, preheat_lin, ns_hybrid, vscan}.json.
Adds:
  - Monte Carlo Var(delta N) from each case's N(c) curve (the Gauss-Hermite value in waterfall_dn.json
    under-weights the ln|chi_L| singularity); analytic check pi^2/(8 s^2);
  - holder on chi written as u = mu_late + G (chi - v): its size, the hidden-sector tadpole it puts on chi
    during inflation (bias/noise), and the late mu in a -v domain;
  - the K at which M1r's own release would reach S/Q = 3.4716 at e = 1e-6 (R_vec scaled as K^2,
    dq_entropy closed form, everything else fixed) [E];
  - the biased-waterfall seed needed to avoid the delta N tail and wrong-sign domains [E].
Run: python3 -B requirements.py -> results/requirements.json
"""
import json
import os
import numpy as np
from scipy.stats import norm

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')


def load(n):
    with open(os.path.join(RES, n + '.json')) as f:
        return json.load(f)


def main():
    g12 = load('g1g2')
    wd = load('waterfall_dn')
    ws = load('waterfall_stoch')
    wb = load('waterfall_bg')
    pl = load('preheat_lin')
    ns = load('ns_hybrid')
    vs = load('vscan')
    out = {}
    # ---- waterfall summary with MC variance ----
    rng = np.random.default_rng(3)
    gs = np.abs(rng.standard_normal(2_000_000))
    rows = []
    for c in wd['cases']:
        if not c.get('ok'):
            rows.append(dict(label=c['label'], ok=False, reason=c.get('reason')))
            continue
        lc = np.log(np.array(c['Ncurve']['c']))
        Nc = np.array(c['Ncurve']['N_after_sigc'])
        good = np.isfinite(Nc)
        lc, Nc = lc[good], Nc[good]
        s = c['s_eff']
        x = np.log(gs)
        Nx = np.interp(x, lc, Nc)
        lo = x < lc[0]
        Nx[lo] = Nc[0] - (x[lo] - lc[0]) / s
        hi = x > lc[-1]
        Nx[hi] = Nc[-1] - (x[hi] - lc[-1]) / s
        dN = Nx - Nx.mean()
        rows.append(dict(label=c['label'], Delta=c['Delta'], sigma_c=c['sigc'], g=c['g_phys'],
                         N_wf_classical=c['N_end_after_sigc'], N_wf_mean=float(Nx.mean()),
                         varN_MC=float(dN.var()), varN_analytic_pi2_8s2=float(np.pi**2 / 8 / s**2),
                         P_tail_067_MC=float((dN > 0.67).mean()), s_eff=s,
                         Pzeta_peak_linear=c['Pzeta_peak'], N_before_end_peak=c['Pzeta_peak_Nbefore_end'],
                         M_PBH_peak_g=c['M_PBH_peak_g'], f_PBH_tail=c['f_PBH_tail_zc067'],
                         f_PBH_gauss=c['f_PBH_gauss_zc067']))
    out['waterfall'] = rows
    okr = [r for r in rows if r.get('ok', True) and 'N_wf_mean' in r]
    out['waterfall_ranges'] = dict(
        N_wf_mean=[min(r['N_wf_mean'] for r in okr), max(r['N_wf_mean'] for r in okr)],
        Pzeta_peak=[min(r['Pzeta_peak_linear'] for r in okr), max(r['Pzeta_peak_linear'] for r in okr)],
        P_tail_067=[min(r['P_tail_067_MC'] for r in okr), max(r['P_tail_067_MC'] for r in okr)],
        M_PBH_g=[min(r['M_PBH_peak_g'] for r in okr), max(r['M_PBH_peak_g'] for r in okr)],
        f_PBH_tail=[min(r['f_PBH_tail'] for r in okr), max(r['f_PBH_tail'] for r in okr)],
        stochastic_tail_067=[min(r['tails']['0.67'] for r in ws['runs']), max(r['tails']['0.67'] for r in ws['runs'])],
        stochastic_varN=[min(r['varN'] for r in ws['runs']), max(r['varN'] for r in ws['runs'])])
    # ---- holder on chi ----
    H = 1.3e11
    MP = 2.435e18
    hub = H / MP
    v = np.sqrt(8.0 / 3.0)
    mu = 7.159456
    u_inf = 10j
    Gv = mu - u_inf                         # u_late(+v) = mu_bare + G v = mu (real), u_inf = mu_bare = 10 i
    g_hold = abs(Gv) * hub / v              # dimensionless: |G| v = 12.3 H
    Y, S = 6.888, 28.576                    # prepared |Y|, |Sigma| in H (M1r at lam = 0.5)
    tad = 2 * g_hold * S * Y**2             # in H^3 (max |dV/dchi| from -(g e^{i a} chi S* Y^2 + h.c.))
    drift = tad / 3.0                        # chi displacement per e-fold in H units (light chi)
    noise = 1 / (2 * np.pi)
    u_minus = u_inf - Gv                     # late mu in a -v domain
    out['holder_on_chi'] = dict(
        form='-(g e^{i alpha} (chi - v) Sigma* Y^2 + h.c.): u = mu_late + G (chi - v); u_inflation = mu_late - G v',
        g=float(g_hold), Gv_over_H=[float(Gv.real), float(Gv.imag)],
        tadpole_bias_over_noise_per_efold=float(drift / noise),
        minus_v_domain=dict(mu_late_over_H=[float(u_minus.real), float(u_minus.imag)], abs_mu=float(abs(u_minus)),
                            CP_phase_arg_uinf_over_ulate=float(np.angle(u_inf / u_minus)),
                            plus_domain_phase=float(np.angle(u_inf / mu))),
        m_chi_over_H_end=wb['m_chi_over_H_end'], release_amplitude_first_overshoot_over_v=wb['representative']['osc_amp_over_v_first'],
        note='label D for the algebra, E for the tadpole (prepared state of M1r used as a stand-in)')
    # ---- S/Q threshold in K on M1r's own release (e = 1e-6) ----
    EEc = 15.60
    W = ((3.4716 + 0.00428) / 0.190255)**(1 / 0.75013) + 0.97494
    Rtot = W / EEc - 1
    Rvec = Rtot - 0.379 - 0.019
    K = 1.986 * np.sqrt(Rvec / 3.608)
    tK = np.interp(K, [wb['K_if_sudden_decay_at_t_Hend'][k] for k in ['0.5', '1', '2', '4', '8', '16']], [0.5, 1, 2, 4, 8, 16])
    out['SQ_threshold'] = dict(W_crit=float(W), R_total_crit=float(Rtot), R_vec_crit=float(Rvec), K_crit=float(K),
                               sudden_decay_time_for_K_crit_Hend=float(tK),
                               label='E: R_vec ~ K^2 at fixed release (WORKING_THEORY 4b); the release itself is not M1r''s here')
    # ---- biased waterfall: seed to remove the tail and the wrong-sign domains ----
    s_rate = 1.098
    chi_PBH = 1 / (2 * np.pi * s_rate * np.sqrt(1e-3))     # P_zeta <= 1e-3 at onset (H units)
    Nsig = 36.0
    nsig_needed = norm.isf(np.exp(-3 * Nsig) * 1e-2)          # < 1e-2 wrong-sign patches in e^{3 N_sigma} volumes
    sig_noise = np.sqrt(5.0) / (2 * np.pi)                    # ~ 5 e-folds of light-field noise before the onset
    chi_dom = nsig_needed * sig_noise
    chib = max(chi_PBH, chi_dom)
    Nwf_b = np.log(0.6466 * v / (chib * hub)) / s_rate
    out['biased_waterfall'] = dict(chi_seed_for_Pzeta_1em3_H=float(chi_PBH), sigmas_for_no_wrong_domain=float(nsig_needed),
                                   chi_seed_for_no_wrong_domain_H=float(chi_dom), chi_seed_H=float(chib),
                                   chi_seed_GeV=float(chib * H), N_wf_classical=float(Nwf_b),
                                   Pzeta_onset=float((1 / (2 * np.pi * s_rate * chib))**2), label='E')
    # ---- requirements matrix ----
    G3 = g12['G3_check']['H=1.3e+11,Gamma=0.5H']
    out['requirements'] = {
        'm_chi = 3 H_i with lam v^4/4 = V0': dict(status='forces v = 1.633 M_P, lam = 4.8e-15 (H = 1.3e11); m_chi/H_end = 4.21',
                                                   label='D/C'),
        'Gamma_eff/H_i >= 0.5': dict(
            status='fails for every channel with a radiatively acceptable coupling',
            perturbative_yN=f"M_N = y v = {g12['G2']['H=1.3e+11']['y=2.0']['MN_vac_over_mchi_half']:.2e} x m_chi/2 (y = 2); CW residual {g12['G2']['H=1.3e+11']['y=2.0']['CW_residual_over_V0_Z2']:.1e} V0; y_cap = {g12['G2']['H=1.3e+11']['y_cap_CW_1pct']:.1e} -> Gamma <= {g12['G2']['H=1.3e+11']['Gamma_over_H_at_ycap_Majorana_unblocked']:.1e} H",
            fermion_preheating_tuned_M0=f"GK energy per crossing {g12['preheating_caps']['y=cap (shifted, CW<=1%)']['GK_energy_fraction_per_crossing']:.1e} at the loop cap (Pauli limited)",
            bosonic_preheating_cap=f"fraction by 2 H_end^-1 <= {pl['bosonic']['cap']['frac_at_2Hend']['max']:.1e}; 50% at {pl['bosonic']['cap']['t_frac_0p5_Hend']['max']:.1f}-{pl['bosonic']['cap']['t_frac_0p5_Hend']['wkb']:.1f} H_end^-1 (optimistic, linear)",
            higgs_portal=f"Gamma = 0.5 H needs lam_chiH = {g12['Higgs_portal']['Gamma=0.5H']['lam_chiH']:.2e}: Higgs tachyonic in inflation (|m^2|/H^2 = {g12['Higgs_portal']['Gamma=0.5H']['abs_mH2_inflation_over_H2']:.1e}) or tachyonic windows (exponent {g12['Higgs_portal']['Gamma=0.5H']['tachyonic_window_exponent_per_half_cycle']:.0f}); lam_H(full) = {g12['Higgs_portal']['Gamma=0.5H']['lamH_full']:.2f}, Landau pole {g12['Higgs_portal']['Gamma=0.5H']['LandauPole_GeV']:.1e} GeV",
            tachyonic_preheating_of_chi='the tachyonic rate is <= 1.1 H: Hubble-limited, no sub-horizon fragmentation; energy stays in massive chi (w -> 0)',
            label='C/D/E'),
        'H_i in [1.25, 1.55]e11 GeV': dict(status='scales trivially; r = 2.5-3.9e-7', label='C'),
        'n_s 0.967-0.973': dict(status='needs non-convex V_SR: hilltop |eta| = 0.0136-0.0166 (alpha_s ~ -4e-9) or log with N_sigma = 30-37, i.e. N_wf <= 23 (alpha_s ~ -(0.7-1.1)e-3); convex V_SR gives n_s > 1',
                                label='C'),
        'no defects': dict(status='fails for real chi: walls with v = 1.633 M_P ~ the topological-inflation threshold 0.33 m_Pl = 1.655 M_P; equivalently the delta N tail (P(zeta > 0.67) = 0.20-0.29) and PBH overproduction (f_PBH 1e6-1e16)',
                           label='C/Q'),
        'holder working': dict(status='possible on chi as u = mu_late + G(chi - v) (on during inflation, quiet m = 4.2 H_end), but +-v domains carry opposite CP phases and different late mu; release = 15% overshoot after a 17-40 e-fold relaxation, not M1r''s',
                               label='D/A'),
    }
    out['minimal_changes'] = dict(
        heavy_waterfall_G3=dict(lam=G3['lam'], m_over_H=G3['m_over_H'], v_over_MP=G3['v_over_MP'], y_max=G3['y_max'],
                                loop_over_lam=G3['y4_16pi2_over_lam']),
        note='requires moving the holder off chi (spectator, G5), plus a bias or complex/gauged chi for defects')
    with open(os.path.join(RES, 'requirements.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in out if k != 'waterfall'}, indent=1))
    for r in rows:
        print(r)


if __name__ == '__main__':
    main()
