"""Diagnostics for the exported benchmarks: defects (bias vs quantum spread), waterfall curvature, R spectrum and the
Y resonance, decay channels, gravitational production, Higgs sign structure, spectator energy, isocurvature [C/E].
Usage: python3 -B diagnostics.py HFa HFb  -> ../results/diagnostics.json
"""
import math, json, os, sys
import numpy as np
import mfh_common as mc
import inflation as I
import waterfall as W
import reheat as RH
import bench as B

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')


def rebuild(name):
    d = json.load(open(os.path.join(RES, 'bench_%s.json' % name)))
    spec = B.SPECS[name]
    p = mc.derived(spec['H_GeV'], spec['m_chi'], spec['m_sig'])
    pr = d['params']
    inf = mc.Inflation(p, pr['c1_Hi3'], pr['m2_Hi2'], c_bias=spec['c_bias'], mu2=pr['mu_R2_Hi2'])
    tab = B.tables(spec)
    cs = [RH.Channel(spec['m_sig'], M, f, tab['%.6f' % (M / spec['m_sig'])]) for M, f in spec['sig_ch']]
    cc = [RH.Channel(spec['m_chi'], M, f, tab['%.6f' % (M / spec['m_chi'])]) for M, f in spec['chi_ch']]
    return d, spec, p, inf, cs, cc


def defects(name):
    d, spec, p, inf, cs, cc = rebuild(name)
    ks = (0.0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0)
    model = W.Model(p, inf, spec['c_bias'], smech=spec['smech'], m_s=3.0, xi_s=-0.75)
    model.chans_sig, model.chans_chi = cs, cc
    rec, _, _ = W.run(model, I.sigma_of_N(inf, 1.0), 0.6, ks=ks, nrec=1500, t_modes=0.6)
    t = rec['t']
    m0 = t >= 0
    kphys = np.array(ks[1:]) * p['m_chi']
    out = []
    for tt in (0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6):
        j = int(np.searchsorted(t, tt))
        if j >= len(t):
            j = len(t) - 1
        G = np.array([abs(rec['modes'][str(k)][j]) for k in ks[1:]])
        if np.any(np.isnan(G)):
            continue
        # vacuum spectrum at the slice: massless-like (M_chi^2(sigma_c) = 0), |chi_k|^2 = 1/(2k); grid integral in k
        kk = np.linspace(kphys[0], kphys[-1], 400)
        Gk = np.exp(np.interp(kk, kphys, np.log(G)))
        dchi2 = np.trapezoid(kk ** 2 * Gk ** 2 / (2 * kk), kk) / (2 * math.pi ** 2)
        # add the IR part k < kphys[0] with G(k -> 0) ~ G(0)
        G0 = abs(rec['modes']['0.0'][j])
        dchi2 += G0 ** 2 * kphys[0] ** 2 / (4 * math.pi ** 2)
        chi_cl = float(rec['chi'][j])
        out.append(dict(t=float(t[j]), chi_cl_over_v=chi_cl / p['v'], dchi_rms_over_v=math.sqrt(dchi2) / p['v'],
                        ratio=chi_cl / math.sqrt(dchi2), G0=float(G0)))
    # wall probability: Gaussian tail with the minimum ratio over the linear stage (before chi_cl ~ v)
    rmin = min(o['ratio'] for o in out if o['chi_cl_over_v'] < 0.9) if out else None
    # number of waterfall correlation volumes in today's horizon ~ exp(3 N_*) (m_chi/H)^3
    lnV = 3 * d['cmb']['Nstar'] + 3 * math.log(spec['m_chi'])
    lnP = -0.5 * rmin ** 2 - math.log(rmin * math.sqrt(2 * math.pi)) if rmin else None
    # waterfall duration and an upper estimate of the curvature from the time delay: zeta ~ H dchi/chi_dot
    jw = int(np.argmax(m0 & (rec['chi'] > 0.5 * p['v'])))
    tw = float(t[jw])
    zeta_est = (1.0 / rmin) * tw if rmin else None
    # R spectrum in the homogeneous stage
    mR = (t > 0.2) & (t < 0.6)
    tR, RR = t[mR], rec['R'][mR]
    tu = np.linspace(tR[0], tR[-1], 4096)
    Ru = np.interp(tu, tR, RR) - np.mean(RR)
    F = np.abs(np.fft.rfft(Ru * np.hanning(len(Ru))))
    w = 2 * math.pi * np.fft.rfftfreq(len(tu), tu[1] - tu[0])
    top = [float(w[i]) for i in np.argsort(F)[::-1][:6]]
    band = (w > 10) & (w < 60)
    frac_band = float(np.sum(F[band] ** 2) / np.sum(F[w > 5] ** 2))
    xiH_R_min = float(0.1 * np.min(rec['R'][m0]))
    return dict(linear_stage=out, ratio_min=rmin, ln_P_wrong_vacuum_per_volume=lnP, ln_volumes_today=lnV,
                walls=('none expected: ln P + ln N_vol = %.1f' % (lnP + lnV)) if lnP else None,
                t_waterfall_chi_half_v=tw, zeta_waterfall_upper_E=zeta_est,
                R_spectrum_top_omega_Hi=top, R_power_fraction_10_60_Hi=frac_band,
                xiH_R_min_Hi2=xiH_R_min, R_min=float(np.min(rec['R'][m0])), R_max_after_slice=float(np.max(rec['R'][m0])))


def channels(name):
    d, spec, p, inf, cs, cc = rebuild(name)
    rho0 = 3 * p['mP'] ** 2
    rows = []
    for nm, chans, m, frac in (('sigma', cs, spec['m_sig'], 0.15), ('chi', cc, spec['m_chi'], 0.85)):
        for c in chans:
            GN = RH.MNU * (c.M * spec['H_GeV']) ** 2 / (8 * math.pi * RH.V_EW ** 2) / spec['H_GeV']
            Sig0 = math.sqrt(2 * frac * rho0) / m
            rows.append(dict(field=nm, m=m, M_N=c.M, f=c.f, y_p=c.M / c.f, Gamma_pert=c.Gp,
                             Gamma_nonpert_at_slice_amplitude=c.rate(frac * rho0), A_slice=Sig0 / (2 * c.f),
                             mu5_max_over_4pi_f=(m * Sig0 / (2 * c.f)) / (4 * math.pi * c.f),
                             f_edge=math.sqrt(m * Sig0 / (8 * math.pi)), M_N_GeV=c.M * spec['H_GeV'],
                             ynu2_seesaw=RH.MNU * c.M * spec['H_GeV'] / RH.V_EW ** 2, Gamma_N_rest=GN))
    # gravitational: X pairs from xi_X R_osc at 2 m_chi (perturbative, V = (A/2) cos(w t) X^2, A = |xi_X| R_osc,
    # energy rate A^2 m / (64 pi) per real dof); R_osc = 3 rho_chi/M_P^2 ~ 3 * 0.85 * 3 H^2 initially
    A = (5.0 / 6.0) * 3 * 0.85 * 3
    dE = 2 * A ** 2 * spec['m_chi'] / (64 * math.pi)          # H_i^5 per real dof x 2 (complex X)
    grav = dict(xiX_Rosc_amplitude_Hi2=A, X_energy_rate_Hi5=dE, X_energy_per_Hubble_Hi4=dE,
                compare_M1r_EX_conformal=1.887e5, ratio=dE / 1.887e5,
                note='E: X pairs of energy m_chi from the 2 m_chi oscillation of R; Y, Sigma, A minimally coupled: '
                     'conformal-breaking production ~ (m/M_P)^4, negligible. Background frequencies 2 m_chi = %.0f, '
                     '2 m_sig = %.0f, m_s = 3 H_i: none near omega_YL = 22.25 or 2 omega_YL = 44.5' % (2 * spec['m_chi'], 2 * spec['m_sig']))
    # spectator
    sst = spec['s_star_MP'] * p['mP']
    gh = 12.2987 / sst
    ys = 4.6e-3
    a_unb = (sst * 2 * ys / 3.0) ** (2.0 / 3.0)
    Gs = ys ** 2 * 3.0 / (16 * math.pi)
    K = d['reheating']['K']
    t_dec = 1 / Gs
    a_dec = math.sqrt(2 * K * t_dec)
    f_dec = 1.5 * spec['s_star_MP'] ** 2 * a_dec
    spect = dict(s_star_Hi=sst, s_star_GeV=sst * spec['H_GeV'], g_holder=gh, energy_fraction=1.5 * spec['s_star_MP'] ** 2,
                 decay='E: y_s s Nbar_s N_s (N_s -> i gamma5 N_s keeps s -> -s), y_s = %.2g (s stays light: one-loop '
                       'mass <= 0.01 H^2), unblocked at a = %.1f, Gamma_s = %.2g H_i' % (ys, a_unb, Gs),
                 Gamma_s=Gs, a_decay=a_dec, energy_fraction_at_decay=f_dec, entropy_injection_E=0.75 * f_dec,
                 into_hidden='none: s -> Y Y closed (m_s = 3 < 2 m_YL), s -> X via loops only',
                 iso_amplitude=1 / (2 * math.pi * sst), beta_iso_per_p2=(1 / (2 * math.pi * sst)) ** 2 / mc.AS,
                 one_loop_holder_mass2_Hi2=gh ** 2 * 28.6 ** 2 / (16 * math.pi ** 2),
                 holder_tadpole_force_Hi3=2 * gh * 28.6 * 6.9 ** 2, drift_50efolds_Hi=2 * gh * 28.6 * 6.9 ** 2 / 3 * 50)
    return dict(channels=rows, gravitational=grav, spectator=spect)


if __name__ == '__main__':
    path = os.path.join(RES, 'diagnostics.json')
    out = json.load(open(path)) if os.path.exists(path) else {}
    for nm in sys.argv[1:]:
        out[nm] = dict(defects=defects(nm), **channels(nm))
        print(nm, json.dumps(out[nm]['defects'], indent=0)[:1500])
        print(json.dumps(out[nm]['channels'], indent=0)[:2500])
        print(json.dumps(out[nm]['spectator'], indent=0))
    json.dump(out, open(path, 'w'), indent=1)
