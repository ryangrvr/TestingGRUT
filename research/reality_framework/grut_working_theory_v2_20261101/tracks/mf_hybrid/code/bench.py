"""mf_hybrid benchmark driver: inflation (one-loop) -> waterfall (homogeneous) -> axial reheating fluid -> export.

Usage: python3 -B bench.py HFa | HFb      -> ../results/bench_<name>.json (+ ../results/axial_tables.json cache)
Slice convention (M1r): a = 1 and H = H_i at the critical point sigma = sigma_c, the end of slow roll (eps_H ~ 1e-7,
R = 12 H_i^2: the static-slice analogue). The eps_H = 1 time is reported too. K = a^2 H / H_i in the radiation era.
"""
import math, json, os, sys, time
import numpy as np
import mfh_common as mc
import inflation as I
import waterfall as W
import reheat as RH

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
T_SW = 0.6     # H_i^-1 after the slice: end of the homogeneous stage

SPECS = {
    'HFa': dict(H_GeV=1.30e11, m_chi=250.0, m_sig=600.0, kappa=10.0, ns=mc.NS_PLANCK, c_bias=1e5, smech='a',
                s_star_MP=3e-3, sig_ch=[(250.0, 1200.0)], chi_ch=[(400.0, 2000.0), (200.0, 2000.0), (120.0, 2000.0)],
                note='Planck-centred n_s; spectator mass from the waterfall, m_s = g_s chi'),
    'HFb': dict(H_GeV=1.30e11, m_chi=250.0, m_sig=600.0, kappa=10.0, ns=mc.NS_ACT, c_bias=1e5, smech='b',
                s_star_MP=3e-3, sig_ch=[(250.0, 1200.0)], chi_ch=[(400.0, 2000.0), (200.0, 2000.0), (120.0, 2000.0)],
                note='ACT-centred n_s; spectator mass from the curvature coupling, m_s^2 = 9 H_i^2 + xi_s R, xi_s = -3/4'),
}


def tables(spec):
    path = os.path.join(RES, 'axial_tables.json')
    cache = json.load(open(path)) if os.path.exists(path) else {}
    need = [(M / spec['m_sig']) for M, f in spec['sig_ch']] + [(M / spec['m_chi']) for M, f in spec['chi_ch']]
    for Mh in need:
        key = '%.6f' % Mh
        if key not in cache:
            t0 = time.time()
            cache[key] = RH.eps_table(Mh)
            print('table Mh=%.4f %.1fs' % (Mh, time.time() - t0), flush=True)
            json.dump(cache, open(path, 'w'))
    return cache


def run(name, Nback=46.0, iters=3, As_t=None, ns_t=None):
    spec = SPECS[name]
    p = mc.derived(spec['H_GeV'], spec['m_chi'], spec['m_sig'])
    mu2 = spec['kappa'] * p['g'] ** 2 * (2 * p['sigma_c']) ** 2
    tab = tables(spec)
    chans_sig = [RH.Channel(spec['m_sig'], M, f, tab['%.6f' % (M / spec['m_sig'])]) for M, f in spec['sig_ch']]
    chans_chi = [RH.Channel(spec['m_chi'], M, f, tab['%.6f' % (M / spec['m_chi'])]) for M, f in spec['chi_ch']]
    hist = []
    for it in range(iters):
        inf, o, _ = I.calibrate(p, ns_t or spec['ns'], Nback, mu2=mu2, c_bias=spec['c_bias'], As_t=As_t)
        s6 = I.sigma_of_N(inf, 6.0)
        model = W.Model(p, inf, spec['c_bias'], smech=spec['smech'], m_s=3.0, xi_s=-0.75)
        model.chans_sig, model.chans_chi = chans_sig, chans_chi
        rec, _, _ = W.run(model, s6, T_SW, nrec=3000, t_modes=T_SW)
        Hi = float(np.interp(0.0, rec['t'], rec['H']))
        j1 = int(np.argmax((rec['t'] > 0) & (rec['eps'] >= 1.0)))
        t_end = float(rec['t'][j1])
        jsw = len(rec['t']) - 1                     # switch to fluids at T_SW (after the linear modes go non-linear)
        Esig, Etot = float(rec['E_sig'][jsw]), float(rec['E_tot'][jsw])
        # E_N: production energy scale (3/4) mu5_max of the dominant chi channel at the switch
        Sig_chi = math.sqrt(2 * (Etot - Esig)) / spec['m_chi']
        EN0 = 0.75 * max(spec['m_chi'] * Sig_chi / (2 * f) for M, f in spec['chi_ch'])
        gs = model.g_s
        v, mchi = p['v'], spec['m_chi']
        if spec['smech'] == 'a':
            s_model = lambda t, rc, Rb: gs ** 2 * (v * v + rc / mchi ** 2)
        else:
            s_model = lambda t, rc, Rb: 9.0 - 0.75 * Rb
        st = dict(t=float(rec['t'][jsw]), lna=float(rec['lna'][jsw]), rho_sig=Esig, rho_chi=Etot - Esig,
                  rho_r=float(rec['rho_r'][jsw]),
                  s=float(rec['s_over_sstar'][jsw]), sd=float(rec['sdot_over_sstar'][jsw]), E_N=EN0)
        fl = RH.run_fluid(st, p['mP'], chans_sig, chans_chi, spec['H_GeV'], s_model, t_end=3e4)
        K = float(fl['K'][-1] / Hi)
        Hstar = o['H']
        Nst = mc.Nstar_from_K(Hstar * spec['H_GeV'], Hi * spec['H_GeV'], K)
        hist.append(dict(iter=it, Nback=Nback, K=K, Nstar=Nst))
        print('%s iter %d Nback %.3f -> K %.4f N* %.3f  t_end %.3f  Esig/E %.3f' % (name, it, Nback, K, Nst, t_end, Esig / Etot),
              flush=True)
        if abs(Nst - Nback) < 0.02:
            break
        Nback = Nst
    msr, _ = I.ms(inf, Nback)
    if As_t is None and (abs(msr['As'] / mc.AS - 1) > 3e-3 or abs(msr['ns'] - spec['ns']) > 3e-4):
        # second pass: retarget the slow-roll calibration so that the Mukhanov-Sasaki A_s and n_s hit the targets
        return run(name, Nback=Nback, iters=2, As_t=mc.AS ** 2 / msr['As'], ns_t=spec['ns'] + (spec['ns'] - msr['ns']))
    return dict(spec=spec, p=p, mu2=mu2, inf=inf, o=o, ms=msr, rec=rec, fl=fl, Hi=Hi, t_end=t_end, K=K, hist=hist,
                model=model, chans_sig=chans_sig, chans_chi=chans_chi)


def summarize(R, name):
    spec, p, o, ms, rec, fl, Hi = R['spec'], R['p'], R['o'], R['ms'], R['rec'], R['fl'], R['Hi']
    HG = spec['H_GeV'] * Hi
    mP = p['mP']
    # fractions and timings in the fluid
    rho = fl['rho_sig'] + fl['rho_chi'] + fl['rho_N'] + fl['rho_SM']
    frad = (fl['rho_N'] + fl['rho_SM']) / rho
    fsm = fl['rho_SM'] / rho
    t = fl['t']

    def first(mask):
        idx = np.argmax(mask)
        return float(t[idx]) if mask.any() else None
    t_eq = first(frad >= 0.5)
    t_99 = first(frad >= 0.99)
    t_sm3 = first(fsm >= 1e-3)
    t_sm50 = first(fsm >= 0.5)
    rhoGeV = lambda x: x * spec['H_GeV'] ** 4
    T_eq = mc.T_of_rho(rhoGeV(np.interp(t_eq, t, fl['rho_N'] + fl['rho_SM']))) if t_eq else None
    T_max = float(np.max(mc.T_of_rho(rhoGeV(fl['rho_N'] + fl['rho_SM']))))
    TK = mc.T_of_rho(3 * R['K'] ** 2 * mP ** 2 * Hi ** 2 * spec['H_GeV'] ** 4)   # radiation-extrapolated T at a = 1
    # spectator
    sst = spec['s_star_MP'] * mP
    ms_late = float(fl['m_s'][-1])
    frac_s_slice = 0.5 * ms_late ** 2 * sst ** 2 / (3 * mP ** 2 * Hi ** 2)
    # s fraction growth in radiation ~ a (until it decays)
    G_s = 1.3e-6
    out = dict(name=name, label='C (given the model choices) unless a field says E', note=spec['note'],
               slice='a = 1, H = H_i at sigma = sigma_c (end of slow roll; eps_H ~ 1e-7; R = 12 H_i^2), M1r static-slice analogue',
               H_i_GeV=HG, H_ref_GeV=spec['H_GeV'], Hi_over_Href=Hi,
               params=dict(lam=p['lam'], v_MP=p['v_MP'], g=p['g'], sigma_c_MP=p['sigma_c_MP'], X=p['X'],
                           m_chi_over_Hi=spec['m_chi'] / Hi, m_sig_over_Hi=spec['m_sig'] / Hi, c_bias_Hi3=spec['c_bias'],
                           c1_Hi3=o['c1'], m2_Hi2=o['m2'], mu_R2_Hi2=R['mu2'], kappa=spec['kappa'],
                           V0_quarter_GeV=(3 * mP ** 2) ** 0.25 * spec['H_GeV'], g_s=R['model'].g_s if spec['smech'] == 'a' else 0.0,
                           xi_s=-0.75 if spec['smech'] == 'b' else 0.0, m_s0_over_Hi=3.0 if spec['smech'] == 'b' else 0.0,
                           s_star_MP=spec['s_star_MP'], s_star_GeV=sst * spec['H_GeV'],
                           holder='-(g_h e^{i alpha} s Sigma* Y^2 + h.c.), |g_h| s_* = |10i - mu| = 12.30 H_i -> |g_h| = %.4g' % (12.2987 / sst),
                           axial_sigma=[dict(M_over_Hi=M, f_over_Hi=f, Gamma_pert=c.Gp) for (M, f), c in zip(spec['sig_ch'], R['chans_sig'])],
                           axial_chi=[dict(M_over_Hi=M, f_over_Hi=f, Gamma_pert=c.Gp) for (M, f), c in zip(spec['chi_ch'], R['chans_chi'])],
                           smech=spec['smech']),
               cmb=dict(Nstar=R['hist'][-1]['Nstar'], Nback_used=R['hist'][-1]['Nback'], As_sr=o['As'], ns_sr=o['ns'],
                        alpha_s_sr=o['alpha_s'], r_sr=o['r'], ms=ms,
                        ns_pull_planck=(ms['ns'] - mc.NS_PLANCK) / mc.NS_PLANCK_ERR, ns_pull_ACT=(ms['ns'] - mc.NS_ACT) / mc.NS_ACT_ERR,
                        sigma_star_MP=o['sigma_star_MP'], Delta_sigma_MP=o['Delta_sigma_MP'], monotonic=o['monotonic'],
                        tuning=o['tuning'], eta=o['eta']),
               waterfall=dict(t_eps1=R['t_end'], N_eps1=float(np.interp(R['t_end'], rec['t'], rec['lna'])),
                              H_eps1_over_Hi=float(np.interp(R['t_end'], rec['t'], rec['H'])) / Hi,
                              Esig_frac_at_eps1=float(np.interp(R['t_end'], rec['t'], rec['E_sig'] / rec['E_tot'])),
                              mode_growth_at_t_modes={k: float(np.nanmax(np.abs(v))) for k, v in rec['modes'].items()},
                              t_modes_end=rec['t_modes_end']),
               reheating=dict(K=R['K'], K_history=R['hist'], t_rad_eq=t_eq, t_rad_99=t_99, T_at_eq_GeV=T_eq, T_max_GeV=T_max,
                              T_K_GeV=TK, t_SM_1em3=t_sm3, t_SM_half=t_sm50, GammaN_rest=fl['GammaN_rest_avg'],
                              label='E: non-perturbative axial rates from the [C] half-period table, no back-reaction on the modes'),
               spectator=dict(m_s_late=ms_late, energy_fraction_after_turn_on=frac_s_slice,
                              iso_dlnS_ds=1.0 / (2 * math.pi * sst), beta_iso_per_p2=(1 / (2 * math.pi * sst)) ** 2 / mc.AS,
                              Gamma_s_assumed=G_s))
    out['extra_params'] = dict(
        higgs='no tree-level portal of the Higgs to sigma, chi or s (lam_chiH = lam_sigmaH = lam_sH = 0, model choice); xi_H = 0.1 (M1r)',
        N_singlets='four Majorana singlets, -(M_j/2) Nbar_j exp(i gamma5 phi/f_j) N_j with phi = sigma (j=1) or chi (j=2,3,4); '
                   'seesaw y_nu^2 = m_nu M/v_EW^2 with m_nu = 0.05 eV, v_EW = 174.1 GeV; N -> L H',
        s_decay='y_s s Nbar_s N_s with y_s = 4.6e-3 (light N_s, N_s -> i gamma5 N_s keeps s -> -s), Gamma_s = 1.26e-6 H_i after unblocking at a ~ 31 [E]',
        hidden_sector='M1r at lambda_Y = lambda_Sigma = 0.5 (BRIEF_DECISIVE), holder on s instead of the inflaton; no sigma/chi-hidden portals (sequestered by choice)',
        cosmological_constant='constant added so that V(vacuum) = 0',
        slice_note='t = 0, N = 0 at sigma = sigma_c; eps_H = 1 at t = %.4f with H/H_i = %.4f' % (out['waterfall']['t_eps1'], out['waterfall']['H_eps1_over_Hi']))
    out['status'] = 'final for the stated model choices (homogeneous + fluid; the reheating rates are [E]); see STATUS.md'
    return out


def export(R, name):
    S = summarize(R, name)
    rec, fl, Hi = R['rec'], R['fl'], R['Hi']
    mP = R['p']['mP']
    sst = R['spec']['s_star_MP'] * mP
    # homogeneous part up to t_end
    m1 = rec['t'] <= T_SW + 1e-12
    idx = np.where(m1)[0]
    idx = idx[np.unique(np.linspace(0, len(idx) - 1, 1600).astype(int))]
    th = rec['t'][idx]
    Etot = rec['E_tot'][idx] + rec['rho_r'][idx]
    ms2h = np.array([R['model'].ms2(x, Rr) for x, Rr in zip(rec['chi'][idx], rec['R'][idx])])
    rho_s_h = sst ** 2 * (0.5 * rec['sdot_over_sstar'][idx] ** 2 + 0.5 * ms2h * rec['s_over_sstar'][idx] ** 2)
    part1 = dict(t=th, N=rec['lna'][idx], H=rec['H'][idx], eps_H=rec['eps'][idx], R=rec['R'][idx],
                 sigma_MP=rec['sigma'][idx] / mP, chi_MP=rec['chi'][idx] / mP, s_over_sstar=rec['s_over_sstar'][idx],
                 m_s=np.sqrt(np.maximum(ms2h, 0)), f_sig=rec['E_sig'][idx] / Etot, f_chi=(rec['E_tot'][idx] - rec['E_sig'][idx]) / Etot,
                 f_rad=rec['rho_r'][idx] / Etot, f_SM=0 * th, f_s=rho_s_h / Etot)
    # fluid part
    jf = np.unique(np.concatenate([np.linspace(0, 2999, 900).astype(int), np.arange(3000, len(fl['t']), 2)]))
    rho = fl['rho_sig'] + fl['rho_chi'] + fl['rho_N'] + fl['rho_SM']
    s_out = np.where(np.isnan(fl['s_over_sstar']), fl['s_env'], fl['s_over_sstar'])
    rho_s = 0.5 * fl['m_s'] ** 2 * sst ** 2 * np.where(np.isnan(fl['s_over_sstar']), fl['s_env'] ** 2, fl['s_over_sstar'] ** 2)
    part2 = dict(t=fl['t'][jf], N=fl['lna'][jf], H=fl['H'][jf], eps_H=fl['eps'][jf], R=fl['R_mean'][jf],
                 sigma_MP=np.full(len(jf), np.nan), chi_MP=np.full(len(jf), np.nan), s_over_sstar=s_out[jf], m_s=fl['m_s'][jf],
                 f_sig=(fl['rho_sig'] / rho)[jf], f_chi=(fl['rho_chi'] / rho)[jf], f_rad=((fl['rho_N'] + fl['rho_SM']) / rho)[jf],
                 f_SM=(fl['rho_SM'] / rho)[jf], f_s=(rho_s / rho)[jf])
    H = dict()
    for k in part1:
        H[k] = np.concatenate([part1[k], part2[k][1:]])
    # H_i units
    H['t'] = H['t'] * Hi
    H['H'] = H['H'] / Hi
    H['R'] = H['R'] / Hi ** 2
    H['a'] = np.exp(H['N'])
    H['K_running'] = H['a'] ** 2 * H['H']
    H['R_osc_sig'] = np.concatenate([np.full(len(th), np.nan), (fl['R_osc_sig'][jf] / Hi ** 2)[1:]])
    H['R_osc_chi'] = np.concatenate([np.full(len(th), np.nan), (fl['R_osc_chi'][jf] / Hi ** 2)[1:]])
    H['s_is_envelope'] = np.concatenate([np.zeros(len(th)), np.isnan(fl['s_over_sstar'][jf][1:]).astype(float)])
    hist = {k: [None if (isinstance(x, float) and math.isnan(x)) else float('%.7g' % x) for x in np.asarray(v, dtype=float)]
            for k, v in H.items()}
    S['history_columns'] = dict(t='cosmic time in H_i^-1 from the slice', N='ln a (a = 1 at the slice)', a='scale factor',
                                H='H/H_i', eps_H='-Hdot/H^2', R='Ricci scalar / H_i^2 (exact in the homogeneous part; cycle mean in the fluid part)',
                                R_osc_sig='amplitude of the 2 m_sig oscillation of R (fluid part)', R_osc_chi='amplitude at 2 m_chi',
                                sigma_MP='clock (homogeneous part)', chi_MP='waterfall (homogeneous part)',
                                s_over_sstar='holder spectator s/s_* (exact to t ~ 60, then its WKB envelope; see s_is_envelope)',
                                m_s='spectator mass / H_i', f_sig='energy fraction in the sigma condensate', f_chi='in the chi condensate',
                                f_rad='in relativistic N + SM radiation', f_SM='in SM radiation (N -> L H)', f_s='in the spectator',
                                K_running='a^2 H / H_i (converges to K in radiation)')
    S['history'] = hist
    path = os.path.join(RES, 'bench_%s.json' % name)
    json.dump(S, open(path, 'w'), separators=(',', ':'))
    print('wrote', path, os.path.getsize(path) / 1e6, 'MB')
    return S


if __name__ == '__main__':
    for nm in sys.argv[1:]:
        t0 = time.time()
        R = run(nm)
        S = export(R, nm)
        print(json.dumps({k: S[k] for k in ['H_i_GeV', 'cmb', 'waterfall', 'reheating', 'spectator']}, indent=1, default=str)[:4000])
        print('%.0fs' % (time.time() - t0))
