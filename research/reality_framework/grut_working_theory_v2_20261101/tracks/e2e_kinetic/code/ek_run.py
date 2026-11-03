#!/usr/bin/env python3
"""e2e_kinetic driver: a full post-waterfall history (coherent mode stage -> dephased kinetic stage) for HFa / HFb.

Pipeline (all model choices; labels in NOTES.md):
  1. start state at t0 (default 0.19 H_i^-1, the first chi turning point after the waterfall, as the mf_hybrid skeptic's
     fluid): energy fractions of the exported homogeneous bench history; the 0.9% the track's unblocked friction had
     moved to radiation by t0 is returned to chi, and the N modes start in the (superadiabatic) vacuum; or a state
     supplied by the lattice (--state, interface in NOTES.md section 6);
  2. COHERENT stage t0 -> t_coh: every N mode evolved exactly (ek_modes.c), Pauli blocking automatic (n <= 1),
     collisions (decay + in-medium) at the parent's turning points, back-reaction on the condensates;
  3. re-binning of the occupations (number-conserving) on a log comoving grid, DEPHASED stage t_coh -> end: rate
     equations with the exact single half-period transition tables (ek_ptable.py) and the golden-rule shell for the
     perturbative channels;
  4. SM sector throughout: hard SM gas from N removal, relaxation into the thermal bath at the Kurkela-Lu rate.
Usage: python3 -B ek_run.py check            (constant-amplitude regressions vs mf_hybrid skeptic's exact modes)
       python3 -B ek_run.py run NAME TAG [--H 1.3e11] [--mnu 0.05] [--tcoh 3] [--lam SU2] ...
"""
import argparse, json, math, sys, time
import numpy as np
import ek_common as C

QLO, QHI, NQLOG = 3.0, 3.0e7, 700
LEGS = (0.25, 0.35, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0)
T_GLOBAL, BUDGET = [time.time()], [1e9]          # wall-clock budget per process (ek_scan sets it)


class Incomplete(Exception):
    pass


def coh_grid(s, A0, a0, dk_m=0.05, extra=3.0):
    kmax = (1.3 * A0 + 4.0 * s['Mh'] + extra) * s['m']
    nq = int(math.ceil(kmax / (dk_m * s['m'])))
    return C.grid_linear(a0 * kmax, nq)


def amp(rho, m, f):
    return math.sqrt(2 * max(rho, 0)) / m / (2 * f)


def tables(spec):
    return [C.load_ptable(round(s['Mh'], 5)) for s in spec]


# ------------------------------------------------------------------ constant-amplitude regression
def check():
    out = dict(label='C', note='constant amplitude, no expansion, no collisions: cumulative energy fraction F after n '
               'half-periods, coherent (exact modes) and dephased (tables), against the mf_hybrid skeptic\'s exact '
               'multi-period runs (sk_axial.json)')
    sk = json.load(open(C.TRK / 'mf_hybrid' / 'skeptic' / 'results' / 'sk_axial.json'))['multi_period']
    spec, mpar = C.species('HFa')
    sm = C.sm_coeffs()
    mP = C.MP_GEV / C.H_REF
    Aset = {0: 12.34, 1: 42.3}
    rhoP = [(2 * spec[0]['f'] * Aset[0] * mpar[0]) ** 2 / 2, (2 * spec[1]['f'] * Aset[1] * mpar[1]) ** 2 / 2]
    tabs = tables(spec)
    res = {}
    for mode in ('coherent', 'dephased'):
        grids, occ = [], []
        for j, s in enumerate(spec):
            A0 = Aset[s['par']]
            if mode == 'coherent':
                q, w = coh_grid(s, A0, 1.0, dk_m=0.04)
            else:
                q, w = C.grid_linear((1.3 * A0 + 4 * s['Mh'] + 3) * s['m'], 1500)
            grids.append((q, w)); occ.append(np.zeros(2 * len(q)))
        for nh in (2, 4, 8, 16):
            st = dict(mP=mP, t=0.0, t_end=nh * math.pi / mpar[1], a=1.0, rho_sig=rhoP[0], rho_chi=rhoP[1], rho_h=0.0, n_h=0.0, rho_th=0.0)
            inp = C.SP / ('chk_%s_%d.in' % (mode, nh))
            C.write_input(inp, st, spec, mpar, sm, dict(coherent=int(mode == 'coherent'), noexp=1, fixamp=1, nocoll=1),
                          grids, occ, tables=tabs if mode == 'dephased' else None)
            t1 = time.time()
            rec = C.run_engine(inp, C.SP / ('chk_%s_%d' % (mode, nh)))
            stt = C.read_state(C.SP / ('chk_%s_%d' % (mode, nh)), [len(g[0]) for g in grids])
            Fs = {}
            for j, s in enumerate(spec):
                q, w = grids[j]
                om = np.sqrt(q ** 2 + s['M'] ** 2)
                E = float(np.sum(np.r_[w, w] * np.r_[om, om] * stt['occ'][j]))
                Fs['N%d' % (j + 1)] = E / rhoP[s['par']]
            res.setdefault(mode, {})[nh] = dict(F=Fs, seconds=time.time() - t1)
            print(mode, nh, {k: '%.4f' % v for k, v in Fs.items()}, '%.1fs' % (time.time() - t1), flush=True)
    # the skeptic's numbers (sigma half-periods counted in sigma periods: compare at equal chi time for chi channels)
    skn = {}
    for key, nm in (('sigma->N1', 'N1'), ('chi->N2', 'N2'), ('chi->N3', 'N3'), ('chi->N4', 'N4')):
        rows = sk[key]['rows']
        skn[nm] = {r['nhalf']: r['F_tot'] for r in rows}
    out['skeptic_F_tot'] = skn
    out['this_track'] = res
    out['skeptic_F_sat'] = {'N1': 0.14, 'N2': 0.13, 'N3': 0.085, 'N4': 0.055}
    C.dump(C.RES / 'check_constant_amplitude.json', out)


# ------------------------------------------------------------------ full history
def start_state(name, H_GeV, t0=0.19, state=None):
    mP = C.MP_GEV / H_GeV
    if state is not None:
        return dict(state, mP=mP)
    b = C.bench_state(name, t0)
    rho = 3 * mP ** 2 * b['H'] ** 2
    f_chi = b['f_chi'] + b['f_rad']            # the track's 0.9% drain returned to chi (production restarts exactly)
    return dict(mP=mP, t=t0, a=b['a'], H=b['H'], rho_sig=b['f_sig'] * rho, rho_chi=f_chi * rho, rho_h=0.0, n_h=0.0,
                rho_th=0.0, bench=b)


def run_history(name, tag, H_GeV=C.H_REF, mnu=C.MNU_REF, t0=0.19, t_coh=3.0, t_end=1e8, lam_group='SU2', G12=4.3e-4,
                KLscale=1.0, cBBscale=1.0, dk_m=0.05, winv=1.0, state=None, occ0=None, verbose=True, cdt=0.04, dshell=0.0,
                nqlog=None):
    global NQLOG
    if nqlog:
        NQLOG = nqlog
    t_start = time.time()
    spec, mpar = C.species(name, H_GeV, mnu)
    sm = C.sm_coeffs(G12=G12, lam_group=lam_group)
    sm = dict(sm, cBB=sm['cBB'] * cBBscale)
    st = start_state(name, H_GeV, t0, state)
    tabs = tables(spec)
    info = dict(name=name, tag=tag, H_GeV=H_GeV, mnu=mnu, t0=st['t'], t_coh=t_coh, lam_group=lam_group, lam_eff=sm['lam_eff'],
                cBB=sm['cBB'], kH=sm['kH'], kL=sm['kL'], G12=G12, KLscale=KLscale, dk_m=dk_m, winv=winv, dshell=dshell, nqlog=NQLOG,
                cdt=cdt, cBBscale=cBBscale,
                species=[{k: s[k] for k in ('par', 'M', 'f', 'm', 'h2', 'Grest', 'kres_m', 'Gpert', 'mnu')} for s in spec])
    recs = []
    # ---- coherent stage
    if t_coh > st['t']:
        grids, occ = [], []
        for j, s in enumerate(spec):
            rhoP = st['rho_sig'] if s['par'] == 0 else st['rho_chi']
            A0 = amp(rhoP, s['m'], s['f'])
            q, w = coh_grid(s, A0, st['a'], dk_m)
            grids.append((q, w))
            occ.append(np.zeros(2 * len(q)) if occ0 is None else C.rebin(occ0[j][0], occ0[j][1], occ0[j][2], np.r_[0.0, 0.5 * (q[1:] + q[:-1]), q[-1] + 0.5 * (q[1] - q[0])]))
        info['coherent_grid'] = [dict(nq=len(g[0]), qmax=float(g[0][-1])) for g in grids]
        # legs with restart files (event times + Bloch vectors): a run can be continued in a later call (--resume)
        legs = [x for x in LEGS if st['t'] < x < t_coh] + [t_coh]
        ck = C.SP / ('ckpt_%s.npz' % tag)
        nqs = [len(g[0]) for g in grids]
        cur = dict(st); bloch = None; frozen = None; ileg = 0; crec = []; last = None
        if ck.exists():
            z = np.load(ck, allow_pickle=True)
            ileg = int(z['ileg']); cur = z['cur'].item(); bloch = list(z['bloch']); frozen = list(z['frozen']); occ = list(z['occ'])
            crec = list(z['crec']); last = z['last'].item()
            if verbose:
                print(tag, 'resuming coherent stage at leg %d (t = %.4g)' % (ileg, cur['t']), flush=True)
        t1 = time.time()
        while ileg < len(legs):
            inp = C.SP / ('%s_coh.in' % tag)
            C.write_input(inp, dict(cur, t_end=legs[ileg]), spec, mpar, sm, dict(coherent=1), grids, occ, bloch=bloch, frozen=frozen,
                          KLscale=KLscale, winv=winv, cdt=cdt)
            crec.append(C.run_engine(inp, C.SP / ('%s_coh' % tag), timeout=3600))
            s1 = C.read_state(C.SP / ('%s_coh' % tag), nqs)
            ev, bloch = C.read_bloch(C.SP / ('%s_coh' % tag), nqs)
            occ, frozen = s1['occ'], s1['frozen']
            cur = dict(mP=st['mP'], t=s1['t'], a=s1['a'], rho_sig=s1['rho_sig'], rho_chi=s1['rho_chi'], rho_h=s1['rho_h'], n_h=s1['n_h'],
                       rho_th=s1['rho_th'], events=ev)
            ileg += 1
            last = {k: s1[k] for k in ('t', 'a', 'H', 'rho_sig', 'rho_chi', 'rho_h', 'n_h', 'rho_th', 'rho_tot')}
            np.savez(ck, ileg=ileg, cur=np.array(cur, dtype=object), last=np.array(last, dtype=object), bloch=np.array(bloch, dtype=object),
                     frozen=np.array(frozen, dtype=object), occ=np.array(occ, dtype=object), crec=np.array(crec, dtype=object))
            if verbose:
                print(tag, 'coherent leg to t = %.4g done (%.0fs)' % (cur['t'], time.time() - t1), flush=True)
            if ileg < len(legs) and time.time() - T_GLOBAL[0] > BUDGET[0]:
                raise Incomplete(tag)
        # the first record of every later leg repeats the previous leg's last row without collision diagnostics: drop it
        recs.extend([crec[0]] + [{k: v[1:] for k, v in r.items()} for r in crec[1:]])
        s1 = dict(last, occ=occ)
        info['coherent_seconds'] = time.time() - t1
        if verbose:
            print(tag, 'coherent stage done %.0fs' % (time.time() - t1), flush=True)
        # re-bin to the log grid
        edges = np.geomspace(QLO, QHI, NQLOG + 1)
        ql, wl = C.grid_log(QLO, QHI, NQLOG)
        occ2 = []
        for j in range(4):
            q, w = grids[j]
            occ2.append(C.rebin(q, w, s1['occ'][j], edges))
        st2 = dict(mP=st['mP'], t=s1['t'], a=s1['a'], rho_sig=s1['rho_sig'], rho_chi=s1['rho_chi'], rho_h=s1['rho_h'], n_h=s1['n_h'],
                   rho_th=s1['rho_th'])
        # number and energy check of the re-binning
        chk = []
        for j, s in enumerate(spec):
            q, w = grids[j]
            omc = np.sqrt((q / s1['a']) ** 2 + s['M'] ** 2); oml = np.sqrt((ql / s1['a']) ** 2 + s['M'] ** 2)
            chk.append(dict(N_ratio=float(np.sum(np.r_[wl, wl] * occ2[j]) / max(np.sum(np.r_[w, w] * s1['occ'][j]), 1e-300)),
                            E_ratio=float(np.sum(np.r_[wl * oml, wl * oml] * occ2[j]) / max(np.sum(np.r_[w * omc, w * omc] * s1['occ'][j]), 1e-300))))
        info['rebin_check'] = chk
        info['state_tcoh'] = {k: s1[k] for k in ('t', 'a', 'H', 'rho_sig', 'rho_chi', 'rho_h', 'n_h', 'rho_th', 'rho_tot')}
    else:
        ql, wl = C.grid_log(QLO, QHI, NQLOG)
        occ2 = [np.zeros(2 * NQLOG) for _ in range(4)] if occ0 is None else [C.rebin(o[0], o[1], o[2], np.geomspace(QLO, QHI, NQLOG + 1)) for o in occ0]
        st2 = dict(st)
    # ---- dephased stage
    inp = C.SP / ('%s_dep.in' % tag)
    C.write_input(inp, dict(st2, t_end=t_end), spec, mpar, sm, dict(coherent=0), [(ql, wl)] * 4, occ2, tables=tabs,
                  KLscale=KLscale, winv=winv, dshell=dshell)
    t1 = time.time()
    rd = C.run_engine(inp, C.SP / ('%s_dep' % tag), timeout=3600)
    recs.append({k: v[1:] for k, v in rd.items()} if recs else rd)
    info['dephased_seconds'] = time.time() - t1
    rec = {k: np.concatenate([r[k] for r in recs]) for k in recs[0]}
    o = np.argsort(rec['t'], kind='stable')
    rec = {k: v[o] for k, v in rec.items()}
    summ = summarise(rec, st['mP'], H_GeV, spec)
    info['seconds'] = time.time() - t_start
    ckf = C.SP / ('ckpt_%s.npz' % tag)
    if ckf.exists():
        ckf.unlink()
    if verbose:
        print(tag, json.dumps({k: summ[k] for k in ('K', 't_SM_half', 't_rad_half', 't_thermal', 'T_RH_GeV', 'T_max_GeV')}, default=float),
              '%.0fs' % info['seconds'], flush=True)
    return dict(info=info, summary=summ, rec=rec)


def first(t, y, thr):
    i = np.where(y >= thr)[0]
    return float(t[i[0]]) if len(i) else None


def summarise(rec, mP, H_GeV, spec):
    t = rec['t']
    K = rec['a'] ** 2 * rec['H']
    f_SM = rec['f_hard'] + rec['f_th']
    f_rad = f_SM + rec['f_N']
    therm_frac = np.where(f_SM > 0, rec['f_th'] / np.maximum(f_SM, 1e-300), 0.0)
    Tsm = rec['T_SM']
    # "thermal equilibrium of the SM bath": thermal part >= 90% of the SM energy (KL relaxation) and SM >= 50% of rho
    i_th = np.where((therm_frac >= 0.9) & (f_SM >= 0.5))[0]
    t_th = float(t[i_th[0]]) if len(i_th) else None
    i_th1 = np.where(therm_frac >= 0.9)[0]
    t_th_bath = float(t[i_th1[0]]) if len(i_th1) else None
    tsmh = first(t, f_SM, 0.5)
    TRH = float(np.interp(tsmh, t, Tsm)) * H_GeV if tsmh else None
    Tmax_i = int(np.argmax(rec['T_th']))
    # Gamma_th(t): energy-weighted removal rate of the N's (all species) and per species band rates
    fN = np.array([rec['fN%d' % j] for j in (1, 2, 3, 4)])
    Gall = np.array([rec['Gall%d' % j] for j in (1, 2, 3, 4)])
    Gth = np.sum(fN * Gall, axis=0) / np.maximum(np.sum(fN, axis=0), 1e-300)
    sel = lambda x: float(np.interp(x, t, Gth))
    out = dict(K=float(K[-1]), t_end=float(t[-1]), condensate_left=float(rec['f_sig'][-1] + rec['f_chi'][-1]),
               t_SM_1em3=first(t, f_SM, 1e-3), t_SM_half=tsmh, t_rad_half=first(t, f_rad, 0.5), t_rad_99=first(t, f_rad, 0.99),
               t_SM_99=first(t, f_SM, 0.99), t_thermal=t_th, t_bath_thermal=t_th_bath,
               T_RH_GeV=TRH, T_RH_definition='T of the SM bath (temperature of its energy) when rho_SM = rho_tot/2',
               T_max_GeV=float(rec['T_th'][Tmax_i] * H_GeV), t_Tmax=float(t[Tmax_i]),
               T_K_GeV=float((30 * 3 * K[-1] ** 2 * mP ** 2 / (math.pi ** 2 * C.GSTAR)) ** 0.25 * H_GeV),
               T_at_thermal_GeV=float(np.interp(t_th, t, Tsm) * H_GeV) if t_th else None,
               Gamma_th_at={'%g' % x: sel(x) for x in (0.2, 0.3, 0.5, 1, 3, 10, 30, 100, 1000)},
               f_rad_at={'%g' % x: float(np.interp(x, t, f_rad)) for x in (0.3, 0.45, 0.6, 1, 3, 10, 100, 1000)},
               f_SM_at={'%g' % x: float(np.interp(x, t, f_SM)) for x in (0.3, 0.45, 0.6, 1, 3, 10, 100, 1000)},
               K_at={'%g' % x: float(np.interp(x, t, K)) for x in (1, 3, 10, 100, 1000, 1e4)})
    return out


def export_fluid(res, name, path, t_switch_note):
    """fluid_<name>.json: bench-format history rows (t >= t0) + the kinetic columns; see NOTES.md section 6."""
    rec, info, summ = res['rec'], res['info'], res['summary']
    d = C.load_bench(name)
    pr = d['params']
    mP = C.MP_GEV / info['H_GeV']
    t, a, H = rec['t'], rec['a'], rec['H']
    # thin to <= 1500 rows, log-spaced in t
    tt = np.unique(np.r_[t[0], np.geomspace(max(t[0], 1e-3), t[-1], 1400), t[-1]])
    idx = np.unique(np.clip(np.searchsorted(t, tt), 0, len(t) - 1))
    g = lambda k: rec[k][idx]
    rho = 3 * mP ** 2 * g('H') ** 2
    f_sig, f_chi = g('f_sig'), g('f_chi')
    f_SM = g('f_hard') + g('f_th')
    f_rad = f_SM + g('f_N')
    # cycle-mean R with w = 0 condensates and w = 1/3 radiation (N relativistic in the band; NR N tail neglected)
    R = (f_sig + f_chi) * rho / mP ** 2
    m_chi, m_sig = pr['m_chi_over_Hi'], pr['m_sig_over_Hi']
    if pr['smech'] == 'a':
        g_s = 3.0 / (pr['v_MP'] * mP)                     # m_s = g_s sqrt(v^2 + <dchi^2>), mf_chain/mf_hybrid convention
        m_s = g_s * np.sqrt((pr['v_MP'] * mP) ** 2 + f_chi * rho / m_chi ** 2)
    else:
        m_s = np.sqrt(np.maximum(9.0 - 0.75 * R, 0.0))
    N = np.log(g('a'))
    eps = 1.5 * (1 + (f_rad / 3.0))
    hist = dict(t=g('t'), N=N, a=g('a'), H=g('H'), eps_H=eps, R=R, R_osc_sig=3 * f_sig * rho / mP ** 2, R_osc_chi=3 * f_chi * rho / mP ** 2,
                m_s=m_s, f_sig=f_sig, f_chi=f_chi, f_rad=f_rad, f_SM=f_SM, f_N=g('f_N'), f_SM_hard=g('f_hard'), f_SM_thermal=g('f_th'),
                T_SM=g('T_SM'), T_thermal=g('T_th'), Tcal2=g('Tcal2'), m_H2_medium=C.sm_coeffs()['kH'] * g('Tcal2'),
                tau_KL=np.where(g('tau_KL') > 1e299, np.nan, g('tau_KL')), K_running=g('a') ** 2 * g('H'),
                Gamma_th=np.sum([g('fN%d' % j) * g('Gall%d' % j) for j in (1, 2, 3, 4)], axis=0) / np.maximum(np.sum([g('fN%d' % j) for j in (1, 2, 3, 4)], axis=0), 1e-300),
                **{'Gamma_band_N%d' % j: g('Gband%d' % j) for j in (1, 2, 3, 4)},
                **{'nbar_band_N%d' % j: g('nband%d' % j) for j in (1, 2, 3, 4)},
                **{'f_N%d' % j: g('fN%d' % j) for j in (1, 2, 3, 4)},
                **{'A_N%d' % j: g('A%d' % j) for j in (1, 2, 3, 4)})
    cols = dict(t='H_i^-1 from the slice', N='ln a', a='scale factor (a = 1 at sigma = sigma_c)', H='H/H_i',
                eps_H='-Hdot/H^2 = (3/2)(1 + w), w = f_rad/3 (cycle mean)', R='cycle-mean Ricci scalar / H_i^2 = (rho_m)/M_P^2',
                R_osc_sig='amplitude of the 2 m_sig oscillation of R (3 rho_sig/M_P^2)', R_osc_chi='amplitude at 2 m_chi',
                m_s='spectator mass / H_i (HFa: g_s sqrt(v^2 + rho_chi/m_chi^2); HFb: sqrt(9 - 3R/4))',
                f_sig='sigma condensate energy fraction', f_chi='chi condensate', f_rad='N + SM', f_SM='SM (hard + thermal)',
                f_N='singlet N (all species)', f_SM_hard='SM not yet thermalised', f_SM_thermal='SM thermalised (Kurkela-Lu relaxation)',
                T_SM='temperature of the SM energy (H_i units)', T_thermal='temperature of the thermalised part', Tcal2='medium scale^2 (H_i^2)',
                m_H2_medium='Higgs asymptotic medium mass^2 kH Tcal^2 (H_i^2)', tau_KL='hard-SM relaxation time (H_i^-1)',
                K_running='a^2 H', Gamma_th='N-removal rate, energy weighted over all N (H_i): decay + in-medium',
                Gamma_band_N1='occupation-weighted removal rate of in-band N1 (H_i)', nbar_band_N1='mean occupation in the band k < mu5,max',
                A_N1='mu5,max/m of the parent for N1')
    out = dict(name='%s_kinetic' % name, label='C on the stated rates (model choices throughout); e2e_kinetic',
               note=t_switch_note, H_i_GeV=info['H_GeV'], slice='a = 1, H = H_i at sigma = sigma_c (as the mf_hybrid bench)',
               params=dict(pr, H_i_GeV=info['H_GeV'], s_star_GeV=pr['s_star_MP'] * C.MP_GEV), info=info, summary=summ,
               spectator=dict(Gamma_s_assumed=d['spectator'].get('Gamma_s_assumed', 1.3e-6)),
               history_columns=cols,
               history={k: [None if not np.isfinite(x) else float('%.7g' % x) for x in np.asarray(v, float)] for k, v in hist.items()})
    C.dump(path, out)
    return out


def jsonable(res, keep_rec=False):
    o = dict(info=res['info'], summary=res['summary'])
    if keep_rec:
        r = res['rec']
        idx = np.unique(np.linspace(0, len(r['t']) - 1, 300).astype(int))
        o['rec'] = {k: [float('%.6g' % v[i]) for i in idx] for k, v in r.items()}
    return o


if __name__ == '__main__':
    if sys.argv[1] == 'check':
        check()
        sys.exit()
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('name'); ap.add_argument('tag')
    ap.add_argument('--H', type=float, default=C.H_REF); ap.add_argument('--mnu', type=float, default=0.05)
    ap.add_argument('--tcoh', type=float, default=3.0); ap.add_argument('--t0', type=float, default=0.19)
    ap.add_argument('--tend', type=float, default=1e8); ap.add_argument('--lam', default='SU2')
    ap.add_argument('--KL', type=float, default=1.0); ap.add_argument('--cBB', type=float, default=1.0)
    ap.add_argument('--dk', type=float, default=0.05); ap.add_argument('--cdt', type=float, default=0.04)
    ap.add_argument('--dshell', type=float, default=0.0); ap.add_argument('--nqlog', type=int, default=0)
    ap.add_argument('--state', default=None, help='lattice hand-over JSON (NOTES.md section 6): t, a, H, rho_sig, rho_chi, '
                    'rho_h, n_h, rho_th and per-species E_N, n_N lists (physical, H_i units); t0 is then the state t')
    ap.add_argument('--out', default=None, help='export: path of the fluid json (default ../results/fluid_<name>.json)')
    a = ap.parse_args()
    state, occ0, t0 = None, None, a.t0
    if a.state:
        import ek_scan as SC
        z = json.load(open(a.state))
        spec_, _ = C.species(a.name, a.H, a.mnu * 1e-9)
        state = {k: float(z[k]) for k in ('t', 'a', 'H', 'rho_sig', 'rho_chi', 'rho_h', 'n_h', 'rho_th')}
        occ0 = []
        for j, s in enumerate(spec_):
            q, w = C.grid_linear(state['a'] * 3e4, 6000)
            occ0.append((q, w, SC.occ_from_moments(float(z['E_N'][j]), float(z['n_N'][j]), q, w, state['a'], s['M'])))
        t0 = state['t']
    r = run_history(a.name, a.tag, a.H, a.mnu * 1e-9, t0, max(a.tcoh, t0 + 2.0) if a.state else a.tcoh, a.tend, a.lam, KLscale=a.KL,
                    cBBscale=a.cBB, dk_m=a.dk, cdt=a.cdt, dshell=a.dshell, nqlog=a.nqlog or None, state=state, occ0=occ0)
    C.dump(C.SP / ('run_%s.json' % a.tag), jsonable(r, keep_rec=True))
    if a.cmd == 'export':
        note = ('kinetic history from t0 = %.4g H_i^-1 (%s); coherent mode stage to t_coh = %.4g, dephased rate equations after. '
                'Lattice splice: see NOTES.md section 6 (switch at t_sw in [0.45, 3] H_i^-1 by moments, or splice this file at '
                't >= t_sw when the lattice fractions agree with the state_tcoh/history rows).' %
                (t0, 'lattice hand-over state' if a.state else 'homogeneous bench state at the first chi turning point', r['info']['t_coh']))
        export_fluid(r, a.name, a.out or (C.RES / ('fluid_%s.json' % a.name)), note)
