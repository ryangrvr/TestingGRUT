"""e2e_lattice job registry and driver.  Each job is one lattice run; `python3 -B jobs.py TAG` creates it (first call),
runs/resumes it within ~560 s of wall time, and when it is done writes the run summary results/runs/<TAG>.json
(time series thinned, percolation, spectra, GW).  Every configuration is a model choice.

Base: HFa, design H = 1.30e11 GeV, 64^3, L = 0.42 H^-1 (comoving at the slice), dt = dx/25, k_cut = 1.2 m_chi,
seed 1, start t0 = -0.02 (sigma just above sigma_c), Pauli closure with local blocking, Gamma_th = 0, tails free.
"""
import json, math, sys, glob
import numpy as np
import latcommon as LC
import latpost as LP

GWS = (0.3, 0.45, 0.5, 0.6, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0)
BASE = dict(name='HFa', H_GeV=1.30e11, closure='pauli', Gth=0.0, N=64, L=0.42, kcut_mchi=1.2, seed=1, t_end=3.0,
            Bglobal=0, gw=1, dtfac=25.0)
JOBS = {
    # closures (HFa, design H)
    'P0':      dict(),
    'P0G':     dict(Bglobal=1, gw=0),
    'P0GM':    dict(Bglobal=2),
    'P10':     dict(Gth=10.0),
    'P10GM':   dict(Gth=10.0, Bglobal=2, gw=0),
    'TR':      dict(closure='track'),
    'HU':      dict(closure='hubble'),
    # HFb
    'B_P0':    dict(name='HFb', gw=0),
    'B_P0GM':  dict(name='HFb', Bglobal=2, gw=0),
    # strong-holder self-consistent scale (mf_critic: c0 = 1.40 at 8.40e10, Gamma_th = 0; c0 = 1.25-1.30 at 5.0e10, Gamma_th = 10)
    'H84_P0':  dict(H_GeV=8.40e10),
    'H84_P0GM': dict(H_GeV=8.40e10, Bglobal=2, gw=0),
    'H50_P10': dict(H_GeV=5.0e10, Gth=10.0, gw=0),
    'H50_P10GM': dict(H_GeV=5.0e10, Gth=10.0, Bglobal=2, gw=0),
    'H240_HU': dict(H_GeV=2.40e11, closure='hubble', t_end=1.5, gw=0),
    # convergence of the base (P0)
    'C_dt2':   dict(dtfac=50.0, t_end=1.5, gw=0),
    'C_L28':   dict(L=0.28, t_end=1.5),
    'C_L63':   dict(L=0.63, t_end=1.5),
    'C_N96':   dict(N=96, t_end=1.5),
    'C_N96L63': dict(N=96, L=0.63, t_end=1.5),
    'C_kc08':  dict(kcut_mchi=0.8, t_end=1.5),
    'C_kc16':  dict(kcut_mchi=1.6, t_end=1.5),
    'C_s2':    dict(seed=2, t_end=1.5, gw=0),
    'C_s3':    dict(seed=3, t_end=1.5, gw=0),
    'HU_s2':   dict(closure='hubble', seed=2, t_end=1.5, gw=0),
    'HU_N96':  dict(closure='hubble', N=96, t_end=1.5, gw=0),
    'HU_L63':  dict(closure='hubble', L=0.63, t_end=1.5, gw=0),
    # long runs
    'P0_long': dict(t_end=8.0),
    'P0GM_long': dict(Bglobal=2, t_end=8.0, gw=0),
}


def cfg(tag):
    c = dict(BASE)
    c.update(JOBS[tag])
    return c


def create(tag):
    c = cfg(tag)
    dx = c['L'] / c['N']
    gws = [g for g in GWS if g <= c['t_end'] + 1e-9] if c['gw'] else []
    snaps = sorted(set([round(x, 4) for x in list(np.arange(0.05, 1.0001, 0.025)) +
                        list(np.arange(1.1, min(c['t_end'], 3.0) + 1e-9, 0.1)) +
                        list(np.arange(3.25, c['t_end'] + 1e-9, 0.25))] + list(gws)))
    LC.make_run(tag, c['name'], c['H_GeV'], c['closure'], c['Gth'], c['N'], c['L'], dt=dx / c['dtfac'],
                kcut_mchi=c['kcut_mchi'], seed=c['seed'], t0=-0.02, t_end=c['t_end'], snaps=snaps, gw_snaps=gws,
                Bglobal=c['Bglobal'], gw=c['gw'])


def summarize(tag, keep_snaps=False):
    c = cfg(tag)
    rd = LC.RUNS / tag
    info = json.load(open(rd / 'info.json'))
    P = info['P']
    q = (P['v'] / P['mP']) ** 2
    D = LC.load_diag(tag)
    t = D['t']
    frad = D['fN1'] + D['fN2'] + D['fN3'] + D['fN4'] + D['fSM']
    # thinned series (every ~0.005 early, ~0.02 late)
    sel = np.unique(np.concatenate([np.searchsorted(t, np.arange(t[0], 1.0, 0.004)),
                                    np.searchsorted(t, np.arange(1.0, t[-1] + 1e-9, 0.02)), [len(t) - 1]]))
    sel = sel[sel < len(t)]
    ser = {k: [float('%.7g' % x) for x in D[k][sel]] for k in LC.COLS}
    ser['frad'] = [float('%.7g' % x) for x in frad[sel]]
    # cycle-averaged R and w over 0.05 windows
    def runmean(y, w=0.05):
        out = np.empty_like(y)
        for i in range(len(y)):
            m = (t >= t[i] - w / 2) & (t <= t[i] + w / 2)
            out[i] = y[m].mean()
        return out
    Rbar = runmean(D['R_mean'])
    wbar = runmean(D['w'])
    ser['R_cyc'] = [float('%.7g' % x) for x in Rbar[sel]]
    ser['w_cyc'] = [float('%.7g' % x) for x in wbar[sel]]
    # snapshots
    snaps = []
    for f in sorted(glob.glob(str(rd / 'snap_*.bin'))):
        S = LC.load_snap(f)
        e = dict(t=S['t'], a=S['a'])
        e['neg'] = LP.clusters(S['chi'] < 0)
        e['neg_av'] = LP.clusters(S['chiav'] < 0)
        k, Pk, _ = LP.power(S['chi'], S['L'])
        e['chi_kpeak'] = float(k[np.argmax(Pk)])
        e['chi_var'] = float(S['chi'].var())
        e['R_rms'] = float(S['R'].std()); e['R_mean'] = float(S['R'].mean())
        e['sa_mean'] = float(S['sa'].mean()); e['sa_std'] = float(S['sa'].std())
        e['sb_mean'] = float(S['sb'].mean()); e['sb_std'] = float(S['sb'].std())
        e['ms2a_std'] = float((9.0 * S['chi'].astype(np.float64) ** 2).std())
        if 'hdot' in S:
            kg, sp, tot, raw = LP.gw_spectrum(S['hdot'], S['L'], q)
            e['gw'] = dict(k=[float('%.6g' % x) for x in kg], dlnk=[float('%.5g' % x) for x in sp], rho_tt=tot, rho_raw=raw,
                           kpeak=float(kg[np.argmax(sp)]))
            e['chi_spec'] = dict(k=[float('%.6g' % x) for x in k], P=[float('%.5g' % x) for x in Pk])
        snaps.append(e)
    # headline numbers
    i_coh = int(np.argmax(D['fcoh'] < 0.1)) if (D['fcoh'] < 0.1).any() else None
    head = dict(t_end=float(t[-1]), a_end=float(D['a'][-1]), H_end=float(D['H'][-1]),
                t_eps1=float(t[np.argmax(D['eps'] >= 1)]) if (D['eps'] >= 1).any() else None,
                t_fcoh_lt_0p1=float(t[i_coh]) if i_coh is not None else None,
                t_fcoh_lt_0p01=float(t[np.argmax(D['fcoh'] < 0.01)]) if (D['fcoh'] < 0.01).any() else None,
                fminus_max=float(D['fminus'].max()), t_fminus_max=float(t[np.argmax(D['fminus'])]),
                fminus_av_max=float(D['fminus_av'].max()), area_av_max=float(D['area_av'].max()),
                area_inst_max=float(D['area_inst'].max()),
                neg_wraps_ever=bool(any(s['neg']['wraps'] for s in snaps)),
                neg_av_wraps_ever=bool(any(s['neg_av']['wraps'] for s in snaps)),
                neg_largest_max=float(max([s['neg']['largest'] for s in snaps] + [0.0])),
                fgrad_max=float(D['fgrad'].max()), t_fgrad_max=float(t[np.argmax(D['fgrad'])]),
                frad_end=float(frad[-1]), fgrad_end=float(D['fgrad'][-1]), fchi_end=float(D['fchi'][-1]),
                fsig_end=float(D['fsig'][-1]), R_rms_max=float(D['R_rms'].max()),
                budget_end=float(D['budget'][-1]), budget_absmax=float(np.abs(D['budget']).max()))
    for tt in (0.3, 0.45, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0):
        if tt <= t[-1] + 0.01:
            tt_ = tt
            tt = min(tt, float(t[-1]))
            head['frad_t%g' % tt_] = float(np.interp(tt, t, frad))
            head['fgrad_t%g' % tt_] = float(np.interp(tt, t, D['fgrad']))
            head['w_cyc_t%g' % tt_] = float(np.interp(tt, t, wbar))
            head['R_cyc_t%g' % tt_] = float(np.interp(tt, t, Rbar))
            head['R_rms_t%g' % tt_] = float(np.interp(tt, t, D['R_rms']))
            head['sa_std_t%g' % tt_] = float(np.interp(tt, t, D['sa_std']))
            head['sb_std_t%g' % tt_] = float(np.interp(tt, t, D['sb_std']))
            head['ms2a_std_over_mean_t%g' % tt_] = float(np.interp(tt, t, D['ms2a_std'] / np.maximum(D['ms2a_mean'], 1e-30)))
            head['ms2b_std_t%g' % tt_] = float(np.interp(tt, t, D['ms2b_std']))
            head['chi_std_t%g' % tt_] = float(np.interp(tt, t, D['chi_std']))
    out = dict(tag=tag, config=c, info={k: info[k] for k in info if k not in ('P',)}, P=P, head=head, series=ser,
               snaps=snaps)
    LC.dump(LC.RES / 'runs' / ('%s.json' % tag), out)
    return out


if __name__ == '__main__':
    tag = sys.argv[1]
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 560.0
    rd = LC.RUNS / tag
    if not (rd / 'run.par').exists():
        create(tag)
    done = LC.run(tag, budget=budget)
    if done and not (LC.RES / 'runs' / ('%s.json' % tag)).exists():
        o = summarize(tag)
        print(json.dumps(o['head'], indent=0))
    print('DONE' if done else 'NOT DONE')
