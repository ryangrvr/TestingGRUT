#!/usr/bin/env python3
"""Track inf_higgs driver: Higgs survival after inflation on the inflaton benchmarks (own code).

Usage (one thread each, foreground):
  python3 -B ih_run.py bench              # xi_EW = 0.1 on every benchmark history -> results/bench_D.json
  python3 -B ih_run.py xiscan <name> <hist> [iH]   # lower edge (xi_EW 0.02 ... 0.12)   -> results/xiscan_*.json
  python3 -B ih_run.py upper  <name> <hist>        # large xi (upper edge)              -> results/upper_*.json
  python3 -B ih_run.py variants           # Hartree at the critical M_t, portal / A-term variants -> results/variants.json
Histories ('hist', every one a model choice):
  MD    Gamma_I = 0 over the computed window (tau <= 80): the realistic slow decay of inf_reheat (Gamma/m_I <~ 1e-6) is a
        matter era here; homogeneous oscillation of the benchmark's own potential.
  MDf   as MD, but the inflaton fragments into a w = 0 fluid at tau_end + 0.43 (T1) / 0.35 (E1) (inf_potential [Q]).
  G05   Gamma_I = 0.5 H_i (inf_potential's / the engine's convention; kinematically blocked per the inf_potential skeptic).
Mode runs are cached in the scratchpad (npz) so the danger function can be re-evaluated (M_t grid, plasma masses).
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np

import ih_common as C
import ih_modes as M

H0C_LN = math.log(0.05 / 2.248e-4)       # ln(k_*/a0H0), h = 0.674 [R]
MT_PDG = (172.28, 172.57, 172.86)
MT_GRID = sorted(set(MT_PDG + tuple(np.round(np.arange(172.70, 173.46, 0.05), 3))))
MH_VARIANTS = (125.09, 125.13, 125.31)
MT_COARSE = (172.75, 172.85, 172.95, 173.05, 173.15, 173.25, 173.35)
T_RH_SLOW = 1e11                          # GeV, slow-decay reference for N_obs (y ~ y_nat; model choice)


def N_obs(name, hist, iH=0):
    """e-folds from the exit of today's horizon to S1 (a = 1)"""
    if name == 'M1r':
        return 60.0
    B = C.load_bench(name, iH)
    r = B['rec']
    if hist == 'G05':
        Ns = r['Nstar']
    else:
        inst = [x for x in r['Nstar_alternatives'] if x['history'] == ['inst']][0]['Nstar']
        Ns = inst - math.log(r['T_inst_GeV'] / T_RH_SLOW) / 3
    return round(Ns + H0C_LN - 1.0, 2)


def background(name, hist, iH=0, tau_end=80.0, Nb0=5.0):
    gam = 0.5 if hist == 'G05' else 0.0
    if name == 'M1r':
        return C.Background(None, gam=0.5, b=3.0, tau_end=40.0)
    B = C.load_bench(name, iH)
    frag = None
    if hist == 'MDf':
        b0 = C.Background(B, gam=gam, Nb0=Nb0, tau_end=3.0)
        frag = b0.tau_endinf + C.TFRAG_AFTER_END[name]
    return C.Background(B, gam=gam, Nb0=Nb0, tau_end=tau_end, frag=frag)


def run_cached(name, hist, xi, iH=0, tag='', Nb0=5.0, n=260, Kmax=None, lamIH=0.0, Aterm=0.0, sm_mt=172.57, rerun=False,
               tau_end=80.0, bg=None):
    key = f'{name}_{hist}_iH{iH}_xi{xi:.6g}_Nb{Nb0}_n{n}_K{Kmax}_lam{lamIH:.3g}_A{Aterm:.3g}_mt{sm_mt}_te{tau_end}{tag}'
    path = C.SP / 'runs' / (key + '.npz')
    path.parent.mkdir(exist_ok=True)
    if path.exists() and not rerun:
        z = np.load(path, allow_pickle=True)
        R = dict(res=json.loads(str(z['res'])), rows=z['rows'], spec=z['spec'], ks=z['ks'], wq=z['wq'], Rb=z['Rb'],
                 m2xbar=z['m2xbar'], key=key)
        if 'rbv' not in z.files:          # Rbar from the first post-inflation period on (fixed 2 Nov); recompute
            bg = background(name, hist, iH, tau_end=tau_end, Nb0=Nb0) if bg is None else bg
            R['Rb'] = M.rbar_rows(bg, R['rows'])
            R['m2xbar'] = M.mean_extra(bg, R['rows'], lamIH, Aterm)
            np.savez_compressed(path, res=json.dumps(R['res'], default=float), rows=R['rows'], spec=R['spec'], ks=R['ks'],
                                wq=R['wq'], Rb=R['Rb'], m2xbar=R['m2xbar'], rbv=2)
        return R
    bg = background(name, hist, iH, tau_end=tau_end, Nb0=Nb0) if bg is None else bg
    H = 1.29e11 if name == 'M1r' else C.load_bench(name, iH)['H_GeV']
    sm = C.SMTab(sm_mt)
    res, rows, spec, ks, wq = M.run(bg, xi, sm, H, n=n, Kmax=Kmax, lamIH=lamIH, Aterm=Aterm, nsave=900)
    Rb = M.rbar_rows(bg, rows)
    m2xbar = M.mean_extra(bg, rows, lamIH, Aterm)
    res['windows'] = bg.windows(nmax=6, xi=xi) if bg.kind == 'bench' else None
    np.savez_compressed(path, res=json.dumps(res, default=float), rows=rows, spec=spec, ks=ks, wq=wq, Rb=Rb,
                        m2xbar=m2xbar, rbv=2)
    return dict(res=res, rows=rows, spec=spec, ks=ks, wq=wq, Rb=Rb, m2xbar=m2xbar, key=key)


_SM = {}


def smtab(Mt, Mh=C.PDG24['Mh']):
    k = (round(Mt, 4), round(Mh, 3))
    if k not in _SM:
        _SM[k] = C.SMTab(Mt, Mh=Mh)
    return _SM[k]


def evalD(R, xi, H_GeV, Mt, Mh=C.PDG24['Mh'], Ns=(50, 55, 60), tau_min=None, mp2=None, **kw):
    sm = smtab(Mt, Mh)
    tm = R['res']['tau_endinf'] if tau_min is None else tau_min
    return M.danger(R['rows'], R['spec'], R['ks'], R['wq'], R['Rb'], xi, sm, H_GeV, Ns=Ns, tau_min=tm, mp2=mp2,
                    m2xbar=R['m2xbar'], **kw)


def crit_Mt(Ds, Mts):
    """first zero crossing of D(M_t) (increasing M_t) by linear interpolation"""
    for i in range(len(Mts) - 1):
        a, b = Ds[i], Ds[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a < 0 <= b:
            return float(Mts[i] + (0 - a) * (Mts[i + 1] - Mts[i]) / (b - a))
        if (not np.isfinite(a)) and np.isfinite(b) and b >= 0:
            return float(Mts[i + 1])
    return None


def summarise(name, hist, iH, R, xi, H_GeV, Nob):
    Ns = (50, 55, 60, Nob)
    out = dict(name=name, hist=hist, iH=iH, H_GeV=H_GeV, xi_H=xi, N_obs=Nob, run=R['res'], key=R['key'])
    # D on the PDG points with criterion variants
    pdg = {}
    for Mt in MT_PDG:
        d = evalD(R, xi, H_GeV, Mt, Ns=Ns)
        var = {}
        for cH in (0.5, 2.0):
            var[f'cH={cH}'] = evalD(R, xi, H_GeV, Mt, Ns=(60, Nob), cH=cH)['best'][60]['D']
        for f in (0.5, 2.0):
            var[f'kc_f={f}'] = evalD(R, xi, H_GeV, Mt, Ns=(60, Nob), kc_f=f)['best'][60]['D']
        pdg[f'{Mt}'] = dict(best=d['best'], variants_D60=var, deadline_curve=d['deadline_curve'],
                            deadline_60=M.deadline(d['deadline_curve'], 60), deadline_Nobs=M.deadline(d['deadline_curve'], Nob),
                            h_max_canon_over_H=d['h_max_canon_over_H'],
                            sqrt_negV2_at_hmax_over_H=d['sqrt_negV2_at_hmax_over_H'])
    out['PDG'] = pdg
    # D(M_t) grid and critical M_t (N = 60 and N_obs), deadline above it
    grid = {}
    for Mt in MT_GRID:
        d = evalD(R, xi, H_GeV, Mt, Ns=(60, Nob))
        band = [evalD(R, xi, H_GeV, Mt, Ns=(60,), cH=c)['best'][60]['D'] for c in (0.5, 2.0)] + \
               [evalD(R, xi, H_GeV, Mt, Ns=(60,), kc_f=f)['best'][60]['D'] for f in (0.5, 2.0)]
        grid[f'{Mt}'] = dict(D60=d['best'][60]['D'], DNobs=d['best'][Nob]['D'], tau_c60=d['best'][60].get('tau_c'),
                             D60_band=[min(band + [d['best'][60]['D']]), max(band + [d['best'][60]['D']])],
                             deadline_60=M.deadline(d['deadline_curve'], 60),
                             deadline_Nobs=M.deadline(d['deadline_curve'], Nob))
    Mts = list(MT_GRID)
    out['Mt_grid'] = grid
    out['Mt_crit_60'] = crit_Mt([grid[f'{m}']['D60'] for m in Mts], Mts)
    out['Mt_crit_Nobs'] = crit_Mt([grid[f'{m}']['DNobs'] for m in Mts], Mts)
    out['Mt_crit_60_band'] = [crit_Mt([grid[f'{m}']['D60_band'][1] for m in Mts], Mts),
                              crit_Mt([grid[f'{m}']['D60_band'][0] for m in Mts], Mts)]
    # M_h variants (coarse)
    mh = {}
    for Mh in MH_VARIANTS:
        Ds = [evalD(R, xi, H_GeV, Mt, Mh=Mh, Ns=(60,))['best'][60]['D'] for Mt in MT_COARSE]
        mh[f'{Mh}'] = dict(D60=Ds, Mt_crit_60=crit_Mt(Ds, list(MT_COARSE)))
    out['Mh_variants'] = dict(Mt=MT_COARSE, rows=mh)
    return out


def cmd_bench():
    out = {}
    path = C.RES / 'bench_D.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    jobs = [('M1r', 'G05', 0), ('HT6t', 'MD', 0), ('HT6t', 'G05', 0), ('T1', 'MD', 0), ('T1', 'MDf', 0), ('T1', 'G05', 0),
            ('E1', 'MDf', 0), ('E1', 'MD', 0), ('E1', 'G05', 0), ('HT6t', 'MD', 1), ('HT6t', 'G05', 1), ('T1', 'MDf', 1),
            ('E1', 'MDf', 1)]
    for name, hist, iH in jobs:
        k = f'{name}|{hist}|iH{iH}'
        if k in out:
            continue
        t0 = time.time()
        H = 1.29e11 if name == 'M1r' else C.load_bench(name, iH)['H_GeV']
        xi = smtab(172.57).xi(H, 0.1)
        R = run_cached(name, hist, xi, iH)
        Nob = N_obs(name, hist, iH)
        o = summarise(name, hist, iH, R, xi, H, Nob)
        if name == 'M1r':
            o['at_H_1.50e11'] = summarise(name, hist, 1, R, smtab(172.57).xi(1.50e11, 0.1), 1.50e11, Nob)
        out[k] = o
        p = o['PDG']
        print(k, 'xi %.4f' % xi, 'Nobs', Nob, 'D60', [round(p[f'{m}']['best'][60]['D'], 1) for m in MT_PDG],
              'DNobs', [round(p[f'{m}']['best'][Nob]['D'], 1) for m in MT_PDG],
              'Mt_crit60 %s Nobs %s band %s' % (o['Mt_crit_60'], o['Mt_crit_Nobs'], o['Mt_crit_60_band']),
              'run %.0fs total %.0fs' % (o['run']['seconds'], time.time() - t0), flush=True)
        C.jdump(out, path)


def cmd_xiscan(name, hist, iH=0):
    path = C.RES / f'xiscan_{name}_{hist}_iH{iH}.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    H = 1.29e11 if name == 'M1r' else C.load_bench(name, iH)['H_GeV']
    Nob = N_obs(name, hist, iH)
    bg = None
    for xe in (0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.14):
        if f'{xe}' in out:
            continue
        xi = smtab(172.57).xi(H, xe)
        if bg is None:
            bg = background(name, hist, iH)
        R = run_cached(name, hist, xi, iH, bg=bg)
        row = dict(xi_EW=xe, xi_H=xi)
        for Mt in MT_PDG + (172.95, 173.05, 173.15):
            d = evalD(R, xi, H, Mt, Ns=(60, Nob))
            row[f'{Mt}'] = dict(D60=d['best'][60]['D'], DNobs=d['best'][Nob]['D'], tau_c=d['best'][60].get('tau_c'),
                                band=[evalD(R, xi, H, Mt, Ns=(60,), cH=c)['best'][60]['D'] for c in (0.5, 2.0)] +
                                     [evalD(R, xi, H, Mt, Ns=(60,), kc_f=f)['best'][60]['D'] for f in (0.5, 2.0)])
        out[f'{xe}'] = row
        print(name, hist, iH, xe, 'xiH %.4f' % xi, {m: round(row[f'{m}']['D60'], 1) for m in MT_PDG}, '%.0fs' % R['res']['seconds'],
              flush=True)
        C.jdump(out, path)
    # edges: D60 = 0 crossing in xi_EW at each M_t (D falls with xi)
    xs = sorted(float(k) for k in out if k[0].isdigit())
    edges = {}
    for Mt in MT_PDG + (172.95, 173.05, 173.15):
        for lab, getter in (('D60', lambda r: r['D60']), ('DNobs', lambda r: r['DNobs']),
                            ('band_hi', lambda r: max(r['band'] + [r['D60']])), ('band_lo', lambda r: min(r['band'] + [r['D60']]))):
            Ds = [getter(out[f'{x}'][f'{Mt}']) for x in xs]
            e = None
            for i in range(len(xs) - 1):
                if Ds[i] >= 0 > Ds[i + 1]:
                    e = xs[i] + (0 - Ds[i]) * (xs[i + 1] - xs[i]) / (Ds[i + 1] - Ds[i])
            if Ds and Ds[0] < 0:
                e = f'< {xs[0]}'
            edges[f'{Mt} {lab}'] = e
    out['edges_xiEW'] = edges
    out['N_obs'] = Nob
    C.jdump(out, path)
    print(json.dumps(edges, indent=0))


def cmd_upper(name, hist, iH=0, xis=(1.0, 3.0, 10.0, 30.0, 100.0, 200.0, 300.0, 1000.0)):
    """large xi(H_i): mean field (2 sigma_max vs h_max) and D; start 2 e-folds before the end (heavy Higgs)."""
    path = C.RES / f'upper_{name}_{hist}_iH{iH}.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    H = 1.29e11 if name == 'M1r' else C.load_bench(name, iH)['H_GeV']
    Nob = N_obs(name, hist, iH)
    bgs = {}
    for xi in xis:
        if f'{xi}' in out:
            continue
        Nb0 = 3.0 if xi <= 30 else 2.0
        te = 150.0
        if Nb0 not in bgs:
            bgs[Nb0] = background(name, hist, iH, tau_end=te, Nb0=Nb0)
        n = 320 if xi <= 100 else 400
        R = run_cached(name, hist, xi, iH, Nb0=Nb0, n=n, tau_end=te, bg=bgs[Nb0])
        rows = R['rows']
        te_inf = R['res']['tau_endinf']
        post = rows[:, 0] > te_inf
        s2 = rows[:, 5] + rows[:, 6]
        sp = M.sigma2_produced(R['res'], rows, R['Rb'], R['ks'], R['wq'], xi)
        late = rows[:, 0] > te_inf + 0.5          # after the first oscillations (production epoch)
        row = dict(xi_H=xi, xi_EW=smtab(172.57).xi_EW(H, xi), runaway=R['res']['runaway'],
                   max_sigma2_post_instantaneous=float(np.max(s2[post])),
                   max_sigma2_produced=float(np.max(sp[late])) if late.any() else None,
                   tau_at_max_produced=float(rows[late][np.argmax(sp[late]), 0]) if late.any() else None,
                   max_nk=float(np.max(R['res']['nk_end'])), k_at_max_nk=float(R['ks'][int(np.argmax(R['res']['nk_end']))]),
                   tau_end=R['res']['tau_end'], seconds=R['res']['seconds'])
        for Mt in MT_PDG:
            sm = smtab(Mt)
            hmc = math.exp(sm.G(H)) * sm.flat['h_max'] / H
            d = evalD(R, xi, H, Mt, Ns=(60, Nob))
            row[f'{Mt}'] = dict(two_sigma_prod_over_hmax=(2 * math.sqrt(row['max_sigma2_produced']) / hmc
                                                          if row['max_sigma2_produced'] is not None else None),
                                D60=d['best'][60]['D'], tau_c=d['best'][60].get('tau_c'))
        out[f'{xi}'] = row
        print(name, hist, xi, {k: (v if not isinstance(v, dict) else {kk: round(vv, 2) if isinstance(vv, float) else vv
                                                                      for kk, vv in v.items()}) for k, v in row.items()},
              flush=True)
        C.jdump(out, path)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'bench':
        cmd_bench()
    elif cmd == 'xiscan':
        cmd_xiscan(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    elif cmd == 'upper':
        xs = tuple(float(x) for x in sys.argv[5].split(',')) if len(sys.argv) > 5 else None
        if xs:
            cmd_upper(sys.argv[2], sys.argv[3], int(sys.argv[4]), xis=xs)
        else:
            cmd_upper(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 0)
