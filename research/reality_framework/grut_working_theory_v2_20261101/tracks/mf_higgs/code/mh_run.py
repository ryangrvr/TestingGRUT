#!/usr/bin/env python3
"""Track mf_higgs driver: Higgs survival after inflation on HFa / HFb (own code; one process, foreground).

  python3 -B mh_run.py bench [names] [hists]   # xi_EW = 0.1: D(M_t) grav-only and with the N -> L H plasma; M_t,crit;
                                               # criterion band; M_h, alpha_s variants; deadlines -> results/post_D.json
  python3 -B mh_run.py plasma                  # plasma timing table (fluid) -> results/plasma.json
  python3 -B mh_run.py xiscan <name> <hist>    # lower xi edge after inflation -> results/xiscan_<name>_<hist>.json
  python3 -B mh_run.py upper <name> <hist>     # large xi (upper edge) -> results/upper_<name>_<hist>.json
  python3 -B mh_run.py Hi <name> <hist>        # H_i dependence (fixed dimensionless history) -> results/Hi_<name>_<hist>.json
  python3 -B mh_run.py checks                  # numerical controls -> results/controls.json
Mode runs are cached in the scratchpad (npz).  Every history and criterion is a model choice (mh_common, mh_modes).
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np

import mh_common as C
import mh_modes as M

MT_PDG = (172.28, 172.57, 172.86)
MT_EXT = (172.28, 172.57, 172.86, 173.15, 173.45, 173.75, 174.05, 174.35, 174.6)
MT_FINE = tuple(np.round(np.arange(172.60, 173.31, 0.05), 3))
TAU_THS = (0.25, 0.3, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17, 20, 25, 30, 40, 60, 80)
_SMT = {}


def smt(Mt, Mh=C.PDG24['Mh'], als=C.PDG24['als']):
    k = (round(Mt, 4), round(Mh, 3), round(als, 5))
    if k not in _SMT:
        _SMT[k] = C.IHC.SMTab(round(Mt, 4), Mh=Mh, als=als)
    return _SMT[k]


_BG = {}


def bgc(name, hist, Hi=None, osc=True):
    k = (name, hist, Hi, osc)
    if k not in _BG:
        gs = 1.0 if Hi is None else Hi / C.H_REF
        _BG[k] = C.Background(name, hist, osc=osc, grest_scale=gs, Hi_GeV=Hi)
    return _BG[k]


def run_cached(name, hist, xi, Hi=None, osc=True, t0=-3.0, t_end=80.0, n=260, Kmax=None, sm_mt=172.57, tag='',
               m2x=None, cK=0.05, rerun=False):
    key = f'{name}_{hist}_Hi{Hi}_osc{int(osc)}_xi{xi:.6g}_t0{t0}_te{t_end}_n{n}_K{Kmax}_mt{sm_mt}_cK{cK}{tag}'
    p = C.SP / 'runs' / (key + '.npz')
    p.parent.mkdir(exist_ok=True)
    if p.exists() and not rerun:
        z = np.load(p, allow_pickle=True)
        return dict(res=json.loads(str(z['res'])), rows=z['rows'], spec=z['spec'], ks=z['ks'], wq=z['wq'], key=key)
    bg = bgc(name, hist, Hi, osc)
    Hi_ = bg.Hi_GeV
    R = M.run(bg, xi, smt(sm_mt) if Hi is None else smt(sm_mt), t0=t0, t_end=t_end, n=n, Kmax=Kmax, m2x=m2x, cK=cK)
    R['res']['Hi_GeV'] = Hi_
    np.savez_compressed(p, res=json.dumps(R['res'], default=float), rows=R['rows'], spec=R['spec'], ks=R['ks'],
                        wq=R['wq'])
    R['key'] = key
    return R


def xi_at(Hi, xi_EW=0.1, Mt=172.57):
    return smt(Mt).xi(Hi, xi_EW)


def D_of(R, Mt, N, plasma=None, Mh=C.PDG24['Mh'], als=C.PDG24['als'], **kw):
    """plasma: None (gravitational only), or a factor f multiplying the stored hard-gas m_ind^2 column."""
    s = smt(Mt, Mh, als)
    mp2 = None if plasma is None else plasma * R['rows'][:, 8]
    return M.danger(R, s, R['res']['xi'], N, mp2=mp2, **kw)


def deadline_curve(R, Mt, N):
    out = {}
    for tt in TAU_THS:
        out[str(tt)] = D_of(R, Mt, N, tau_th=tt)['D']
    return out


def deadline(curve):
    items = sorted((float(k), v) for k, v in curve.items())
    prev = None
    for t, D in items:
        if np.isfinite(D) and D >= 0:
            if prev is None or not np.isfinite(prev[1]):
                return t
            return prev[0] + (0 - prev[1]) * (t - prev[0]) / (D - prev[1])
        prev = (t, D)
    return None


def fmt(x):
    return round(x, 2) if isinstance(x, float) and np.isfinite(x) else x


def cmd_bench(names=('HFa', 'HFb'), hists=('EXP', 'UNBL', 'PAULI0', 'PAULI10', 'PAULIb')):
    path = C.RES / 'post_D.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    for name in names:
        for hist in hists:
            k = f'{name}|{hist}'
            if k in out:
                continue
            t0 = time.time()
            bg = bgc(name, hist)
            xi = xi_at(bg.Hi_GeV)
            R = run_cached(name, hist, xi)
            Nob = R['res']['N_obs']
            o = dict(name=name, hist=hist, xi_H=xi, N_obs=Nob, K=R['res']['K'], run={kk: v for kk, v in R['res'].items()
                                                                                    if kk != 'nk_end'})
            # gravitational only and with the plasma, on the extended grid
            ext = {}
            for Mt in MT_EXT:
                g = D_of(R, Mt, Nob)
                g60 = D_of(R, Mt, 60.0)['D']
                pl = D_of(R, Mt, Nob, plasma=1.0)
                pl2 = D_of(R, Mt, Nob, plasma=1e-2)['D']
                pl4 = D_of(R, Mt, Nob, plasma=1e-4)['D']
                ext[f'{Mt}'] = {'grav': g, 'grav_D60': g60, 'plasma': pl, 'plasma_c1e-2': pl2, 'plasma_c1e-4': pl4}
            o['Mt_ext'] = ext
            fine = {}
            for Mt in MT_FINE:
                g = D_of(R, Mt, Nob)['D']
                band = [D_of(R, Mt, Nob, cH=c)['D'] for c in (0.5, 2.0)] + [D_of(R, Mt, Nob, kc_f=f)['D'] for f in (0.5, 2.0)]
                fine[f'{Mt}'] = dict(D=g, band=[min(band + [g]), max(band + [g])])
            o['Mt_fine'] = fine
            Ms = list(MT_FINE)
            o['Mt_crit_grav'] = M.crit(Ms, [fine[f'{m}']['D'] for m in Ms])
            o['Mt_crit_grav_band'] = [M.crit(Ms, [fine[f'{m}']['band'][1] for m in Ms]),
                                      M.crit(Ms, [fine[f'{m}']['band'][0] for m in Ms])]
            o['Mt_crit_grav_N60'] = M.crit(list(MT_EXT), [ext[f'{m}']['grav_D60'] for m in MT_EXT])
            Me = list(MT_EXT)
            o['Mt_crit_plasma'] = M.crit(Me, [ext[f'{m}']['plasma']['D'] for m in Me])
            o['Mt_crit_plasma_c1e-2'] = M.crit(Me, [ext[f'{m}']['plasma_c1e-2'] for m in Me])
            o['Mt_crit_plasma_c1e-4'] = M.crit(Me, [ext[f'{m}']['plasma_c1e-4'] for m in Me])
            # M_h and alpha_s variants (gravitational only)
            var = {}
            for lab, Mh, als in (('Mh-1s', 125.09, 0.1180), ('Mh+1s', 125.31, 0.1180), ('als-1s', 125.20, 0.1171),
                                 ('als+1s', 125.20, 0.1189)):
                Mts = (172.28, 172.57, 172.86, 173.15, 173.45) if 'als' in lab else (172.57, 172.86, 173.15)
                Ds = [D_of(R, Mt, Nob, Mh=Mh, als=als)['D'] for Mt in Mts]
                var[lab] = dict(Mt=Mts, D=Ds, Mt_crit=M.crit(list(Mts), Ds))
            o['variants'] = var
            # deadlines (sharp plasma onset) at the PDG points and above
            dl = {}
            for Mt in (172.86, 173.0, 173.15, 173.45, 173.75):
                cv = deadline_curve(R, Mt, Nob)
                dl[f'{Mt}'] = dict(curve=cv, deadline=deadline(cv))
            o['deadlines'] = dl
            out[k] = o
            C.jdump(out, path)
            print(k, 'xi %.4f N_obs %.2f K %.2f' % (xi, Nob, o['K']),
                  'D grav', [fmt(ext[f'{m}']['grav']['D']) for m in MT_PDG],
                  'D plasma', [fmt(ext[f'{m}']['plasma']['D']) for m in MT_EXT],
                  'crit grav %s band %s | plasma %s (c1e-2 %s, c1e-4 %s)' % (o['Mt_crit_grav'], o['Mt_crit_grav_band'],
                                                                            o['Mt_crit_plasma'], o['Mt_crit_plasma_c1e-2'],
                                                                            o['Mt_crit_plasma_c1e-4']),
                  'als', {kk: v['Mt_crit'] for kk, v in var.items()},
                  'deadlines', {kk: fmt(v['deadline']) for kk, v in dl.items()},
                  'run %.0fs total %.0fs' % (R['res']['seconds'], time.time() - t0), flush=True)


def thermal_time(f, mPr, alpha=1 / 40.0):
    """inf_reheat's pessimistic (LPM-type) thermalisation: first t with t - t_start >= alpha^-2 (E/T_*)^{1/2}/T_*,
    T_* = (30 rho_SM/(pi^2 g*))^{1/4} [E]."""
    t = f['t']
    rho = f['f_SM'] * 3 * f['H'] ** 2 * mPr ** 2
    Ts = (30 * np.maximum(rho, 0) / (math.pi ** 2 * 106.75)) ** 0.25
    E = np.where(np.isfinite(f['E_SM']), f['E_SM'], np.inf)
    tth = np.where(Ts > 0, alpha ** -2 * np.sqrt(np.maximum(E / np.maximum(Ts, 1e-300), 1.0)) / np.maximum(Ts, 1e-300), np.inf)
    i = np.where((t - t[0]) >= tth)[0]
    return dict(t=float(t[i[0]]) if len(i) else None, T_GeV=float(Ts[i[0]] * C.H_REF) if len(i) else None,
                kappaT2_over_Hi2=float(0.11 * Ts[i[0]] ** 2) if len(i) else None)


def cmd_plasma():
    """plasma timing from the actual decay chain (fluid), for every history and the m-tilde variant."""
    path = C.RES / 'plasma.json'
    out = {}
    cases = [('UNBL', dict(blocking=False), 1.0), ('PAULI0', dict(blocking=True, Gth=0.0), 1.0),
             ('PAULI1', dict(blocking=True, Gth=1.0), 1.0), ('PAULI10', dict(blocking=True, Gth=10.0), 1.0),
             ('PAULIb', dict(blocking=True, Gth=0.0, pert_free=False), 1.0),
             ('UNBL_mt0.015', dict(blocking=False), 0.3), ('PAULI0_mt0.015', dict(blocking=True, Gth=0.0), 0.3)]
    mPr = C.MPL / C.H_REF
    for name in ('HFa', 'HFb'):
        o = {}
        for lab, kw, gs in cases:
            if name == 'HFb' and lab not in ('UNBL', 'PAULI0'):
                continue
            f = C.fluid_cached(name, grest_scale=gs, **kw)
            t = f['t']
            m2 = np.where(np.isfinite(f['E_SM']) & (f['E_SM'] > 0),
                          C.C_HARD * f['f_SM'] * 3 * f['H'] ** 2 * mPr ** 2 / np.where(f['E_SM'] > 0, f['E_SM'], 1) ** 2, 0)
            g = lambda key, x: float(np.interp(x, t, f[key]))
            o[lab] = dict(K=f['meta']['K_final'], start=f['meta']['start'],
                          f_rad={f'{x}': g('f_rad', x) for x in (0.3, 0.6, 1, 3, 10, 100)},
                          t_fSM={f'{thr:g}': C.first_time(t, f['f_SM'], thr) for thr in (1e-9, 1e-6, 1e-3, 0.1, 0.5)},
                          t_m2ind={f'{thr:g}': C.first_time(t, m2, thr) for thr in (1e-3, 1.0, 10.0, 100.0, 1e3)},
                          m2ind_at={f'{x}': float(np.interp(x, t, m2)) for x in (0.195, 0.2, 0.25, 0.3, 0.5, 1, 3, 10, 30, 80)},
                          m2ind_over_H2_at={f'{x}': float(np.interp(x, t, m2) / np.interp(x, t, f['H']) ** 2)
                                            for x in (0.2, 0.3, 1, 10, 80)},
                          E_SM_at={f'{x}': g('E_SM', x) for x in (0.2, 0.3, 1, 10)},
                          E_N_at={f'{x}': g('E_N', x) for x in (0.2, 0.3, 1)},
                          T_rad_equiv_GeV_at_10=float((30 * g('f_SM', 10) * 3 * g('H', 10) ** 2 * mPr ** 2 /
                                                       (math.pi ** 2 * 106.75)) ** 0.25 * C.H_REF),
                          t_thermal_pessimistic=thermal_time(f, mPr))
            print(name, lab, json.dumps({kk: o[lab][kk] for kk in ('K', 't_fSM', 't_m2ind')}, default=fmt), flush=True)
        out[name] = o
        C.jdump(out, path)


def cmd_xiscan(name, hist):
    path = C.RES / f'xiscan_{name}_{hist}.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    bg = bgc(name, hist)
    for xe in (0.02, 0.04, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.14):
        if f'{xe}' in out:
            continue
        xi = xi_at(bg.Hi_GeV, xe)
        R = run_cached(name, hist, xi)
        Nob = R['res']['N_obs']
        row = dict(xi_EW=xe, xi_H=xi)
        for Mt in (172.28, 172.57, 172.86, 173.15):
            row[f'{Mt}'] = dict(grav=D_of(R, Mt, Nob)['D'], plasma=D_of(R, Mt, Nob, plasma=1.0)['D'],
                                band=[D_of(R, Mt, Nob, cH=c)['D'] for c in (0.5, 2.0)] +
                                     [D_of(R, Mt, Nob, kc_f=f)['D'] for f in (0.5, 2.0)])
        out[f'{xe}'] = row
        print(name, hist, xe, {m: (fmt(row[f'{m}']['grav']), fmt(row[f'{m}']['plasma'])) for m in (172.28, 172.57, 172.86, 173.15)},
              flush=True)
        C.jdump(out, path)
    xs = sorted(float(k) for k in out if k[0].isdigit())
    edges = {}
    for Mt in (172.28, 172.57, 172.86, 173.15):
        for lab, get in (('grav', lambda r: r['grav']), ('band_hi', lambda r: max(r['band'] + [r['grav']])),
                         ('band_lo', lambda r: min(r['band'] + [r['grav']])), ('plasma', lambda r: r['plasma'])):
            Ds = [get(out[f'{x}'][f'{Mt}']) for x in xs]
            e = None
            for i in range(len(xs) - 1):
                if np.isfinite(Ds[i]) and Ds[i] >= 0 > Ds[i + 1]:
                    e = xs[i] + (0 - Ds[i]) * (xs[i + 1] - xs[i]) / (Ds[i + 1] - Ds[i])
                elif Ds[i] >= 0 and not np.isfinite(Ds[i + 1]):
                    e = f'between {xs[i]} and {xs[i + 1]}'
            if Ds and (Ds[0] < 0 or not np.isfinite(Ds[0])):
                e = f'< {xs[0]}'
            edges[f'{Mt} {lab}'] = e
    out['edges_xiEW'] = edges
    C.jdump(out, path)
    print(json.dumps(edges, indent=0))


def cmd_upper(name, hist, xis=(1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0)):
    path = C.RES / f'upper_{name}_{hist}.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    bg = bgc(name, hist)
    for xi in xis:
        if f'{xi}' in out:
            continue
        n = 320 if xi <= 100 else 400
        R = run_cached(name, hist, xi, t0=-1.0, t_end=120.0, n=n)
        rows = R['rows']
        te = R['res']['tau_endinf']
        post = rows[:, 0] > te
        nk = np.maximum(np.array(R['res']['nk_end']), 0.0)
        # produced variance, mapped back adiabatically (inf_higgs convention)
        sp = []
        for i in np.where(post)[0]:
            a, H, Rb = rows[i, 1], rows[i, 2], rows[i, 4]
            w = np.sqrt(R['ks'] ** 2 / a ** 2 + max(xi * Rb, 0.0) + H ** 2)
            sp.append(float(np.sum(R['wq'] * nk / (a ** 3 * w))))
        sp = np.array(sp)
        late = rows[post, 0] > te + 0.5
        row = dict(xi_H=xi, xi_EW=smt(172.57).xi_EW(bg.Hi_GeV, xi), runaway=R['res']['runaway'],
                   max_sigma2_inst_post=float(np.max(rows[post, 6])),
                   max_sigma2_produced=float(np.max(sp[late])) if late.any() else None,
                   max_nk=float(np.max(nk)), k_at_max_nk=float(R['ks'][int(np.argmax(nk))]),
                   tau_end=R['res']['t_end'], seconds=R['res']['seconds'])
        for Mt in MT_PDG:
            s = smt(Mt)
            hmc = math.exp(s.G(bg.Hi_GeV)) * s.flat['h_max'] / bg.Hi_GeV
            row[f'{Mt}'] = dict(two_sigma_prod_over_hmax=(2 * math.sqrt(row['max_sigma2_produced']) / hmc
                                                          if row['max_sigma2_produced'] else None),
                                grav=D_of(R, Mt, R['res']['N_obs'])['D'], plasma=D_of(R, Mt, R['res']['N_obs'], plasma=1.0)['D'])
        out[f'{xi}'] = row
        print(name, hist, xi, json.dumps(row, default=fmt), flush=True)
        C.jdump(out, path)


def cmd_Hi(name, hist, His=(6e10, 1.0e11, 2.1e11, 3e11, 5e11)):
    path = C.RES / f'Hi_{name}_{hist}.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    for Hi in His:
        if f'{Hi:.3g}' in out:
            continue
        bg = bgc(name, hist, Hi)
        xi = xi_at(Hi)
        R = run_cached(name, hist, xi, Hi=Hi)
        Nob = R['res']['N_obs']
        Ds = [D_of(R, Mt, Nob)['D'] for Mt in MT_FINE]
        Dp = [D_of(R, Mt, Nob, plasma=1.0)['D'] for Mt in MT_EXT]
        o = dict(H_i=Hi, xi_H=xi, N_obs=Nob, Mt_crit_grav=M.crit(list(MT_FINE), Ds), D_grav_fine=Ds,
                 D_plasma_ext=Dp, Mt_crit_plasma=M.crit(list(MT_EXT), Dp),
                 tau_m2ind_1=C.first_time(bg.pl_t, bg.pl_m2, 1.0))
        if o['Mt_crit_grav'] is None and Ds[0] >= 0:
            Dl = [D_of(R, Mt, Nob)['D'] for Mt in (171.5, 171.8, 172.1, 172.4, 172.6)]
            o['Mt_crit_grav'] = M.crit([171.5, 171.8, 172.1, 172.4, 172.6], Dl)
            o['D_grav_low'] = Dl
        out[f'{Hi:.3g}'] = o
        print(name, hist, 'H_i %.3g' % Hi, 'crit grav', o['Mt_crit_grav'], 'crit plasma', o['Mt_crit_plasma'],
              'tau(m2>=1)', o['tau_m2ind_1'], flush=True)
        C.jdump(out, path)


if __name__ == '__main__':
    c = sys.argv[1]
    if c == 'bench':
        names = tuple(sys.argv[2].split(',')) if len(sys.argv) > 2 else ('HFa', 'HFb')
        hists = tuple(sys.argv[3].split(',')) if len(sys.argv) > 3 else ('EXP', 'UNBL', 'PAULI0', 'PAULI10', 'PAULIb')
        cmd_bench(names, hists)
    elif c == 'plasma':
        cmd_plasma()
    elif c == 'xiscan':
        cmd_xiscan(sys.argv[2], sys.argv[3])
    elif c == 'upper':
        cmd_upper(sys.argv[2], sys.argv[3])
    elif c == 'Hi':
        cmd_Hi(sys.argv[2], sys.argv[3])
