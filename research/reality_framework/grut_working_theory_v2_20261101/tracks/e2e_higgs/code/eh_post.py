#!/usr/bin/env python3
"""e2e_higgs: Higgs modes and the danger exponent D after inflation on the final histories (mf_higgs code, read-only;
only the background and the medium column change).  [C] on the argued criterion [E]; every choice a model choice.

  python3 -B eh_post.py main     -> results/post.json   (D grid, M_t,crit + criterion band, N_obs variants, deadlines,
                                                         alpha_s and M_h variants, R-fluctuation variants, media)
  python3 -B eh_post.py xiscan   -> results/xiscan.json (lower xi edge after inflation, gravitational only)
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np

import eh_common as E

MM, MC = E.MM, E.MC
HISTS = ('KIN', 'LAT_P0', 'LAT_GM', 'LATF_P0')
HIS = (5e10, 1.3e11, 2.4e11)
MT_COARSE = (172.28, 172.57, 172.86, 173.15, 173.45, 173.75, 174.05, 174.35, 174.6)
MT_FINE = tuple(np.round(np.arange(172.40, 173.351, 0.05), 3))
TAU_THS = (0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 60, 80)
K_LAT_OWN = {'LAT_P0': 3.59, 'LAT_GM': 5.65, 'LATF_P0': 3.59}    # e2e_lattice HB.stage2 K (its digest: 3.6-4 / 5.7-6)


def N_of(K, Hi):
    return MC.N_obs('HFa', K) + 0.5 * math.log(Hi / E.H_BENCH)


def D(R, Mt, N, Mh=E.PDG['Mh'], als=E.PDG['als'], **kw):
    return MM.danger(R, E.smt(Mt, Mh, als), R['res']['xi'], N, **kw)


def crit_scan(R, N, Ms=MT_FINE, Mh=E.PDG['Mh'], als=E.PDG['als'], **kw):
    Ds = [D(R, m, N, Mh, als, **kw)['D'] for m in Ms]
    return MM.crit(list(Ms), Ds), Ds


def with_R(R, Rb_new):
    R2 = dict(R)
    rows = R['rows'].copy()
    rows[:, 4] = Rb_new
    R2['rows'] = rows
    return R2


def deadline(curve):
    prev = None
    for t, Dv in curve:
        if np.isfinite(Dv) and Dv >= 0:
            if prev is None or not np.isfinite(prev[1]):
                return t
            return prev[0] + (0 - prev[1]) * (t - prev[0]) / (Dv - prev[1])
        prev = (t, Dv)
    return None


def fm(x):
    return None if x is None else (round(float(x), 3) if np.isfinite(x) else str(x))


def cmd_main():
    path = E.RES / 'post.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    tw = time.time()
    for hist in HISTS:
        for Hi in HIS:
            key = f'{hist}@{Hi:.3g}'
            if key in out:
                continue
            R = E.run_cached(hist, Hi)
            bg = E.FinalBG(hist, Hi)
            res = R['res']
            Kkin = E.K_KIN[Hi]
            # N_obs choices: the history's own K (spliced/own), e2e_kinetic's K(H_i), the lattice's own fluid K
            Ns = {'K_hist': (bg.K_hist if Hi == 1.3e11 or hist != 'KIN' else Kkin), 'K_kinetic': Kkin}
            if hist in K_LAT_OWN:
                Ns['K_lattice_fluid'] = K_LAT_OWN[hist]
            Nv = {k: N_of(K, Hi) for k, K in Ns.items()}
            Nc = Nv['K_hist'] if hist != 'KIN' else N_of(Kkin, Hi)
            o = dict(hist=hist, H_i=Hi, xi_H_at_Hi=res['xi'], K_choices=Ns, N_obs=Nv, N_central=Nc,
                     K_hist_dimless_bg=bg.K_hist, wronskian_err_end=res['wronskian_err_end'], lna_mismatch=bg.lna_mismatch,
                     max_sigma2=res['max_sigma2'])
            o['D_grav_coarse'] = {f'{m}': fm(D(R, m, Nc)['D']) for m in MT_COARSE}
            mc, Ds = crit_scan(R, Nc)
            o['Mt_crit_grav'] = mc
            o['D_grav_fine'] = {f'{m}': fm(d) for m, d in zip(MT_FINE, Ds)}
            # criterion band (as mf_higgs): cH = 0.5 / 2, k_c factor 0.5 / 2
            bands = []
            for kw in (dict(cH=0.5), dict(cH=2.0), dict(kc_f=0.5), dict(kc_f=2.0)):
                bands.append(crit_scan(R, Nc, **kw)[0])
            bb = [b for b in bands + [mc] if b is not None]
            o['Mt_crit_band'] = [min(bb), max(bb)] if bb else None
            o['Mt_crit_by_N'] = {k: crit_scan(R, N)[0] for k, N in Nv.items()}
            o['Mt_crit_N60'] = crit_scan(R, 60.0)[0]
            info = D(R, 172.86, Nc)
            o['roll_172.86'] = {k: fm(v) if isinstance(v, (float, np.floating)) else v for k, v in info.items()}
            # deadlines: D(tau_th) = 0 at several M_t (the latest time a Higgs-coupled medium must appear)
            o['deadline'] = {}
            for m in (173.0, 173.15, 173.45, 173.75, 174.05, 174.6):
                cur = [(tt, D(R, m, Nc, tau_th=tt)['D']) for tt in TAU_THS]
                o['deadline'][f'{m}'] = fm(deadline(cur))
            # R fluctuations (lattice R_rms): pessimistic local variant  Rbar -> Rbar - n R_rms  [A]
            t = R['rows'][:, 0]
            Rb = R['rows'][:, 4]
            rms = bg.Rrms(t)
            o['R_rms_variants'] = {}
            for nsig in (1.0, 2.0):
                R2 = with_R(R, Rb - nsig * rms)
                o['R_rms_variants'][f'minus_{nsig:g}sigma'] = dict(Mt_crit=crit_scan(R2, Nc)[0],
                                                                    D_172p86=fm(D(R2, 172.86, Nc)['D']))
            # media (N gas with e2e_sm Yukawas + SM: e2e_kinetic's kH Tcal^2 on KIN, hard gas on the lattice f_SM)
            mN = bg.col('m2_N', t)
            mS = bg.col('m2_SM', t)
            mNy = bg.col('m2_N_ymin', t)
            med = {}
            for lab, m2 in (('N+SM', mN + mS), ('N_gas', mN), ('N_gas_ymin', mNy), ('SM', mS)):
                Dmax = max(D(R, m, Nc, mp2=m2)['D'] for m in (174.05, 174.6))
                # smallest factor f (decades) on m2 with D(174.6) < 0
                fsafe = None
                for f in 10.0 ** np.arange(0, -9, -0.5):
                    if D(R, 174.6, Nc, mp2=f * m2)['D'] < 0:
                        fsafe = f
                    else:
                        break
                med[lab] = dict(D_max_174_to_174p6=fm(Dmax), smallest_safe_factor_at_174p6=fsafe)
            o['media'] = med
            o['medium_m2_at'] = {f'{x}': dict(N=float(bg.col('m2_N', x)), SM=float(bg.col('m2_SM', x)))
                                 for x in (0.195, 0.2, 0.25, 0.5, 1.0, 3.0, 10.0, 30.0)}
            # M_h and alpha_s (central history entries; finer grids on the design H_i only)
            if Hi == 1.3e11 or hist == 'KIN':
                o['Mh'] = {}
                for Mh in (125.09, 125.31):
                    o['Mh'][f'{Mh}'] = crit_scan(R, Nc, Ms=tuple(np.round(np.arange(172.50, 173.351, 0.05), 3)), Mh=Mh)[0]
                o['alpha_s'] = {}
                for als in (0.1171, 0.1189):
                    Ms = tuple(np.round(np.arange(172.20, 173.651, 0.05), 3))
                    c_, Ds_ = crit_scan(R, Nc, Ms=Ms, als=als)
                    o['alpha_s'][f'{als}'] = dict(Mt_crit=c_, D_at_172p57=fm(Ds_[Ms.index(172.6)] if 172.6 in Ms else None),
                                                  D_central_Mt=fm(D(R, 172.57, Nc, als=als)['D']))
            out[key] = o
            E.jdump(out, path)
            print(key, 'xi %.4f N %.3f crit %s band %s byN %s | Rrms-1/-2 %s %s | media %s | %.0fs' % (
                res['xi'], Nc, fm(mc), o['Mt_crit_band'], {k: fm(v) for k, v in o['Mt_crit_by_N'].items()},
                fm(o['R_rms_variants']['minus_1sigma']['Mt_crit']), fm(o['R_rms_variants']['minus_2sigma']['Mt_crit']),
                {k: (v['D_max_174_to_174p6'], v['smallest_safe_factor_at_174p6']) for k, v in med.items()},
                time.time() - tw), flush=True)
            if 'alpha_s' in o:
                print('    Mh', o['Mh'], 'als', o['alpha_s'], 'deadline', o['deadline'], flush=True)


def cmd_xiscan():
    path = E.RES / 'xiscan.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    xis = (0.05, 0.065, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3)
    for hist, Hi in [(h, 1.3e11) for h in HISTS] + [('KIN', 5e10), ('KIN', 2.4e11)]:
        key = f'{hist}@{Hi:.3g}'
        if key in out:
            continue
        bgK = E.FinalBG(hist, Hi).K_hist
        Nc = N_of(E.K_KIN[Hi] if hist == 'KIN' else bgK, Hi)
        rows = {}
        for xe in xis:
            R = E.run_cached(hist, Hi, xi_EW=xe)
            rows[f'{xe}'] = {f'{m}': fm(D(R, m, Nc)['D']) for m in (172.57, 172.86, 173.15)}
        edge = {}
        for m in ('172.57', '172.86', '173.15'):
            ys = [rows[f'{x}'][m] for x in xis]
            e = None
            for i in range(len(xis) - 1):
                a, b = ys[i], ys[i + 1]
                if isinstance(a, float) and isinstance(b, float) and a >= 0 > b:
                    e = xis[i] + (0 - a) * (xis[i + 1] - xis[i]) / (b - a)
                    break
            if e is None and isinstance(ys[0], float) and ys[0] < 0:
                e = f'< {xis[0]}'
            edge[m] = e
        out[key] = dict(N_central=Nc, D=rows, xi_EW_lower_edge_grav=edge)
        E.jdump(out, path)
        print(key, edge, flush=True)


if __name__ == '__main__':
    {'main': cmd_main, 'xiscan': cmd_xiscan}[sys.argv[1]]()
