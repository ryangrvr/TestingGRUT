#!/usr/bin/env python3
"""Skeptic of inf_higgs: independent D at xi_EW = 0.1 on the benchmarks (own code sk_common.py).

  python3 -B sk_bench.py run <name> <hist> [iH]     # one history -> scratch npz + results/sk_bench.json entry
     hist: MD (Gamma = 0), G05 (Gamma = 0.5), MDf (w = 0 fluid at the track's fragmentation time),
           RDf (radiation-like w = 1/3 fluid at the fragmentation time: a 'void' bracket for the fragmented T1/E1),
           MDnoH (MD with Hartree off)
  python3 -B sk_bench.py Mtgrid                      # D(M_t), M_t,crit for every stored run (N = 60 and the track's N_obs)
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np

import sk_common as C

MT_GRID = sorted(set([172.28, 172.57, 172.86] + [round(float(x), 3) for x in np.arange(172.70, 173.46, 0.05)]))
TFRAG = {'T1': 0.43, 'E1': 0.35}         # inf_potential's linear-backreaction times (track's MDf convention)
H0C_LN = math.log(0.05 / 2.248e-4)


def N_obs_track(name, hist, iH=0, T_RH=1e11):
    """the track's N_obs (e-folds from today's horizon exit to S1) for comparison"""
    r = json.load(open(C.TRACKS / 'inf_potential' / 'results' / f'bench_{name}.json'))['by_Hi'][iH]
    if hist == 'G05':
        Ns = r['Nstar']
    else:
        inst = [x for x in r['Nstar_alternatives'] if x['history'] == ['inst']][0]['Nstar']
        Ns = inst - math.log(r['T_inst_GeV'] / T_RH) / 3
    return Ns + H0C_LN - 1.0


def key(name, hist, iH, H_GeV=None, A=0.0, lam=0.0, tag=''):
    return f'{name}_{hist}_iH{iH}_H{(H_GeV or 0):.4g}_A{A:.4g}_lam{lam:.3g}{tag}'


def run(name, hist, iH=0, H_GeV=None, Aterm=0.0, lamIH=0.0, xi_EW=0.1, tag='', save=True, **kw):
    pot = C.Quad(3.0) if name == 'M1r' else C.Pot(name, iH)
    H = pot.H_GeV if H_GeV is None else H_GeV
    sm0 = C.sm_in(172.57)
    xi = sm0.xi_at(H, xi_EW)
    gam = 0.5 if (hist == 'G05' or name == 'M1r') else 0.0
    frag, fw = None, 0.0
    if hist in ('MDf', 'RDf'):
        frag = TFRAG[name]
        fw = 0.0 if hist == 'MDf' else 1 / 3
    hart = hist != 'MDnoH'
    t0 = time.time()
    ev = C.evolve(pot, xi, sm0, H_GeV=H, gam=gam, tau_end=40.0 if name == 'M1r' else 80.0, hartree=hart,
                  Aterm=Aterm, lamIH=lamIH, frag_after_end=frag, frag_w=fw, **kw)
    ev['seconds'] = time.time() - t0
    ev['Rbar'] = C.period_mean(ev, 1, fixed=(math.pi / 3 if name == 'M1r' else None))
    ev['m2xbar'] = C.period_mean(ev, 5, fixed=(math.pi / 3 if name == 'M1r' else None)) if (Aterm or lamIH) else None
    if save:
        k = key(name, hist, iH, H_GeV, Aterm, lamIH, tag)
        np.savez_compressed(C.SP / f'{k}.npz', **{kk: (v if v is not None else np.array(np.nan))
                                                  for kk, v in ev.items() if not isinstance(v, dict)})
    return ev


def load(k):
    z = np.load(C.SP / f'{k}.npz', allow_pickle=True)
    ev = {kk: z[kk] for kk in z.files}
    for kk in ('tau_end_inf', 'tau0', 'xi', 'H_GeV', 'frag'):
        if kk in ev:
            v = float(ev[kk])
            ev[kk] = None if (kk == 'frag' and (not np.isfinite(v) or v == 0.0)) else v
    if 'm2xbar' in ev and ev['m2xbar'].shape == ():
        ev['m2xbar'] = None
    return ev


def Dgrid(ev, Ns=(60.0,), H_GeV=None, Mh=125.20, als=0.1180, Mts=MT_GRID, **kw):
    out = {}
    for Mt in Mts:
        smi = C.sm_in(Mt, Mh, als)
        b = C.danger(ev, smi, ev['Rbar'], Ns=Ns, m2xbar=ev.get('m2xbar'), H_GeV=H_GeV, **kw)
        out[Mt] = {N: b[N] for N in Ns}
    return out


def cmd_run(name, hist, iH=0):
    ev = run(name, hist, iH)
    Nob = 60.0 if name == 'M1r' else N_obs_track(name, 'G05' if hist == 'G05' else 'MD', iH)
    Ns = (60.0, round(Nob, 2))
    g = Dgrid(ev, Ns=Ns)
    Mts = list(g.keys())
    o = dict(name=name, hist=hist, iH=iH, xi_H=ev['xi'], N_obs_track=Nob, seconds=ev['seconds'], nstep=int(ev['nstep']),
             wronskian=[ev['wronskian_err0'], ev['wronskian_err_end']], runaway=bool(ev['runaway']),
             tau_end_inf=ev['tau_end_inf'],
             D60={str(m): g[m][60.0]['D'] for m in Mts}, DNobs={str(m): g[m][Ns[1]]['D'] for m in Mts},
             tau_c60={str(m): g[m][60.0].get('tau_c') for m in Mts},
             Mt_crit_60=C.crit(Mts, [g[m][60.0]['D'] for m in Mts]),
             Mt_crit_Nobs=C.crit(Mts, [g[m][Ns[1]]['D'] for m in Mts]))
    path = C.RES / 'sk_bench.json'
    try:
        allr = json.load(open(path))
    except Exception:
        allr = {}
    allr[f'{name}|{hist}|iH{iH}'] = o
    C.jdump(allr, path)
    print(name, hist, iH, 'xi %.5f' % ev['xi'], 'D60 172.57 %.1f 172.86 %.1f' % (o['D60']['172.57'], o['D60']['172.86']),
          'crit60 %s Nobs %s' % (o['Mt_crit_60'], o['Mt_crit_Nobs']), 'steps', o['nstep'], '%.0fs' % ev['seconds'], flush=True)
    return o


if __name__ == '__main__':
    c = sys.argv[1]
    if c == 'run':
        cmd_run(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 0)
