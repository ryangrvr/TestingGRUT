#!/usr/bin/env python3
"""inf_reheat item 3 (Claude, own code): M1r release on slow-reheating histories with rprh.c.

Sets (all: M1r prepared state at u = 10 i, lambda = 0.5, portals vanish at M_Y, g_I re-fixed for each (b, wT) so that
mu_eff = 10 i H_i at the slice; released at rest at a = 1 with H = H_i):
  repro   (3, 0.5) to 1000: bit-for-bit check against dq_entropy / fx_holder.
  gscan   b = 3, Gamma in {0.3 ... 0.005}: full runs to radiation domination (end = max(1000, 25/Gamma)); yield, D, K.
  early   b = 3, Gamma in {1e-3, 1e-4, 0}: charge and energy to tau = 400 (matter era; the charge is fixed early).
  bscan   b in {1.5, 2, 4, 5}, Gamma = 1e-3, to tau = 400 (charge kick versus inflaton mass at slow decay).
  bhigh   b in {7.153, 11.427, 22, 33.185, 101.687} (the T-model masses at alpha = 0.03, 0.01, -, 1e-3, 1e-4), Gamma = 1e-3,
          to tau = 400, appended to bscan.
  bench   the self-consistent T-model benchmark points (b, Gamma_I) of benchmarks.json (alpha = 0.3, 1, 3, 10), to tau = 400.
  tmodel  b = 3, plateau wT in {1.0, 1.5, 2.5} (alpha = wT^2/6), Gamma = 1e-3, to tau = 400 (anharmonic first swing).
Usage: run_engine.py [set ...]   -> ../results/engine_runs.json (merged).   Labels: [C].
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rh_common as rc
import fx_common as fc

OUT = rc.RES / 'engine_runs.json'


def summarise(rr, lo_frac=0.4):
    end = rr['end']; lo = lo_frac * end
    NX = rc.window_mean(rr, 'NX', lo); EX = rc.window_mean(rr, 'EX', lo)
    NYL = rc.window_mean(rr, 'NYL', lo); NYH = rc.window_mean(rr, 'NYH', lo)
    L = rc.last(rr)
    D, K = rc.dilution(rr)
    y, yr, rr4 = fc.yield_over_target(rr, NX)
    eps = float(abs(NX) / (2 * np.sqrt(2 * 1e-3) * (EX / (3 * 1e-3)) ** 0.75))
    return dict(tag=rr['tag'], b=rr['b'], Gamma=rr['gam'], wT=rr['wT'], end=end, window=[lo, end], NX=NX, EX=EX, NYL=NYL, NYH=NYH,
                eps=eps, E_over_Ec=eps ** (-4 / 3), a_end=L['a'], H_end=L['H'], radfrac_end=L['radfrac'], rhoI_end=L['rho_I'],
                D=D, K=K, yield_over_YQ_at_HM1=yr, I0=rr['I0'], END=rr['END'], seconds=rr['seconds'], G=rr['G'])


def do(tag, b, gam, end, wT=0.0, nsamp=None):
    p, st, G, info = rc.prepare(0.5, b, wT)
    if nsamp is None:
        nsamp = int(round(end / 0.2))
    binw = end / 50.0
    # binw must be a multiple of end/nsamp
    nsamp = int(round(end / binw)) * max(1, int(round(nsamp / round(end / binw))))
    rr = rc.run(tag, p, b, gam, st['z'], G=G, end=end, nsamp=nsamp, binw=binw, wT=wT)
    s = summarise(rr)
    s['prep'] = dict(g_I=info['g_I'], alpha_I=info['alpha_I'], phi_slice_over_MP=info['I0'])
    # charge-fixing time: first tau at which |N_X| is within 10% / 1% of the late value
    S = rr['S']; tau = S[:, 0]; nx = np.abs(S[:, rc.COLS.index('NX')])
    for f in (0.1, 0.01):
        ok = np.abs(nx / abs(s['NX']) - 1) < f
        idx = np.where(~ok)[0]
        s['tau_fix_%g' % f] = float(tau[idx[-1] + 1]) if len(idx) and idx[-1] + 1 < len(tau) else float(tau[0])
    print('%-16s b=%.3g G=%.3g wT=%.3g NX=%.3f EX=%.1f eps=%.5f NYL=%.4f D=%.4f K=%.4f Y/YQ=%.4f rad=%.3g (%.1fs, %d steps)' % (
        tag, b, gam, wT, s['NX'], s['EX'], s['eps'], s['NYL'], s['D'], s['K'], s['yield_over_YQ_at_HM1'], s['radfrac_end'],
        s['seconds'], s['END']['steps']), flush=True)
    return s


def main():
    sets = sys.argv[1:] or ['repro', 'gscan', 'early', 'bscan', 'tmodel']
    res = json.load(open(OUT)) if OUT.exists() else {}
    if 'repro' in sets:
        p, st, G, info = rc.prepare(0.5, 3.0)
        rr = rc.run('repro_3_05', p, 3.0, 0.5, st['z'], G=G, end=1000.0)
        res['repro'] = dict(NX_400_1000=rc.window_mean(rr, 'NX', 400), NYL=rc.window_mean(rr, 'NYL', 400), NYH=rc.window_mean(rr, 'NYH', 400),
                            EX=rc.window_mean(rr, 'EX', 400), targets=dict(NX=-8045.485140302082, NYL=0.566640911258934,
                            NYH=0.010377131381878653, EX=188706.29275778163), D_K=rc.dilution(rr), END=rr['END'],
                            yield_over_YQ=fc.yield_over_target(rr, rc.window_mean(rr, 'NX', 400))[1])
        print('repro', res['repro'])
    if 'gscan' in sets:
        res['gscan'] = [do('g_3_%g' % g, 3.0, g, max(1000.0, 25.0 / g)) for g in (0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)]
    if 'early' in sets:
        res['early'] = [do('e_3_%g' % g, 3.0, g, 400.0) for g in (0.5, 1e-2, 1e-3, 1e-4, 0.0)]
    if 'bscan' in sets:
        res['bscan'] = [do('b_%g_1e-3' % b, b, 1e-3, 400.0) for b in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0)]
    if 'bhigh' in sets:
        res['bscan'] = [r for r in res.get('bscan', []) if r['b'] <= 5.0] + \
            [do('b_%g_1e-3' % b, b, 1e-3, 400.0) for b in (7.153, 11.427, 22.0, 33.185, 101.687)]
    if 'bench' in sets:
        bm = json.load(open(rc.RES / 'benchmarks.json'))
        pts = sorted({(round(b_['solution']['b'], 4), b_['solution']['Gamma_over_Hi']) for b_ in bm['benchmarks']
                      if b_['solution'] and b_['fam']['alpha'] in (0.3, 1.0, 3.0, 10.0)})
        res['bench_points'] = [do('bench_%g' % b, b, g, 400.0) for b, g in pts]
    if 'tmodel' in sets:
        res['tmodel'] = [do('t_3_w%g' % w, 3.0, 1e-3, 400.0, wT=w) for w in (1.0, 1.5, 2.5)]
    res['method'] = __doc__
    json.dump(res, open(OUT, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
