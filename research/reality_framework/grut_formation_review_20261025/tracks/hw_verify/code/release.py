#!/usr/bin/env python3
"""hw_verify (Claude, own code): does the instant-release assumption matter more or less in the heavy regime?

Finite release duration: background 1 of hwrel.c.  At a = 1 the universe is inflaton-matter dominated with
H = H_i, the inflaton decays into radiation at width Gamma_phi, and the Hubble-induced X mass is tied to the
inflaton energy, -c_H (rho_phi/3)|X|^2, so it turns off over ~1/Gamma_phi (and the matter era lasts until
H ~ Gamma_phi).  Gamma_phi/H_i -> infinity is the instant release (background 0).  Fields start at the
inflationary minimum with zero velocity in both cases.
Late bookkeeping [D]: after reheating H = H_i sqrt(Rr/3)/a^2 with Rr = rho_r a^4 (H_i units; Rr = 3 for the instant
case), so the comoving entropy is (Rr/3)^{3/4} times the instant value (yield_diluted = yield / (Rr/3)^{3/4}) and the
decay-averaged <a> of the ledger is (Rr/3)^{1/4} times larger (R_central_reheat = R_central (Rr/3)^{1/4}).
Points: heavy spot (m = 30, Y-relaxed, delta = pi/2) and the stated regime (Codex candidate 10: stated couplings,
c_H = 10, H_i = 4.3487e11 GeV, delta = pi/6, Codex and Y-relaxed preparations).  a_end = 100, window 30-100.
Usage: python3 -B release.py ../results/release_<group>.json <group>   group = heavy | stated_codex | stated_relaxed
       (~2 min per group on 4 cores)
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hw

H = hw.H_CAND
GAMS = [None, 100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1]


def job(spec):
    name, cp, delta, prep, gam = spec
    t0 = time.time()
    z0, info = (hw.prep_relaxed if prep == 'relaxed' else hw.prep_codex)(cp, H, hw.CH, delta)
    kw = dict(bg=0) if gam is None else dict(bg=1, gam=gam)
    try:
        r = hw.run(cp, H, z0, a_end=100.0, windows=(1, 3, 10, 30, 100), **kw)
    except Exception as exc:
        return dict(name=name, delta=delta, prep=prep, gamma_over_H=gam, status=repr(exc))
    w = r['windows']['30-100']
    Rr = r['final'][-1]
    dil = (Rr / 3.0)**0.75
    o = hw.summarize(cp, H, w, prep_info=info,
                     extra=dict(name=name, delta=delta, prep=prep, gamma_over_H=gam, controls=r['end'],
                                Rr_end=Rr, rhophi_frac_end=r['final'][-2] / r['final'][-3]**3 / (r['final'][-2] / r['final'][-3]**3 + Rr / r['final'][-3]**4) if gam else 0.0,
                                entropy_dilution=dil,
                                windows_all={k: {q: v[q] for q in ('eps', 'NX', 'EX', 'NYH', 'NYL', 'Nsig', 'rhophi_frac', 'absX')}
                                             for k, v in r['windows'].items()},
                                history=[{q: s[q] for q in ('a', 'eps', 'NX', 'EX', 'absX', 'absY', 'rhophi_frac', 'H', 'NYH', 'NYL')} for s in r['samples']],
                                seconds=time.time() - t0))
    o['yield_diluted_over_YQ'] = o['yield_over_YQ'] / dil
    o['R_central_reheat'] = o['R_central'] * (Rr / 3.0)**0.25
    o['R_L_at_diluted_yield_note'] = 'R_L is the normalised-yield table value at the actual eps'
    print(name, prep, gam, 'eps %.5f NX %.1f EX %.4g NY %.4g R %.4g (reheat %.4g) RL %s Y/YQ %.3g (diluted %.3g) Rr %.3f %.0fs'
          % (o['eps'], o['NX'], o['EX'], o['NYH'] + o['NYL'], o['R_central'], o['R_central_reheat'], o['R_L'],
             o['yield_over_YQ'], o['yield_diluted_over_YQ'], Rr, o['seconds']), flush=True)
    return o


if __name__ == '__main__':
    out = Path(sys.argv[1]); group = sys.argv[2]
    if group == 'heavy':
        sp = [('heavy_spot', hw.heavy(), np.pi / 2, 'relaxed', g) for g in GAMS]
    elif group == 'stated_codex':
        sp = [('stated_cand10', dict(hw.STATED), np.pi / 6, 'codex', g) for g in GAMS]
    elif group == 'stated_relaxed':
        sp = [('stated_cand10', dict(hw.STATED), np.pi / 6, 'relaxed', g) for g in GAMS]
    else:
        raise SystemExit('group = heavy | stated_codex | stated_relaxed')
    rows = []
    with Pool(4) as p:
        for r in p.imap_unordered(job, sp):
            rows.append(r)
            out.write_text(json.dumps(dict(method=__doc__, H_GeV=H, cH=hw.CH, group=group, rows=rows), indent=1, default=float))
    print('rows', len(rows))
