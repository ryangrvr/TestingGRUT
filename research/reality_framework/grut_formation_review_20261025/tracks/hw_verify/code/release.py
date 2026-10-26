#!/usr/bin/env python3
"""hw_verify (Claude, own code): does the instant-release assumption matter more or less in the heavy regime?

Finite release duration: background 1 of hwrel.c.  At a = 1 the universe is inflaton-matter dominated with
H = H_i, the inflaton decays into radiation at width Gamma_phi, and the Hubble-induced X mass is tied to the
inflaton energy, -c_H (rho_phi/3)|X|^2, so it turns off over ~1/Gamma_phi (and the matter era lasts until
H ~ Gamma_phi).  Gamma_phi/H_i -> infinity is the instant release (background 0).  Fields start at the
inflationary minimum with zero velocity in both cases.
Points: heavy spot (m = 30, Y-relaxed, delta = pi/2) and the stated regime (Codex candidate 10: stated couplings,
c_H = 10, H_i = 4.3487e11 GeV, delta = pi/6, Codex and Y-relaxed preparations).  a_end = 300, windows 100-300.
Usage: python3 -B release.py ../results/release.json     (~10 min on 4 cores)
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
        r = hw.run(cp, H, z0, a_end=300.0, windows=(1, 10, 30, 100, 300), **kw)
    except Exception as exc:
        return dict(name=name, delta=delta, prep=prep, gamma_over_H=gam, status=repr(exc))
    w = r['windows']['100-300']
    o = hw.summarize(cp, H, w, prep_info=info,
                     extra=dict(name=name, delta=delta, prep=prep, gamma_over_H=gam, controls=r['end'],
                                windows_all={k: {q: v[q] for q in ('eps', 'NX', 'EX', 'NYH', 'NYL', 'Nsig', 'rhophi_frac')}
                                             for k, v in r['windows'].items()},
                                history=[{q: s[q] for q in ('a', 'eps', 'NX', 'EX', 'absX', 'rhophi_frac', 'H')} for s in r['samples']],
                                seconds=time.time() - t0))
    print(name, prep, gam, 'eps %.5f NX %.1f EX %.4g R %.4g RL %s Y/YQ %.3g %.0fs' % (o['eps'], o['NX'], o['EX'], o['R_central'],
                                                                                   o['R_L'], o['yield_over_YQ'], o['seconds']), flush=True)
    return o


if __name__ == '__main__':
    out = Path(sys.argv[1])
    sp = [('heavy_spot', hw.heavy(), np.pi / 2, 'relaxed', g) for g in GAMS]
    sp += [('stated_cand10', dict(hw.STATED), np.pi / 6, 'codex', g) for g in GAMS]
    sp += [('stated_cand10', dict(hw.STATED), np.pi / 6, 'relaxed', g) for g in GAMS]
    rows = []
    with Pool(4) as p:
        for r in p.imap_unordered(job, sp):
            rows.append(r)
            out.write_text(json.dumps(dict(method=__doc__, H_GeV=H, cH=hw.CH, rows=rows), indent=1, default=float))
    print('rows', len(rows))
