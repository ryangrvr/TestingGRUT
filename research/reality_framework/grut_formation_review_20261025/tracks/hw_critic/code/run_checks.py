#!/usr/bin/env python3
"""hw_critic (Claude, own code): the release recomputations.
 A  spot point, Y-relaxed and Codex preparations (reproduction of ls_scaling / hw_* with this integrator);
 B  held-phase de Sitter relaxation from the Codex state for N e-folds, then release carrying the velocities
    (hw_verify item 1: R 385 / 88 / 7.9 at N = 0.5 / 1 / 2);
 C  the disputed number: free (unheld) de Sitter evolution from the Y-relaxed state for N e-folds, then release
    (i) carrying the acquired velocities (hw_verify's method) and (ii) with the velocities zeroed (closer to
    hw_gates' zero-velocity reading of the phase scan); records the charge density present at the end of the stage;
 D  finite release: inflaton decaying at Gamma_phi/H_i = 10, 1, 0.1 (hw_verify item 4);
 E  window control a_end = 100, window 30-100.
Usage: python3 -B run_checks.py ../results/checks.json      (~5 min on 4 cores)
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hwc

P = hwc.spot(); H = hwc.H_SPOT; D0 = np.pi / 2


def job(spec):
    kind, arg = spec
    t0 = time.time()
    extra = {}
    if kind == 'spot':
        z, info = hwc.prep(arg, D0, P); v = None; kw = {}
    elif kind == 'spot_a100':
        z, info = hwc.prep('relaxed', D0, P); v = None; kw = dict(a_end=100.0, win=(30.0, 100.0), nwin=281)
    elif kind == 'codex_held':
        z0, info = hwc.prep('codex', D0, P); z, v = hwc.desitter(z0, arg, P, held_delta=D0); kw = {}
    elif kind in ('free_vel', 'free_novel'):
        z0, info = hwc.prep('relaxed', D0, P); z, v = hwc.desitter(z0, arg, P, held_delta=None); kw = {}
        X, Y, S = z
        extra = dict(stage_end=dict(thX=float(np.angle(X)), thS=float(np.angle(S)), delta=float(np.angle(S) - 6 * np.angle(X)),
                                    absX=float(abs(X)), absY=float(abs(Y)), absS=float(abs(S)),
                                    nX_phys=float(-2 * np.imag(np.conj(X) * v[0])), nY_phys=float(-2 * np.imag(np.conj(Y) * v[1])),
                                    nS_phys=float(-2 * np.imag(np.conj(S) * v[2]))))
        if kind == 'free_novel':
            v = None
    elif kind == 'gamma':
        z, info = hwc.prep('relaxed', D0, P); v = None; kw = dict(bg=1, gam=arg)
    else:
        raise ValueError(kind)
    w = hwc.release(z, P, v0=v, **kw)
    o = hwc.summarize(P, w, H, tag='%s_%s' % (kind, arg), N_or_gamma=arg, seconds=time.time() - t0,
                      hist=w['hist'], window=w['window'], **extra)
    print('%-18s %-8s eps %+.5f NX %9.2f EX %9.1f NYH %9.4g NYL %9.4g Nsig %8.2g R %9.4g RL %s Y/YQ %.4g %4.0fs'
          % (kind, arg, o['eps'], o['NX'], o['EX'], o['NYH'], o['NYL'], o['Nsig'], o['R_central'],
             ('%.1f' % o['R_L_table']) if o['R_L_table'] else 'n/a', o['yield_over_YQ'], o['seconds']), flush=True)
    return o


if __name__ == '__main__':
    out = Path(sys.argv[1])
    jobs = [('spot', 'relaxed'), ('spot', 'codex'), ('spot_a100', None)]
    jobs += [('codex_held', N) for N in (0.5, 1.0, 2.0)]
    jobs += [('free_vel', N) for N in (0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0)]
    jobs += [('free_novel', N) for N in (0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0)]
    jobs += [('gamma', g) for g in (10.0, 1.0, 0.1)]
    rows = []
    with Pool(4) as pool:
        for o in pool.imap_unordered(job, jobs):
            rows.append(o)
            out.write_text(json.dumps(dict(method=__doc__, H_GeV=H, couplings=P, rows=rows), indent=1, default=float))
    print('done', len(rows))
