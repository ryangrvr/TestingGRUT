#!/usr/bin/env python3
"""ls_ledger (Claude, own code): homogeneous X, Y, Sigma releases for Codex b448f25 candidates 10 and 17.

Releases use our own sp_yledger integrator (tracks/sp_yledger/code/yrelease.c: cosmic time, physical fields,
Dormand-Prince 5(4); Codex uses conformal fields and DOP853), with our own multi-start inflationary minimum
(ymodel.minimum).  The binary is built in this track's scratch directory; nothing is written into sp_yledger.
Preparation (as stated in Codex 373d37e/b448f25): c_H = 10, H_i = 434867798206.6179 GeV, aligned radial minimum,
X > 0, Y < 0, Sigma = -|Sigma| exp(i delta), zero physical velocity, instant radiation release at a = 1.
Ensemble members: delta = delta0 + N(0, 0.001) drawn with numpy default_rng(2026102510 / 2026102517), which
reproduces the 48 published member phases exactly; we run members 0, 3, ..., 45 (16 per candidate) and the
nominal delta0.  Late state: window 600-1000 means from the C code plus our own Python snapshot projection
(ls_project.project) at a = 1000 and the Sigma radial decomposition (ls_project.sigma_decomposition).
Usage: python3 -B ls_release.py <out.json> [members_step] [extras]    (~10 min on 4 cores; extras adds a=3000 and
       rtol controls, ~10 more min)
"""
import json, os, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('SPY_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                      'scratchpad/review/ls_ledger/build')
sys.path.insert(0, str(HERE.parents[1] / 'sp_yledger' / 'code'))
sys.path.insert(0, str(HERE))
import ymodel as Y          # noqa: E402  our sp_yledger code
import ls_project as P      # noqa: E402

H = 434867798206.6179
CANDS = {10: dict(delta0=np.pi / 6, seed=2026102510), 17: dict(delta0=5 * np.pi / 24, seed=2026102517)}
WKEYS = ['eps', 'NX', 'Ncirc', 'EX', 'NYH', 'NYL', 'Nsig', 'JYcodex', 'NY', 'NS', 'Neff']


def job(spec):
    label, cand, member, delta, aend, rtol = spec
    z0, m = Y.initial(H, 10.0, delta)
    tend = 0.5 * aend**2
    r = Y.run(H, z0, t_end=tend, rtol=rtol, atol=rtol * 1e-3, nsamp=120)
    F = r['final']
    z = np.array([F[0] + 1j * F[1], F[2] + 1j * F[3], F[4] + 1j * F[5]])
    zd = np.array([F[6] + 1j * F[7], F[8] + 1j * F[9], F[10] + 1j * F[11]])
    snap = {k: (float(v) if np.ndim(v) == 0 else v) for k, v in P.project(tend, z, zd, H).items()}
    sig = P.sigma_decomposition(tend, z, zd, H)
    hist = [dict(a=s['a'], NYH=s['NYH'], NYL=s['NYL'], Nsig=s['Nsig'], eps=s['eps'], EX=s['EX']) for s in r['samples']]
    return dict(label=label, candidate=cand, member=member, delta=delta, a_end=aend, rtol=rtol,
                minimum=dict(x=float(m['x']), y=float(m['y']), s=float(m['s']), grad=m['grad_norm']),
                windows={k: {c: w[c] for c in WKEYS} for k, w in r['windows'].items()},
                snapshot=snap, sigma_decomposition=sig, history=hist, controls=r['end'],
                final_state=F[:12], t_end=tend)


def main(out, step=3, extras=False):
    specs = []
    for c, d in CANDS.items():
        deltas = d['delta0'] + np.random.default_rng(d['seed']).normal(0, 1e-3, 48)
        specs.append((f'{c}_nominal', c, -1, float(d['delta0']), 1000.0, 1e-11))
        for i in range(0, 48, step):
            specs.append((f'{c}_m{i}', c, i, float(deltas[i]), 1000.0, 1e-11))
        if extras:
            specs.append((f'{c}_nominal_a3000', c, -1, float(d['delta0']), 3000.0, 1e-11))
            specs.append((f'{c}_nominal_rtol1e-12', c, -1, float(d['delta0']), 1000.0, 1e-12))
            specs.append((f'{c}_m0_rtol1e-10', c, 0, float(deltas[0]), 1000.0, 1e-10))
    specs.sort(key=lambda s: -s[4])  # long runs first
    Y.build()
    with Pool(4) as p:
        res = p.map(job, specs)
    Path(out).write_text(json.dumps(dict(method=__doc__, H_GeV=H, runs=res), indent=1, default=float))
    for r in res:
        w = r['windows']['600-1000']
        print(r['label'], 'eps %.8f NYH %.4f NYL %.4f Nsig %.4f | snap NYH %.4f NYL %.4f Nsig %.4f free %.4f'
              % (w['eps'], w['NYH'], w['NYL'], w['Nsig'], r['snapshot']['NYH'], r['snapshot']['NYL'],
                 r['snapshot']['Nsig'], r['sigma_decomposition']['N_free']), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 3, len(sys.argv) > 3 and sys.argv[3] == 'extras')
