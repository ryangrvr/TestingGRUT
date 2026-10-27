#!/usr/bin/env python3
"""st_response track (Claude, own code): run campaign for st_tangent.py.  Saves raw samples to the scratch directory.
Usage: python3 -B runs.py <job> [<job> ...]
Jobs:
  q3_k213      (3, 0.5), k = 2.13, rtol 1e-10, Codex's sample grid to tau 800   (required)
  q3_k3        (3, 0.5), k = 0.02, 2.13, 5.8 together, rtol 1e-10, to tau 800
  q3_k213_tight (3, 0.5), k = 2.13, rtol 1e-12, atol 1e-15, to tau 240            (own tolerance control)
  q3_k213_rk45  (3, 0.5), k = 2.13, RK45 (different method), rtol 1e-11, to tau 160 (own method control)
  q5_k2        (5, 1), k = 2.0, 2.13, rtol 1e-10, to tau 800
  q3_k213_e0   (3, 0.5), k = 2.13, e = 0 (gauge term off) to tau 240             (sensitivity, not a comparison)
  q3_gauge6    (3, 0.5), e = 1e-6, k = 0.02, 2.13, 5.8 to tau 800                 (vs Codex gauge3.npz)
  q3_k002_e0   (3, 0.5), e = 0, k = 0.02 to tau 800                               (size of the gauge term at small k)
  q5_k002      (5, 1), k = 0.02 to tau 800                                        (Codex's quiet5 10% crossing mode)
  q3_old       (3, 0.5), k = 2.122989207841491 on the old PHASE_PROFILES grid to 80 (old phase archive, sign convention)
The sample grid (times) is read from Codex's archived quiet3.npz only so the samples coincide for the comparison.
"""
import json, os, sys, time
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import st_tangent as S

SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
OUT = SP / 'review/st_response/runs'

JOBS = {
    'q3_k213': dict(b=3.0, gam=0.5, ks=[2.13], tmax=800.0),
    'q3_k3': dict(b=3.0, gam=0.5, ks=[0.02, 2.13, 5.8], tmax=800.0),
    'q3_k213_tight': dict(b=3.0, gam=0.5, ks=[2.13], tmax=240.0, rtol=1e-12, atol=1e-15),
    'q3_k213_rk45': dict(b=3.0, gam=0.5, ks=[2.13], tmax=160.0, rtol=1e-11, atol=1e-14, method='RK45'),
    'q5_k2': dict(b=5.0, gam=1.0, ks=[2.0, 2.13], tmax=800.0),
    'q3_k213_e0': dict(b=3.0, gam=0.5, ks=[2.13], tmax=240.0, e=0.0),
    'q3_gauge6': dict(b=3.0, gam=0.5, ks=[0.02, 2.13, 5.8], tmax=800.0, e=1e-6),
    'q3_k002_e0': dict(b=3.0, gam=0.5, ks=[0.02], tmax=800.0, e=0.0),
    'q5_k002': dict(b=5.0, gam=1.0, ks=[0.02], tmax=800.0),
    'q3_old': dict(b=3.0, gam=0.5, ks=[2.122989207841491], tmax=80.0, grid='old'),
}


def main(job):
    spec = dict(JOBS[job]); tmax = spec.pop('tmax')
    if spec.pop('grid', None) == 'old':
        grid = np.load(CODEX.parent / 'grut_reheating_modes_20261026' / 'PHASE_PROFILES.npz')['times']
    else:
        grid = np.load(CODEX / ('quiet3.npz' if spec['b'] == 3.0 else 'quiet5.npz'))['times']
    times = grid[grid <= tmax + 1e-9]
    t = time.time()
    model, sol, info, sec = S.run(spec.pop('b'), spec.pop('gam'), spec.pop('ks'), times, **spec)
    ob = S.observables(model, sol.y.T)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / (job + '.npz'), times=sol.t, Y=sol.y.T, ks=model.ks, **ob)
    meta = dict(job=job, spec=JOBS[job], seconds=sec, nfev=int(sol.nfev), status=sol.status,
                max_gauss=float(np.max(np.abs(ob['gauss']))),
                max_Haux_rel=float(np.max(np.abs(ob['Haux'] / ob['H'] - 1))),
                prep_residuals=[info['prep']['residual'], info['prep_plus_residual'], info['prep_minus_residual']])
    (OUT / (job + '.json')).write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    for j in sys.argv[1:]:
        main(j)
