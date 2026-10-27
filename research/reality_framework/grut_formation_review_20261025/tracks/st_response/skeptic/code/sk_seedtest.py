#!/usr/bin/env python3
"""st_response skeptic (Claude, own code): where does the ~5e-9 inter-code floor come from?

st_response attributes its 5e-9..8e-8 residual against Codex entirely to Codex's integration, because its own
tolerance/method controls move its answer by only ~1e-12.  Our third code (sk_tangent.py, also ~1e-12 by its own
tolerance control) differs from st_response's e = 0 run by ~5e-9 at the same e = 0.  Here we rerun our code at
k = 2.13 on (3, 0.5), e = 0, to tau 240 with three different seed vectors (ours, st_response's, Codex's archived)
and compare with st_response's e = 0 run (q3_k213_e0.npz).  The seeds are read as DATA (st_response's run file,
Codex's archive); no code of theirs is imported.  Writes ../results/sk_seedtest.json.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_tangent as S

SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
STR = SP / 'review/st_response/runs/q3_k213_e0.npz'


def main():
    R = np.load(STR); tR = R['times']
    q_st = R['Y'][0, 17:23]; dz_st = (q_st[0::2] + 1j * q_st[1::2]) / np.sqrt(2)
    C = np.load(CODEX / 'quiet3.npz'); qc = C['selected_vectors'][0, 0, :6]; dz_cx = (qc[0::2] + 1j * qc[1::2]) / np.sqrt(2)
    model = S.Model(3.0, 0.5, [2.13])
    y0, z0, dz_own, info = S.initial(model)
    out = dict(seed_diffs=dict(own_vs_st=float(np.max(np.abs(dz_own - dz_st))), own_vs_codex=float(np.max(np.abs(dz_own - dz_cx))),
                               st_vs_codex=float(np.max(np.abs(dz_st - dz_cx)))), runs={})
    ref = R['amp'][:, 0]; env = np.maximum.accumulate(np.abs(ref))
    C3 = C['amplitude_per_radian'][:, 9]; tC = C['times']; mC = tC <= 240.0
    envC = np.maximum.accumulate(np.abs(C3[mC]))
    for name, dz in (('own', dz_own), ('st_response', dz_st), ('codex', dz_cx)):
        y = y0.copy(); T = y[16:].reshape(1, 12); T[:, 0:6:2], T[:, 1:6:2] = dz.real, dz.imag; T[:, 6:] = 0
        te = np.union1d(tR, tC[mC])
        sol = solve_ivp(model.rhs, (0, 240.0), y, method='DOP853', t_eval=te, rtol=1e-11, atol=1e-12)
        ob = S.observables(model, sol.y.T)
        iR = np.searchsorted(sol.t, tR - 1e-12); iC = np.searchsorted(sol.t, tC[mC] - 1e-12)
        out['runs'][name] = dict(vs_st_response_e0=float(np.max(np.abs(ob['A'][iR, 0] - ref) / env)),
                                 vs_codex_e1e7=float(np.max(np.abs(ob['A'][iC, 0] - C3[mC]) / envC)))
        print(name, out['runs'][name], flush=True)
    (Path(__file__).resolve().parents[1] / 'results' / 'sk_seedtest.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
