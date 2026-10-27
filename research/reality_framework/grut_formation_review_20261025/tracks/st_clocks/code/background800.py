#!/usr/bin/env python3
"""st_clocks (Claude, own code): extend OUR reheating engine (tracks/rp_reheating/code/rp.py + rpreheat.c, Claude's own
engine that reproduced Codex e160356) to tau = H_i t = 800 on the two quiet backgrounds (m_I, Gamma_I)/H_i = (3, 0.5),
(5, 1), at the nominal rtol 1e-10 and a tight control rtol 1e-12.  No Codex code is imported.

Output: $SP/review/st_clocks/rp/runs/<tag>.npy (all 26 diagnostic columns at Delta tau = 0.05) and
results/background800.json (late C, matching, checks, sample a values).
Usage: RP_SCRATCH=$SP/review/st_clocks/rp python3 -B background800.py
"""
import json, os, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRACKS = HERE.parents[1]
sys.path.insert(0, str(TRACKS / 'rp_reheating' / 'code'))
os.environ.setdefault('RP_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                      'scratchpad/review/st_clocks/rp')
import rp  # noqa: E402

CASES = [('q3', 3.0, 0.5), ('q5', 5.0, 1.0)]


def main():
    out = {}
    for tag, b, g in CASES:
        for label, rtol, atol in [('nominal', 1e-10, 1e-13), ('tight', 1e-12, 1e-15)]:
            name = '%s_%s_800' % (tag, label)
            res, S = rp.run(name, mode=0, b=b, gam=g, end=800.0, rtol=rtol, atol=atol, nsamp=16000)
            c = {k: i for i, k in enumerate(rp.COLS)}
            tau, a, H = S[:, c['tau']], S[:, c['a']], S[:, c['H']]
            C = a * a * H
            sel = tau >= 64.0
            out[name] = dict(b=b, gamma_I=g, rtol=rtol, final_a=float(a[-1]), final_H=float(H[-1]),
                             C_at_800=float(C[-1]), C_min_64_800=float(C[sel].min()), C_max_64_800=float(C[sel].max()),
                             late_window_NX_640_800=res['late']['NX'], NX_at_800=float(S[-1, c['NX']]),
                             NYL_late=res['late']['NYL'], NYH_late=res['late']['NYH'], EX_late=res['late']['EX'],
                             radfrac_final=float(S[-1, c['radfrac']]), checks=res['checks'], seconds=res['seconds'],
                             npy=str(rp.SCR / 'runs' / (name + '.npy')))
            print(name, json.dumps({k: v for k, v in out[name].items() if k != 'npy'}), flush=True)
    (HERE.parent / 'results').mkdir(exist_ok=True)
    (HERE.parent / 'results' / 'background800.json').write_text(json.dumps(out, indent=2) + '\n')


if __name__ == '__main__':
    main()
