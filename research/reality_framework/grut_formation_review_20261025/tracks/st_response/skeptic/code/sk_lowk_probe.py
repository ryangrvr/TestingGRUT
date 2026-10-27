#!/usr/bin/env python3
"""st_response skeptic: is the 2e-7..1e-6 low-k (0.1-0.8) difference between our e = 0 code and Codex's e = 1e-7
archive the gauge (longitudinal) term?  Runs st_response's OWN st_tangent.py (their code, used here as a tool, labelled
as such) at e = 0 and e = 1e-7 for k = 0.1, 0.2, 0.8 on (3, 0.5) to tau 800 on Codex's grid.  Usage: sk_lowk_probe.py <e>"""
import sys, json
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'code'))
import st_tangent as ST
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
C = np.load(SP / 'b73/research/reality_framework/grut_source_timing_20261026/quiet3.npz')
e = float(sys.argv[1]); ks = [0.1, 0.2, 0.8]
model, sol, info, sec = ST.run(3.0, 0.5, ks, C['times'], e=e)
ob = ST.observables(model, sol.y.T)
np.savez(SP / ('review/st_response_skeptic/lowk_e%g.npz' % e), times=sol.t, amp=ob['amp'], dN=ob['dN'], ks=np.array(ks))
print('done', e, sec)
