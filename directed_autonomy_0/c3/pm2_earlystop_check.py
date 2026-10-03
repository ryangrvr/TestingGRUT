"""Validation of the exact spanning-tree early stop (ledger PM2-I1) against the full explicit-Euler run (registered
tolerances, no early stop).  Independent code path, not independent reviewer.  Usage: python pm2_earlystop_check.py L seed
(seed 0 = uniform start)."""
import sys, time
import numpy as np
import pm2_stageA
pm2_stageA.IMPL['max_rel_change'] = 0.05   # as run for the logged validation (before PM2-I2 production setting)
from pm2_stageA import lattice, adapt, DELTA, energy, tree_from_edges
L, seed = int(sys.argv[1]), int(sys.argv[2]); N, E = lattice(L); h = np.ones(N); h[0] = -(N - 1)
C0 = np.ones(len(E)) if seed == 0 else 1 + DELTA * np.random.default_rng(seed).uniform(-1, 1, len(E))
t = time.time(); Ce, Qe, se, re = adapt(N, E, h, C0, True); te = time.time() - t
t = time.time(); Cf, Qf, sf, rf = adapt(N, E, h, C0, False); tf = time.time() - t
same = np.array_equal(Ce > 0, Cf > 0); relC = np.max(np.abs(Ce - Cf)) / Cf.max()
print(f"L={L} seed={seed} early: steps {se} {te:.0f}s | full: steps {sf} resid {rf:.1e} {tf:.0f}s active {(Cf>0).sum()} "
      f"tree {tree_from_edges(N,E,Cf>0)[2] and (Cf>0).sum()==N-1} | same topology {same} | max|dC|/maxC {relC:.1e} | "
      f"dE/E {abs(energy(Ce,Qe)-energy(Cf,Qf))/energy(Cf,Qf):.1e}", flush=True)
