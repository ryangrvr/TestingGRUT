"""Validation of the y = sqrt(C) integrator (ledger PM2-I2) against the registered explicit Euler in C (both with the
exact early stop PM2-I1, itself validated against full Euler).  Independent code path, not independent reviewer.
Usage: python pm2_yvar_check.py L seed (seed 0 = uniform start)."""
import sys, time
import numpy as np
import pm2_stageA
pm2_stageA.IMPL['max_rel_change'] = 0.05   # as run for the logged validation (before PM2-I2 production setting)
from pm2_stageA import lattice, adapt, adapt_y, DELTA, energy
L, seed = int(sys.argv[1]), int(sys.argv[2]); N, E = lattice(L); h = np.ones(N); h[0] = -(N - 1)
C0 = np.ones(len(E)) if seed == 0 else 1 + DELTA * np.random.default_rng(seed).uniform(-1, 1, len(E))
t = time.time(); Cy, Qy, sy, _ = adapt_y(N, E, h, C0); ty = time.time() - t
t = time.time(); Cc, Qc, sc, _ = adapt(N, E, h, C0); tc = time.time() - t
dif = int(np.sum((Cy > 0) != (Cc > 0)))
print(f"L={L} seed={seed} y: steps {sy} {ty:.0f}s | C: steps {sc} {tc:.0f}s | same topology {dif == 0} (edges differing {dif}) | "
      f"E_y {energy(Cy, Qy):.6e} E_C {energy(Cc, Qc):.6e}", flush=True)
