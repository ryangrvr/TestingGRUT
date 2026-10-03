"""Step-refinement check (ledger PM2-I3): final topology vs max_rel_change for the registered C-Euler (adapt) and the
y = sqrt(C) integrator (adapt_y), both with exact early stop.  Independent code path, not independent reviewer.
Usage: python pm2_dt_check.py L seed."""
import sys
import numpy as np
import pm2_stageA as M
L, seed = int(sys.argv[1]), int(sys.argv[2]); N, E = M.lattice(L); h = np.ones(N); h[0] = -(N - 1)
C0 = np.ones(len(E)) if seed == 0 else 1 + M.DELTA * np.random.default_rng(seed).uniform(-1, 1, len(E))
res = {}
for name, f, rels in [("C", M.adapt, [0.05, 0.025, 0.0125]), ("y", M.adapt_y, [0.05, 0.025, 0.0125, 0.00625])]:
    for r in rels:
        M.IMPL['max_rel_change'] = r; C = f(N, E, h, C0)[0]; res[(name, r)] = C > 0
ref = res[("y", 0.00625)]
print(f"L={L} seed={seed} edges differing from y@0.00625: " +
      ", ".join(f"{k[0]}@{k[1]}: {int(np.sum(v != ref))}" for k, v in res.items()), flush=True)
