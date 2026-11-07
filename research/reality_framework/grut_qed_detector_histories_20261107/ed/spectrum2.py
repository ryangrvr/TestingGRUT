import json, sys, time, os
import numpy as np
from scipy.sparse.linalg import eigsh
from schwinger import hamiltonian
cache = "spectrum_cache_m025.json"
res = json.load(open(cache)) if os.path.exists(cache) else {}
m = float(sys.argv[1]); xs = [float(v) for v in sys.argv[2].split(",")]; Ns = [int(v) for v in sys.argv[3].split(",")]
for x in xs:
    for N in Ns:
        key = f"m{m}_x{x}_N{N}"
        if key in res: continue
        t = time.time()
        H, L, st = hamiltonian(N, x, m)
        assert abs(H - H.T).max() < 1e-12
        w = eigsh(H, k=2, which="SA", tol=1e-10)[0]
        w.sort()
        omega = (w[1] - w[0]) / (2 * np.sqrt(x))
        res[key] = {"m": m, "x": x, "N": N, "E0_per_N": w[0] / N, "omega1": omega, "binding": omega - 2 * m}
        json.dump(res, open(cache, "w"), indent=1)
        print(key, round(omega - 2*m, 6), f"{time.time()-t:.1f}s", flush=True)
