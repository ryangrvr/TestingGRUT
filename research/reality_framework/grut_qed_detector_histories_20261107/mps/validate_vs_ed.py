import json, sys, time
sys.path.insert(0, "../qed")
import numpy as np
from schwinger_mps import gap
ed0 = json.load(open("../qed/spectrum_cache.json")); ed1 = json.load(open("../qed/spectrum_cache_m025.json"))
cases = [(0.0, 4.0, 16, ed0), (0.0, 16.0, 20, ed0), (0.25, 8.0, 18, ed1), (0.25, 32.0, 24, ed1)]
out = []
for m, x, N, ed in cases:
    t = time.time(); om, info = gap(N, x, m, return_info=True)
    ref = ed[f"m{m}_x{x}_N{N}"]["omega1"]
    out.append({"m": m, "x": x, "N": N, "dmrg": om, "ed": ref, "E0_per_N_dmrg": info["E0"]/N, "E0_per_N_ed": ed[f"m{m}_x{x}_N{N}"]["E0_per_N"], "rel": abs(om - ref) / ref, "sec": round(time.time() - t, 1)})
    print(out[-1], flush=True)
ok = all(o["rel"] < 1e-5 for o in out)
json.dump({"gate": "rel<1e-5 vs exact diagonalization", "pass": ok, "cases": out}, open("validation_vs_ed.json", "w"), indent=1)
print("GATE", "PASS" if ok else "FAIL")
