import json, os, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
from schwinger_mps import gap
m = float(sys.argv[1]); cache = f"production_m{m}.json"
res = json.load(open(cache)) if os.path.exists(cache) else {}
for x in [9, 16, 36, 64]:
    for L in [12, 16, 20, 24]:
        N = int(round(L * np.sqrt(x))); N += N % 2
        key = f"x{x}_L{L}"
        if key in res: continue
        t = time.time(); om, info = gap(N, float(x), m, return_info=True)
        res[key] = {"m": m, "x": x, "L": L, "N": N, "omega1": om, **{k: float(v) for k, v in info.items()}, "sec": round(time.time() - t, 1)}
        json.dump(res, open(cache, "w"), indent=1); print(key, N, res[key]["sec"], flush=True)
# convergence check at largest point
t = time.time(); om160 = gap(192, 64.0, m, chi=(40, 80, 120, 160))
res["chi160_check_x64_L24"] = {"omega1_chi160": om160, "delta_vs_chi120": om160 - res["x64_L24"]["omega1"], "sec": round(time.time() - t, 1)}
json.dump(res, open(cache, "w"), indent=1); print("DONE", flush=True)
