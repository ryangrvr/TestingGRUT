"""Check of the Li HDE CPL projection quoted by the others track (w_a > 0 for c = 1)."""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
out = {}
for c, Om0 in ((1.0, 0.70), (1.0, 0.73), (0.6, 0.70)):
    f = lambda x, O: [O[0] * (1 - O[0]) * (1 + 2 * math.sqrt(max(O[0], 0)) / c)]
    s = solve_ivp(f, (0, math.log(0.3)), [Om0], dense_output=True, rtol=1e-10, atol=1e-12)
    a = np.linspace(0.3, 1, 200)
    O = s.sol(np.log(a))[0]
    w = -1 / 3 - 2 * np.sqrt(O) / (3 * c)
    A = np.vstack([np.ones_like(a), 1 - a]).T
    w0, wa = np.linalg.lstsq(A, w, rcond=None)[0]
    out[f"c={c},Omega_de0={Om0}"] = {"w(a=1)": float(w[-1]), "CPL_w0": float(w0), "CPL_wa": float(wa)}
json.dump(out, open(sys.argv[1], "w"), indent=1); print(json.dumps(out, indent=1))
