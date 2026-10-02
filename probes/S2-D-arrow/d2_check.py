"""D2 large-sample check: the Haar measure on an energy shell is invariant under U_t and K, so
P(dS > 0) = P(dS < 0) exactly; test the 150-sample 0.61/0.39 split."""
import contextlib, io
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import s2_darrow as M   # reruns the main script silently to reuse its model
rng = np.random.default_rng(7)
sel = np.abs(M.E) < 1.0
inc = []; d = []
for _ in range(1500):
    v = rng.normal(size=sel.sum()) + 1j * rng.normal(size=sel.sum()); st = M.V[:, sel] @ v / np.linalg.norm(v)
    a = M.S_red(st, M.n, M.SYS); b = M.S_red(M.ev(st, 1.0), M.n, M.SYS); d.append(b - a)
d = np.array(d)
print(f"energy shell, 1500 samples: <dS> = {d.mean():+.5f} +- {d.std()/np.sqrt(len(d)):.5f}; P(increase) = {np.mean(d>0):.3f} +- {np.sqrt(.25/len(d)):.3f}")
