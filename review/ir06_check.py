"""IR-05 follow-up: every case SCOUT-2 classified as 'incompatible' must be re-tested over the COMPLETE grouping class.
S2-SigmaH 'incompatible' cases: H1 (already IR-05), H5a (MFI mid-spectrum eigenstate), H7 (MFI Haar state).
Diagnostic file: rebuilds the ORIGINAL states by importing the original S2-SigmaH module (allowed for discrepancy
reconstruction); the compatibility test itself uses review code (nr_sigmah.py)."""
import sys, contextlib, io
import numpy as np
from nr_sigmah import *
sys.path.insert(0, "../probes/S2-SigmaH"); sys.path.insert(0, "../probes/S2-Sigma")
from nr2_diag import id_frame_min
with contextlib.redirect_stdout(io.StringIO()):
    import s2_sigma as SG
H = mfi(); ev, V = np.linalg.eigh(H)
# rebuild H7's Haar state exactly as s2_sigmah did: rng(1729), after the mean-field minimizations (6 BFGS restarts)
from scipy.optimize import minimize
rng = np.random.default_rng(1729)
def prod_state_a(angles):
    out = np.ones(1, complex)
    for a in angles: out = np.kron(out, np.array([np.cos(a / 2), np.sin(a / 2)]))
    return out
def energy(a): s = prod_state_a(a); return float(np.real(np.vdot(s, H @ s)))
mf = min((minimize(energy, rng.uniform(0, 2 * np.pi, n), method="BFGS") for _ in range(6)), key=lambda r: r.fun)
haar7 = (lambda v: v / np.linalg.norm(v))(rng.normal(size=N) + 1j * rng.normal(size=N))
eig5 = V[:, N // 2].astype(complex)
for nm, psi in [("H5a mid-spectrum eigenstate", eig5), ("H7 Haar state (original seed)", haar7)]:
    m, g, t = id_frame_min(psi)
    print(f"{nm:<32}: min over ALL 202 id-frame groupings of orbit-min C = {m:.4f} bits (grouping {g}, t = {t:+.2f}) -> "
          f"{'CERTIFIED incompatible in the H-local frame' if m > 0.1 else 'NOT certified (coarse-compatible)'}")
