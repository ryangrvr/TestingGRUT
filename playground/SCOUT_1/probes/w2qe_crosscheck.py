"""W2-QE cross-check: forward particles vs backward-trajectory evaluation on the SAME 2-mode box (firewall setup)."""
import numpy as np
import w2qe_bohm_relaxation as F
import w2qe_firewall as B

rng = np.random.default_rng(2024)
M = 5
phases = rng.uniform(0, 2 * np.pi, (M, M))
amps = np.zeros((M, M)); amps[0, 1] = amps[1, 0] = 1.0
box = F.Box(M, amps, phases)
smooth = lambda x, y: (4 / np.pi ** 2) * np.sin(x) ** 2 * np.sin(y) ** 2
X0 = F.sample(smooth, 40000)
for C in (4, 8):
    print(f"forward particles, C={C}: t=0 Hbar = {F.hbar(box, X0, 0.0, C):.4f}")
out, X = F.evolve(box, X0, np.pi, 0.002, 8, nrec=2)
print("forward particles C=8 at t = 0, pi:", [round(h, 4) for _, h in out], " C=4 at pi:", round(F.hbar(box, X, np.pi, 4), 4))
Hf, Hb, Z = B.measures(box, smooth, np.pi, 0.002)
print("backward-trajectory at t = pi: norm", round(Z, 4), "H_fine", round(Hf, 4), "Hbar", {k: round(v, 4) for k, v in Hb.items()})
Hf, Hb, Z = B.measures(box, smooth, np.pi, 0.0005)
print("backward-trajectory dt=0.0005 at t = pi: norm", round(Z, 4), "H_fine", round(Hf, 4), "Hbar", {k: round(v, 4) for k, v in Hb.items()})
