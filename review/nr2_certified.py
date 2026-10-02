"""NR2 with CERTIFIED incompatibility: psi = U psi_prod with U a deep random Clifford, accepted only if, in the H-local
(id) frame, EVERY grouping has orbit-minimum total correlation Cmin >= 0.5 bit. The U-frame (where psi is exactly
product) is included as candidate 'conf'. Tests: dominance? both frames on the front? rules A vs B disagree?"""
import numpy as np
from nr_sigmah import *
from nr2_diag import id_frame_min
H = mfi(); p = prod_state(1.4)
res = []
for seed in (101, 202, 303):
    rng = np.random.default_rng(seed); tries = 0
    while True:
        tries += 1; U = rand_clifford(rng, depth=40); psi = U @ p
        m, g, t = id_frame_min(psi)
        if m >= 0.5: break
    C = [rand_clifford(rng) for _ in range(2)]; Hu = haar(rng)
    frames = [("id", np.eye(N, dtype=complex)), ("conf", U), ("cl2", C[0]), ("cl3", C[1]), ("haar", Hu),
              ("comm", expm(-1j * H * 0.7)), ("Wpsi", product_frame(psi, rng).conj().T)]
    R = evaluate(H, psi, frames); front, dom = pareto(R, "Cmin"); cls = lu_classes(front, R, frames)
    fr = sorted({R[c[0]]["frame"] for c in cls})
    A = sorted({R[i]["frame"] for i in lex(R, "L", "Cmin")}); B = sorted({R[i]["frame"] for i in lex(R, "C", "Cmin")})
    print(f"seed {seed}: accepted after {tries} tries; id-frame min Cmin = {m:.3f} bits; front {len(front)} -> {len(cls)} LU classes;"
          f" frames {fr}; dominant {'yes' if dom else 'NONE'}; rule A {A} vs rule B {B}")
    res.append((not dom) and A != B)
print(f"NR2-certified: {'REPRODUCED' if all(res) else 'PARTIAL' if any(res) else 'NOT REPRODUCED'} ({sum(res)}/3): "
      "no dominance and priority-order disagreement for certified-incompatible (H-local frame vs psi-product frame) pairs")
