"""Diagnose NR2: is the 'conflict' pair actually incompatible? For psi_c = C psi_prod, compute, in the H-local (id)
frame, the minimum over groupings of Cmin (orbit-minimum total correlation), and which groupings reach it."""
import numpy as np, sys
from nr_sigmah import *
H = mfi(); p = prod_state(1.4); ev, V = np.linalg.eigh(H)
def id_frame_min(psi):
    traj = [V @ (np.exp(-1j * ev * t) * (V.conj().T @ psi)) for t in TG]
    ent = {}
    best = (9, None, None)
    for g in GROUPS:
        Ct = []
        for k, st in enumerate(traj):
            tot = 0
            for b in g:
                key = (k, tuple(b))
                if key not in ent: ent[key] = block_entropies(st, [b])[0]
                tot += ent[key]
            Ct.append(tot)
        m = min(Ct)
        if m < best[0] - 1e-9: best = (m, g, TG[int(np.argmin(Ct))])
    return best
if __name__ == '__main__':
    sys.path.insert(0, "../probes/S2-Sigma")
    import s2_sigma as ORIG    # only to rebuild the ORIGINAL conflict state C_orig psi_prod for comparison
    cands = [("original S2-SigmaH conflict (CLIFFS[0])", ORIG.CLIFFS[0] @ p)]
    for seed in (11, 22, 33):
        rng = np.random.default_rng(seed); C = rand_clifford(rng); cands.append((f"fresh seed {seed} cl1", C @ p))
    for nm, psi in cands:
        m, g, t = id_frame_min(psi)
        print(f"{nm:<42}: min over id-frame groupings of Cmin = {m:.4f} bits at t = {t:+.2f}, grouping {g}")
