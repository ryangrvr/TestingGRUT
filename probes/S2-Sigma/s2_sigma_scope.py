"""S2-Sigma scope check (Sigma-8): which assumptions buy uniqueness?
Restrict candidates step by step and test for a weakly dominant candidate.
  scope 1 (CPR-like): d = 2 fixed, n = 6 factors fixed, H-only criteria (L, P, M); frames: id, Cliffords, Haar,
          commutant frames exp(-iHs) for s in {0.3, 0.7, 1.9}
  scope 2: same, full (L, Q, R, P, M) with a state (Carroll-Singh-like: state-dependent)
  scope 3: d = 4 fixed (all 15 pairings) x frames, H-only
"""
import numpy as np
from scipy.linalg import expm
import s2_sigma as S

QUBITS = [[q] for q in range(S.n)]
PAIRINGS = [g for g in S.GROUPINGS if all(len(b) == 2 for b in g)]
print("=" * 78); print("S2-Sigma scope check: which assumptions buy uniqueness?"); print("=" * 78)
for mname, H in S.MODELS.items():
    frames = [("id", np.eye(S.N))] + [(f"cliff{i+1}", U) for i, U in enumerate(S.CLIFFS)] + [("haar", S.HAAR)] + \
             [(f"commutant s={s}", expm(-1j * H * s)) for s in (0.3, 0.7, 1.9)]
    st = S.states_for(H)["random local product"]
    st_t = S.evolve(st, H, 1.0)
    for scope, groupings, use_state in [("1: d=2, n=6, H-only", [QUBITS], False),
                                        ("2: d=2, n=6, H + state", [QUBITS], True),
                                        ("3: d=4 pairings, H-only", PAIRINGS, False)]:
        rows, meta = [], []
        for fname, U in frames:
            c = S.pauli_coeffs(U.conj().T @ H @ U)
            if use_state:
                E0 = S.entropies(S.in_frame(st, U)); Et = S.entropies(S.in_frame(st_t, U))
            for g in groupings:
                L, P, M = S.h_criteria(c, g)
                v = [L, P, M]
                if use_state:
                    v += list(S.state_criteria(E0, Et, g))
                rows.append(v); meta.append((fname, S.dtype_of(g), g))
        V = np.round(np.array(rows), 9)
        dom = S.dominant(V); front = S.pareto(V)
        print(f"  {mname[:28]:<28} scope {scope:<24}: candidates {len(V):3d}, front {len(front):3d}, "
              f"dominant set = {sorted({S.label(meta, i) for i in dom}) if len(dom) else 'NONE'}")
