"""Independent physical-matrix checks of the Gauss-law star decomposition.

This validates its factors on complete finite graph blocks, including the
degree-four central star. The analytical inequality itself has no electric
cutoff. No frozen Hamiltonian file is edited.
"""
from pathlib import Path
import itertools
import json
import math
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parents[1] / "constructive_physics/attempt12"))
from model import build, COST, ELECTRIC, HOP, MASS


def run(n, edges, L, a):
    pairs = list(itertools.combinations(range(n), 2))
    mask = sum(1 << pairs.index(tuple(sorted(e))) for e in edges)
    q = len(edges)
    H, states, _, _ = build(n, q, L, masks=[mask], rewire=False)
    incident = [[e for e, pair in enumerate(pairs)
                 if mask >> e & 1 and i in pair] for i in range(n)]
    removable = np.array([sum((bits >> (n + e) & 1) and
                                E[2*e] == E[2*e+1] == 0
                                for e in range(len(pairs)) if mask >> e & 1)
                          for _, bits, E in states])
    expected = H + diags(a * removable)
    stars = []
    rows = []
    coo = H.tocoo()
    for i in range(n):
        d = len(incident[i])
        occupancy = np.array([bits >> i & 1 for _, bits, _ in states])
        diag = np.array([
            COST*d/2 + 2*MASS*(bits >> i & 1)
            + ELECTRIC/4*sum(E[2*e]**2 + E[2*e+1]**2
                              for e in incident[i])
            + a/2*sum((bits >> (n+e) & 1) and
                       E[2*e] == E[2*e+1] == 0 for e in incident[i])
            for _, bits, E in states])
        delta = ELECTRIC/4 - a/2
        baseline = (COST + a)*d/2
        slack = diag - baseline - (2*MASS + delta)*occupancy
        assert min(slack) >= -1e-12
        # A frozen hopping entry changes exactly one vertex occupation.
        select = np.array([r != c and
                           ((states[r][1] ^ states[c][1]) >> i & 1)
                           for r, c in zip(coo.row, coo.col)], dtype=bool)
        hopping = coo_matrix((coo.data[select],
                            (coo.row[select], coo.col[select])),
                            shape=H.shape).tocsr()
        star = hopping + diags(diag)
        stars.append(star)
        zero = np.flatnonzero(occupancy == 0)
        one = np.flatnonzero(occupancy == 1)
        off = hopping[zero][:, one].toarray()
        norm = float(np.linalg.norm(off, 2)) if off.size else 0.
        assert norm <= HOP*math.sqrt(d) + 1e-12
        if star.shape[0] <= 200:
            minimum = float(np.linalg.eigvalsh(star.toarray())[0])
        else:
            minimum = float(eigsh(star, k=1, which="SA",
                                 return_eigenvectors=False, tol=1e-12)[0])
        D = 2*MASS + delta
        lower = baseline + (D - math.sqrt(D*D + 4*d*HOP*HOP))/2
        assert minimum >= lower - 1e-10
        rows.append(dict(vertex=i, degree=d, minimum_diagonal_slack=float(min(slack)),
                         hopping_block_norm=norm, hopping_norm_bound=HOP*math.sqrt(d),
                         actual_star_minimum=minimum, analytical_star_lower=lower))
    difference = sum(stars) - expected
    error = float(np.linalg.norm(difference.data))
    assert error < 1e-11
    return dict(N=n, Q=q, L=L, a=a, dimension=len(states),
                decomposition_error=error, stars=rows)


if __name__ == "__main__":
    controls = [
        run(5, [(0, i) for i in range(1, 5)], 3, .1),
        run(4, [(0, 1), (1, 2), (2, 3), (3, 0)], 3, .1),
        run(4, [(0, 1), (1, 2), (2, 3), (3, 0)], 3, .1*(6-4)/6),
    ]
    out = dict(status="GAUSS_STAR_LOWER_BOUND_CHECK_PASS", controls=controls,
               scope="Complete regulated graph blocks validate factors; the proved bound covers all integer electric fields.")
    (P / "star_bound_results.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
