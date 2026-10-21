"""Supplement to item 3/4: with p(+) = p(-) = 1/2 per correlation volume (symmetry of the ensemble), do both vacua
percolate, so that the initial wall configuration contains walls that span the box (not only closed bubbles)?
Uncorrelated site model on a simple-cubic lattice with lattice spacing = freeze-out correlation length, open boundaries,
6-connectivity.  Literature: p_c(site, sc) = 0.3116077(2) (Wang et al. 2013, arXiv:1302.0421, abstract) < 1/2.
Writes results/percolation.json.  Deterministic seeds."""
import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from common import to_builtin

OUT = Path(__file__).resolve().parent.parent / "results"


def analyse(L, seed):
    rng = np.random.default_rng(seed)
    plus = rng.random((L, L, L)) < 0.5
    out = {"L": L, "seed": seed, "fraction_plus": float(plus.mean())}
    labs = {}
    for name, mask in (("plus", plus), ("minus", ~plus)):
        lab, n = ndimage.label(mask)
        sizes = np.bincount(lab.ravel())
        sizes[0] = 0
        big = int(np.argmax(sizes))
        span = [bool(np.isin(big, np.take(lab, 0, axis=a)).item() and np.isin(big, np.take(lab, L - 1, axis=a)).item()) for a in range(3)]
        out[name] = {"n_clusters": int(n), "largest_fraction_of_sign": float(sizes[big] / mask.sum()),
                     "largest_spans_xyz": span, "second_largest_fraction": float(np.sort(sizes)[-2] / mask.sum())}
        labs[name] = (lab, big)
    # wall bonds: nearest-neighbour pairs with opposite sign
    lp, bp = labs["plus"]
    lm, bm = labs["minus"]
    n_wall = 0
    n_inf = 0
    for a in range(3):
        s0 = [slice(None)] * 3
        s1 = [slice(None)] * 3
        s0[a] = slice(0, L - 1)
        s1[a] = slice(1, L)
        A, B = plus[tuple(s0)], plus[tuple(s1)]
        opp = A != B
        n_wall += int(opp.sum())
        # identify plus side and minus side labels of each wall bond
        lpA, lpB = lp[tuple(s0)], lp[tuple(s1)]
        lmA, lmB = lm[tuple(s0)], lm[tuple(s1)]
        plus_lab = np.where(A, lpA, lpB)
        minus_lab = np.where(A, lmB, lmA)
        n_inf += int((opp & (plus_lab == bp) & (minus_lab == bm)).sum())
    out["wall_bonds"] = n_wall
    out["fraction_of_wall_area_between_the_two_largest_clusters"] = n_inf / n_wall
    return out


def main():
    t0 = time.time()
    runs = [analyse(L, 1000 + 10 * L + k) for L in (32, 64, 128) for k in range(3)]
    res = {"model": "uncorrelated sites, p = 1/2, simple cubic, open boundaries, 6-connectivity",
           "p_c_site_sc_literature": 0.3116077, "runs": runs,
           "summary": {L: {"min_largest_fraction_plus": min(r["plus"]["largest_fraction_of_sign"] for r in runs if r["L"] == L),
                           "min_largest_fraction_minus": min(r["minus"]["largest_fraction_of_sign"] for r in runs if r["L"] == L),
                           "all_span": all(all(r[s]["largest_spans_xyz"]) for r in runs if r["L"] == L for s in ("plus", "minus")),
                           "mean_wall_area_fraction_between_largest": float(np.mean([r["fraction_of_wall_area_between_the_two_largest_clusters"] for r in runs if r["L"] == L]))}
                       for L in (32, 64, 128)},
           "caveat": "Initial-condition geometry only.  Says nothing about subsequent coarsening, scaling, expansion, "
                     "backreaction or persistence; those belong to the network modules.",
           "seconds": time.time() - t0}
    (OUT / "percolation.json").write_text(json.dumps(to_builtin(res), indent=1))
    print(json.dumps(to_builtin(res["summary"]), indent=1), "time %.1fs" % res["seconds"])


if __name__ == "__main__":
    main()
