"""Extend the full-H ground-state test to a vanishing density of degree defects.

D(G) = sum_i |degree_i - 4|. Weighted Schur estimates control all rewiring
inside D<=K and across its boundary. A conservative matching-count ratio
bound permits finite, outward-rounded certificates without huge integers.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import itertools
import json
import math

import numpy as np
from scipy.sparse import coo_matrix

from variational_exclusion import regular_lower, frac_dec, sqrt_lower, matching_counts, finite_certificate

P = Path(__file__).resolve().parent


def graph_control(n=6):
    q = 2*n
    pairs = list(itertools.combinations(range(n), 2))
    M = len(pairs)
    masks = [sum(1 << e for e in es)
             for es in itertools.combinations(range(M), q)]
    index = {mask: i for i, mask in enumerate(masks)}
    defects = []
    for mask in masks:
        degree = [0]*n
        for e, (i, j) in enumerate(pairs):
            if mask >> e & 1:
                degree[i] += 1
                degree[j] += 1
        defects.append(sum(abs(d-4) for d in degree))
    defects = np.array(defects)
    rr, cc = [], []
    for i, mask in enumerate(masks):
        for e in range(M):
            if not mask >> e & 1:
                continue
            for f in range(M):
                if mask >> f & 1:
                    continue
                j = index[mask ^ (1 << e) ^ (1 << f)]
                assert abs(int(defects[i]-defects[j])) <= 4
                rr.append(i)
                cc.append(j)
    adjacency = coo_matrix((np.ones(len(rr)), (rr, cc)),
                           shape=(len(masks), len(masks))).tocsr()
    rows = []
    for K in [0, 2, 4, 6]:
        ids = np.flatnonzero(defects <= K)
        other = np.flatnonzero(defects > K)
        A = adjacency[ids][:, ids].toarray()
        B = adjacency[ids][:, other].toarray()
        internal_norm = float(max(abs(np.linalg.eigvalsh(A)))) if len(ids) else 0.
        cross_norm = float(np.linalg.norm(B, 2)) if B.size else 0.
        cross_count_bound = 5*(K+4)*(M-q) + q*(K+4)*(n-1)
        assert cross_norm**2 <= q*(M-q)*cross_count_bound + 1e-8
        incoming = np.sum(B, axis=0)
        assert max(incoming, default=0.) <= cross_count_bound
        # Test several independently supplied positive weighting functions.
        checks = []
        for eta in [.2, .5, .8]:
            weights = eta**defects[ids]
            weighted_rows = A @ weights / weights
            bound = eta**2*q*(M-q) + eta**(-4)*(5*K*(M-q)+q*K*(n-1))
            assert max(weighted_rows, default=0.) <= bound + 1e-8
            assert internal_norm <= bound + 1e-8
            checks.append(dict(eta=eta, actual_max_weighted_row=float(max(weighted_rows, default=0.)),
                               proven_row_bound=bound))
        rows.append(dict(K=K, sector_size=len(ids), internal_adjacency_norm=internal_norm,
                         crossing_adjacency_norm=cross_norm,
                         max_incoming_crossings=float(max(incoming, default=0.)),
                         proven_incoming_bound=cross_count_bound, weight_checks=checks))
    return dict(N=n, Q=q, graph_count=len(masks), rows=rows)


def conservative_trial_upper(n):
    """Replace the exact H_(k+1)/H_k by the proved lower bound N-3.

    H_0=1, H_1=N-1, H_k=(N-1)H_(k-1)-(k-1)H_(k-2).
    For N>=5, N-3 <= H_k/H_(k-1) <= N-1, k<=N.
    """
    assert n >= 5
    q = 2*n
    M = n*(n-1)//2
    center = n//4
    width = math.isqrt(n)
    ids = list(range(max(0, center-width), min(n, center+width)+1))
    w = [width+1-abs(k-center) for k in ids]
    norm = sum(x*x for x in w)
    diag = F(q, 2) + F(17, 8)*F(sum(k*x*x for k, x in zip(ids, w)), norm)
    with localcontext() as c:
        c.prec = 80
        c.rounding = ROUND_FLOOR
        total = Decimal(0)
        for j, k in enumerate(ids[:-1]):
            c2_lower = F(q-k, M-k)*(k+1)*(n-k)*(n-3)
            total += Decimal(w[j]*w[j+1])*sqrt_lower(c2_lower)
        subtraction = 2*total/Decimal(norm)
        with localcontext() as up:
            up.prec = 80
            up.rounding = ROUND_CEILING
            return frac_dec(diag, ROUND_CEILING) - subtraction


def certificate(n, K):
    q = 2*n
    M = n*(n-1)//2
    A1 = F(q*(M-q), 10*M)
    A2 = F(5*K*(M-q)+q*K*(n-1), 10*M)
    # Choosing eta using a float affects efficiency, not validity: every
    # positive rational eta<=1 gives a proved bound evaluated exactly below.
    eta = F(min(1., (2*float(A2/A1))**(1/6))).limit_denominator(1_000_000)
    internal = min(A1, A1*eta**2 + A2/eta**4)
    cross2 = F(q*(M-q)*(5*(K+4)*(M-q)+q*(K+4)*(n-1)), 100*M*M)
    upper = conservative_trial_upper(n)
    with localcontext() as c:
        c.prec = 80
        c.rounding = ROUND_FLOOR
        lower = Decimal(n)*regular_lower(n) - frac_dec(internal, ROUND_CEILING)
        gap = lower - upper
        if gap > 0:
            gap2 = gap*gap
            with localcontext() as up:
                up.prec = 80
                up.rounding = ROUND_CEILING
                probability = min(Decimal(1), frac_dec(cross2, ROUND_CEILING)/gap2)
        else:
            probability = Decimal(1)
        defect_density_lower = Decimal(K+1)*(1-probability)/Decimal(n)
    return dict(N=n, K=K, eta=str(eta), trial_energy_upper=str(upper),
                near_regular_energy_lower=str(lower), gap_lower=str(gap),
                internal_rewiring_norm_upper=str(frac_dec(internal, ROUND_CEILING)),
                cross_coupling_square_upper=str(frac_dec(cross2, ROUND_CEILING)),
                near_regular_ground_probability_upper=str(probability),
                mean_absolute_degree_defect_per_vertex_lower=str(defect_density_lower))


def asymptotic_defect_corollary():
    # g_infinity = 27/40 + sqrt(21/8) - sqrt(128321)/160.
    # These rational comparisons prove g>0.056 without float roots.
    assert F(81009, 50000)**2 < F(21, 8)
    assert F(17911, 50)**2 > 128321
    assert F(27, 40) + F(81009, 50000) - F(17911, 8000) > F(7, 125)
    # beta=3 Gamma 9^(1/3)<5/8, delta^(1/3)<47/1000.
    assert F(25, 12)**3 > 9
    delta = F(1, 10000)
    assert F(47, 1000)**3 > delta
    gap_lower = F(7, 125) - F(5, 8)*F(47, 1000)
    cross2_density_upper = F(9, 50)*delta
    assert gap_lower**2 > 39*cross2_density_upper
    # Hence limsup P(D<=floor(delta N))<1/40, and the mean defect
    # density is at least delta*(39/40) in the thermodynamic limit.
    return dict(defect_threshold_density=str(delta),
                probability_limsup_upper="1/40",
                mean_absolute_degree_defect_density_liminf_lower=str(delta*F(39, 40)),
                comparisons="Exact rational square/cube comparisons; the asymptotic limit follows from the proved trial and Schur bounds.")


if __name__ == "__main__":
    control = graph_control()
    ratio_checks = []
    for n in [5, 6, 16, 32, 64, 128]:
        counts = matching_counts(n)
        ratios = [F(counts[k+1], counts[k])*F(k+1, n-k)
                  for k in range(n)]
        assert all(n-3 <= r <= n-1 for r in ratios)
        if n >= 16:
            assert conservative_trial_upper(n) >= Decimal(finite_certificate(n)["rayleigh_upper_total"])
        ratio_checks.append(dict(N=n, min_ratio=float(min(ratios)), max_ratio=float(max(ratios))))
    certificates = [certificate(n, math.isqrt(math.isqrt(n)))
                    for n in [2**14, 2**16, 2**18, 2**20, 2**22, 2**24]]
    assert all(float(c["gap_lower"]) > 0 for c in certificates[1:])
    assert float(certificates[-1]["near_regular_ground_probability_upper"]) < .001
    out = dict(status="NEAR_REGULAR_GROUND_EXCLUSION_CHECK_PASS",
               full_graph_matrix_control=control, exact_matching_ratio_controls=ratio_checks,
               certificates=certificates,
               asymptotic_defect_corollary=asymptotic_defect_corollary(),
               scope="Full-H analytical bounds for D=sum_i|degree_i-4|; outward finite certificates have no electric or cluster cutoff.")
    (P / "near_regular_results.json").write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps(out, indent=2))
