"""Exact/bounded ensemble test for the *supplied* symmetric link-Gaussian action.

No physical causal-set path measure or quantum history dynamics is derived.
Run from this directory: python gaussian_ensemble_entropy_control.py
Requires numpy and the preserved supplied_probe.py.
"""
import math
import numpy as np
from supplied_probe import four_dimensional_diamond, links, SEED


def gaussian_log_weight(C, z):
    """Positive real Gaussian det(I+z*A_Hasse/2)^(-1/2)."""
    R = links(C)
    A = R + R.T
    assert np.array_equal(A, A.T) and np.all(np.diag(A) == 0)
    eig = np.linalg.eigvalsh(A / 2)
    k = 1 + z * eig
    assert min(k) > 0, "Gaussian kinetic operator not positive definite"
    return float(-0.5 * np.log(k).sum()), int(R.sum())


def bipartite_causal_set(B):
    a, b = B.shape
    C = np.zeros((a + b, a + b))
    C[a:, :a] = B.T
    return C


def test():
    n = 180
    alpha = 1.8
    z = alpha / n
    geom = four_dimensional_diamond(n, SEED)
    count = int(geom.sum())
    assert count == 1413
    rng = np.random.default_rng(71011)
    random_B = np.zeros((90, 90))
    random_B.flat[rng.choice(random_B.size, count, replace=False)] = 1
    random_layered = bipartite_causal_set(random_B)
    full_layered = bipartite_causal_set(np.ones((90, 90)))
    bound = alpha**2 / (32 * (1 - alpha/2))
    for name, C in (
        ("4D_sprinkling", geom),
        ("matched_random_2_layer", random_layered),
        ("complete_balanced_2_layer", full_layered),
    ):
        logw, e = gaussian_log_weight(C, z)
        assert 0 <= e <= n*n//4
        assert -1e-10 <= logw <= bound + 1e-10
        print(name, "relations=", int(C.sum()), "links=", e,
              "log_weight=", format(logw, ".12f"),
              "weight=", format(math.exp(logw), ".12f"))
    full_logw, full_e = gaussian_log_weight(full_layered, z)
    assert full_e == n*n//4
    # Complete bipartite Hasse spectrum of A/2: +/-n/4 and n-2 zeros.
    full_exact = -0.5*math.log(1 - alpha**2/16)
    assert abs(full_logw-full_exact) < 1e-12

    # Exhaust all 2^(3*3) labelled bipartite orders, fixed partition.
    n_small = 6
    z_small = alpha/n_small
    weights = []
    for mask in range(2**9):
        B = np.array([(mask >> j) & 1 for j in range(9)],
                     dtype=float).reshape((3, 3))
        logw, _ = gaussian_log_weight(bipartite_causal_set(B), z_small)
        weights.append(math.exp(logw))
    assert len(weights) == 512
    assert min(weights) >= 1 - 1e-12
    assert max(weights) <= math.exp(bound) + 1e-12
    assert abs(max(weights) - math.exp(full_exact)) < 1e-12
    print("enumerated_two_layer_N6:", len(weights),
          "sum_weights=", format(sum(weights), ".12f"),
          "mean_weight=", format(np.mean(weights), ".12f"))
    print("uniform_logweight_bound=", format(bound, ".12f"),
          "uniform_max_weight=", format(math.exp(bound), ".12f"))
    print("N180_bipartite_log_configuration_count=",
          format((n*n//4)*math.log(2), ".8f"))
    # z fixed at 0.01 eventually fails positivity for complete 2-layer orders.
    fixed_z = 0.01
    assert 1-fixed_z*180/4 > 0 and 1-fixed_z*500/4 < 0
    print("fixed_z=0.01_complete_layered_positive_at_N180_but_not_N500: PASS")
    print("PASS: bounded weights cannot reverse a vanishing uniform fraction")


if __name__ == "__main__":
    test()
