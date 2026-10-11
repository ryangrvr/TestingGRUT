"""Retarded-history determinant no-go and Gaussian action control.

No causal-set dynamics is proposed.  This script tests a tempting but invalid
history weight det(I + z G_R) and distinguishes it from a symmetric Euclidean
Gaussian matter action. Reproduces the 180-event comparison from the reviewed
TestingGRUT causal_state_probe work, using NumPy only.
"""
import numpy as np

SEED = 20261011
N = 180
FAKE_SEED = 71011


def four_dimensional_diamond(n, seed):
    rng = np.random.default_rng(seed)
    t_half = 0.5 * rng.random(n) ** 0.25
    t = np.where(rng.random(n) < 0.5, t_half, 1 - t_half)
    radius = np.minimum(t, 1 - t)
    vec = rng.normal(size=(n, 3))
    vec /= np.linalg.norm(vec, axis=1)[:, None]
    x = vec * (radius * rng.random(n) ** (1 / 3))[:, None]
    dt = t[:, None] - t[None, :]
    dx2 = np.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=-1)
    return ((dt > 0) & (dt * dt > dx2)).astype(float)


def links(C):
    return C * ((C @ C) == 0)


def compare():
    geom = four_dimensional_diamond(N, SEED)
    relations = int(geom.sum())
    assert relations == 1413, relations
    rng = np.random.default_rng(FAKE_SEED)
    block = np.zeros((N // 2, N // 2))
    block.flat[rng.choice(block.size, relations, replace=False)] = 1
    fake = np.zeros((N, N))
    fake[N // 2:, :N // 2] = block
    out = {}
    for name, C in (("4D_sprinkling", geom), ("height_2_fake", fake)):
        R = links(C)
        assert np.all(((C @ C) > 0) <= C.astype(bool))
        assert int(C.sum()) == relations
        assert np.max(np.abs(np.diag(R))) == 0
        # Acyclicity makes R nilpotent: the identity is exact for any real z.
        assert np.max(np.abs(np.diag(R @ R))) == 0
        assert np.linalg.norm(np.diag(R @ R @ R)) == 0
        for z in (0.013, -0.08, 0.4):
            sign, ld = np.linalg.slogdet(np.eye(N) + z * R)
            assert sign == 1 and abs(ld) < 5e-10, (name, z, sign, ld)
        # A quadratic real Euclidean action depends on the symmetric Hessian.
        S = (R + R.T) / 2
        eigen = np.linalg.eigvalsh(S)
        z = 0.01
        assert (1 + z * eigen.min()) > 0
        logdet = np.log1p(z * eigen).sum()
        sign, ld_direct = np.linalg.slogdet(np.eye(N) + z * S)
        assert sign == 1 and abs(logdet - ld_direct) < 1e-10
        # The universal leading small-z term follows from trace expansion.
        nlinks = int(R.sum())
        assert abs(float(np.trace(S @ S)) - nlinks / 2) < 1e-9
        leading = -z * z * nlinks / 4
        # Positive Euclidean Gaussian matter weighting det(K)^(-1/2).
        out[name] = {
            "relations": relations,
            "links": nlinks,
            "retarded_logdet": 0.0,
            "symmetric_logdet": float(logdet),
            "leading_small_z": float(leading),
            "gaussian_log_weight": float(-logdet / 2),
            "symmetric_min_eig_K": float(1 + z * eigen.min()),
        }
    # Positive control: a directed three-cycle is not a causal order.
    cycle = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    z = 0.3
    det_cycle = np.linalg.det(np.eye(3) + z * cycle)
    assert abs(det_cycle - (1 + z ** 3)) < 1e-13
    print("Directed noncausal cycle det(I+zR) =", det_cycle)
    for name, row in out.items():
        print(name, row)
    print("PASS: retarded determinant geometry-blind; symmetric action responds but is supplied")
    return out


if __name__ == "__main__":
    compare()
