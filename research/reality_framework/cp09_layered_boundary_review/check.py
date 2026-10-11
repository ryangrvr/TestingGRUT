"""CP-09 independent layered-family variational and boundary-count audit.

Verifies exact finite two-layer hypercube boundary enumeration, parallel
transport consistency, a normalized matter-product trial and its energy.
No CP-09 frozen source import or generating-law change. Requires numpy.
Run: python check.py
"""
import math
import numpy as np

GAMMA = math.pi / 4
OMEGA = 0.25
J = 1.0


def valid_poset(C):
    """Boolean transitivity and antisymmetry, independent of CP-09 code."""
    return (not np.any(C & C.T)
            and np.all(((C.astype(int) @ C.astype(int)) > 0) <= C))


def boundary_formula(a, b):
    return np.array([a * (a - 1) * (3/4)**b,
                     b * (b - 1) * (3/4)**a,
                     a * b * 2.**(-(a + b - 1))])


def spin_trial(B, gamma):
    """Product matter vector for bipartite relation matrix B.

    Lower spin i acquires exp(-i gamma Z_i) for every outgoing edge;
    upper spins stay in the X=+1 eigenstate. The initial lower Bloch
    phase -b*gamma optimizes the mean X expectation.
    """
    a, b = B.shape
    plus = np.array([1., 1.], complex) / math.sqrt(2)
    ans = np.array([1.], complex)
    for i in range(a):
        theta = gamma * (int(B[i].sum()) - b/2)
        spin = np.array([np.exp(-1j * theta),
                         np.exp(1j * theta)]) / math.sqrt(2)
        ans = np.kron(ans, spin)
    for _ in range(b):
        ans = np.kron(ans, plus)
    return ans


def test_partition(a, b, gamma=GAMMA):
    n = a + b
    num = 1 << (a * b)
    counts = np.zeros(3, dtype=int)
    mater_sum = 0.0
    worst_parallel_err = 0.
    Z = np.diag([1., -1.]).astype(complex)
    X = np.array([[0., 1.], [1., 0.]], complex)
    I = np.eye(2, dtype=complex)
    transports = {}
    for i in range(a):
        for j in range(a, n):
            P = np.array([[1.]], complex)
            for k in range(n):
                P = np.kron(P, Z if k == i else X if k == j else I)
            transports[i, j] = (math.cos(gamma)*np.eye(1 << n)
                                - 1j*math.sin(gamma)*P)
    for mask in range(num):
        B = np.array([(mask >> k) & 1 for k in range(a*b)],
                     dtype=int).reshape(a, b)
        C = np.zeros((n, n), dtype=int)
        C[:a, a:] = B
        assert valid_poset(C)
        trial = spin_trial(B, gamma)
        assert abs(np.vdot(trial, trial)-1) < 1e-12
        # Each internal hypercube edge is covariantly constant.
        for i in range(a):
            for j in range(b):
                if B[i, j]:
                    continue
                changed = B.copy()
                changed[i, j] = 1
                err = np.linalg.norm(
                    spin_trial(changed, gamma)-transports[i, j+a] @ trial)
                worst_parallel_err = max(worst_parallel_err, float(err))
        # Matter expectation: all upper spins remain |X+>, and lower
        # spins acquire relative phase 2 gamma (degree_i - b/2).
        mater_sum += sum(1-math.cos(2*gamma*(int(B[i].sum())-b/2))
                         for i in range(a))
        # External single-relation transitions, all orientations.
        for i in range(n):
            for j in range(n):
                if i == j or C[i, j] or (i < a and j >= a):
                    continue  # existing or still in the bipartite cube
                changed = C.copy()
                changed[i, j] = 1
                if valid_poset(changed):
                    if i < a and j < a:
                        counts[0] += 1
                    elif i >= a and j >= a:
                        counts[1] += 1
                    elif i >= a and j < a:
                        counts[2] += 1
                    else:
                        raise AssertionError((i, j))
    boundary = counts / num
    want_boundary = boundary_formula(a, b)
    matter_energy = OMEGA * mater_sum / num
    want_matter = OMEGA * a * (1-math.cos(gamma)**b)
    assert np.allclose(boundary, want_boundary, atol=2e-13)
    assert abs(matter_energy-want_matter) < 2e-12
    assert worst_parallel_err < 2e-12
    trial_energy = J*sum(boundary)+matter_energy
    return dict(a=a, b=b, orders=num, measured_boundary=boundary.tolist(),
                exact_boundary=want_boundary.tolist(),
                measured_matter_energy=matter_energy,
                exact_matter_energy=want_matter,
                variational_energy=trial_energy,
                max_internal_transport_error=worst_parallel_err)


def analytic_large_n(n):
    assert n % 2 == 0
    a = b = n // 2
    kin = J * sum(boundary_formula(a, b))
    mat = OMEGA*a*(1-math.cos(GAMMA)**b)
    return kin, mat, (kin+mat)/n


def run():
    for a, b in [(1, 1), (1, 2), (2, 2), (2, 3), (3, 3)]:
        result = test_partition(a, b)
        print('EXHAUSTIVE', result)
    print('BALANCED LARGE-N VARIATIONAL UPPER BOUNDS (not spectral ground energies)')
    for n in (10, 20, 30, 40, 50, 60, 100, 200):
        kin, mat, total_per_event = analytic_large_n(n)
        print(f'N={n:3d} kinetic={kin:.10f} matter={mat:.10f}'
              f' total/N={total_per_event:.10f}')
    kin, mat, per = analytic_large_n(200)
    assert abs(per-OMEGA/2) < 1e-7
    assert kin < 1e-6
    print('PASS: exponentially weak layered-hypercube boundary, exact'
          ' parallel transport, and O(Omega*N/2) trial energy')


if __name__ == '__main__':
    run()
