"""Exact-formula and dense-matrix diagnostics of the declared new species.

These are transport controls, not formation simulations or empirical data.
The all-size counterexample is proved in MATTER_COUPLINGS.md.
"""
from fractions import Fraction as Q
from math import pi, sqrt
import numpy as np


def unitary(H, t):
    e, v = np.linalg.eigh(H)
    return (v * np.exp(-1j * e * t)) @ v.conj().T


def star(N, root, kappa):
    h = np.zeros((N, N))
    for j in range(N):
        if j != root:
            h[root, j] = h[j, root] = -kappa / sqrt(N - 1)
    return h


def run():
    rng = np.random.default_rng(901731)
    rows = []
    for N in [8, 16, 32, 64, 128, 256]:
        M, kappa = N - 1, .25
        jlocal = kappa / 1024
        t = pi / (2 * kappa)
        H = star(N, 0, kappa)
        U = unitary(H, t)
        initial = np.zeros(N); initial[:2] = 1 / sqrt(2)
        final = U @ initial
        ell, d = N // 4, N // 8
        far = np.arange(d, ell)
        # Small N=8 has d=1, which includes the other source; skip event formula there.
        if d >= 2:
            expected = .5 * (1 / M + 1 / M**2)
            assert np.max(abs(abs(final[far])**2 - expected)) < 2e-15
            response = np.imag(np.conj(U[far, 1]) * U[far, 0])
            assert np.max(abs(response + M**-1.5)) < 2e-15
            hstep = 2e-5
            def pulsed(theta):
                x = initial.astype(complex); x[0] *= np.exp(-1j * theta)
                return abs(U @ x)**2
            finite_response = (pulsed(hstep) - pulsed(-hstep)) / (2 * hstep)
            assert np.max(abs(finite_response[far] - response)) < 2e-11
        assert abs(np.linalg.norm(H, 2) - kappa) < 3e-15
        perm = rng.permutation(N)
        P = np.eye(N)[perm]
        newroot = int(np.where(perm == 0)[0][0])
        covariance = np.linalg.norm(P @ H @ P.T - star(N, newroot, kappa))
        assert covariance < 1e-14
        # Occupied bonds are a nonzero, fixed local contribution of the new species.
        m = max((N + 1) // 2, 2)
        local = np.zeros_like(H)
        for a in range(m - 1):
            local[a, a + 1] = local[a + 1, a] = -jlocal
        ulocal = unitary(H + local, t)
        deviation = np.linalg.norm((ulocal - U) @ initial)
        assert deviation <= 2 * jlocal * t + 1e-14
        pfar = float(sum(abs(final[far])**2))
        pfar_local = float(sum(abs((ulocal @ initial)[far])**2))
        if d >= 2:
            assert pfar > 1 / 32 and pfar_local > 1 / 64
        # Fully root-blind, permutation-invariant counterexample to exponential tails.
        allpairs = -kappa / M * (np.ones((N, N)) - np.eye(N))
        ua = unitary(allpairs, 1 / kappa)
        formula = np.exp(-1j / M) * (np.exp(1j * N / M) - 1) / N
        assert abs(ua[N - 1, 0] - formula) < 4e-15
        assert abs(ua[N - 1, 0]) >= .9 / N
        rows.append(dict(N=N, recorded_length=m, prefix_length=ell, minimum_recorded_distance=d,
                         star_norm=float(np.linalg.norm(H, 2)), relabel_covariance_error=float(covariance),
                         fixed_probe_time=t, star_far_set_probability=pfar,
                         star_plus_local_far_probability=pfar_local,
                         local_hopping=jlocal, state_deviation=deviation,
                         state_deviation_bound=2 * jlocal * t,
                         root_blind_offdiagonal_amplitude=float(abs(ua[N - 1, 0])),
                         all_size_geometry_motion_included='analytic Duhamel/prefix-tail bound, not this frozen-word matrix'))
    assert Q(1, 32) - Q(1, 128) - Q(1, 256) > Q(1, 64)
    assert Q(11, 10) + Q(1, 4) + Q(1, 2048) < Q(7, 5)
    assert Q(3) * Q(7, 5) * Q(8) < Q(128, 2)
    return dict(status='PASS', matrices=rows, pure_state_source_number_variance=.25,
                scope='new distinct-sector diagnostics; hypotheses and all-size proofs are separate; no WM2 Hamiltonian changed')
