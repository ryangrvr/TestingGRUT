"""Independent finite-pulse and quantum-clock resource controls for frozen CP-09.

Use only the supplied N=2 CP-09 matrices; no imported project code.
This does NOT simulate the full bounded-clock lattice calculation.
Dependencies: numpy, scipy. Run python check.py.
"""
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, expm

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def cp09_two_event():
    gamma = math.pi / 4
    Omega = 1 / 4
    U01 = math.cos(gamma) * np.eye(4) - 1j * math.sin(gamma) * np.kron(Z, X)
    U10 = math.cos(gamma) * np.eye(4) - 1j * math.sin(gamma) * np.kron(X, Z)
    H = np.zeros((12, 12), dtype=complex)
    onsite = Omega * (2 * np.eye(4) - np.kron(X, I) - np.kron(I, X))
    for k, U in ((1, U01), (2, U10)):
        H[0:4, 0:4] += np.eye(4)
        H[4*k:4*k+4, 4*k:4*k+4] += np.eye(4)
        H[4*k:4*k+4, 0:4] -= U
        H[0:4, 4*k:4*k+4] -= U.conj().T
    for k in range(3):
        H[4*k:4*k+4, 4*k:4*k+4] += onsite
    A = np.zeros_like(H)
    B = np.zeros_like(H)
    A[4:8, 4:8] = np.kron(I, Z)
    A[8:12, 8:12] = np.kron(Z, I)
    B[4:8, 4:8] = np.kron(X, I)
    B[8:12, 8:12] = np.kron(I, X)
    e, v = eigh(H)
    assert abs(e[0] - 0.09475528369770218) < 1e-11
    assert np.linalg.norm(H - H.conj().T) < 1e-13
    return H, v[:, 0], A, B


def chiral_square_width_response(H, ground, A, B, width, eps=0.001, zeta=0.001):
    # Both interactions autonomous in the translation-clock model P; this
    # evaluates its exact characteristics, not the dispersive P^2/2m evolution.
    tau = 0.4
    H0 = np.kron(H, I)
    K = np.kron(B, I)
    F = np.kron(A, Y)
    M = np.kron(np.eye(H.shape[0]), X)
    vals = []
    for sign in (1, -1):
        psi = np.kron(ground, np.array([1, 0]))
        psi = expm(-1j * (H0 - sign * eps / width * K) * width) @ psi
        psi = expm(-1j * H0 * (tau - width)) @ psi
        psi = expm(-1j * (H0 - zeta / width * F) * width) @ psi
        vals.append(float(np.vdot(psi, M @ psi).real))
    return -(vals[0] - vals[1])/(4 * eps * zeta)


def free_gaussian_parameters(mass, sigma=0.08, source_flight=0.7, detector_flight=1.1, width=0.15):
    # 1D minimal-uncertainty, initially unchirped Gaussian,
    # <p>=m (speed 1), sigma_p=1/(2 sigma_q), hbar=1.
    sp = 1/(2*sigma)
    Emean = mass/2 + sp**2/(2*mass)
    Evar = sp**2 + sp**4/(2*mass**2)
    st_source = math.hypot(sigma, source_flight/(2*mass*sigma))
    st_detector = math.hypot(sigma, detector_flight/(2*mass*sigma))
    # Operational estimate, NOT an error bound for the detector response:
    # Fourier W_s ~ exp(i k q) with k=1/width and W_d analogously.
    eta = (detector_flight-source_flight)/(mass*width**2)
    weyl_commutator_norm = 2*abs(math.sin(eta/2))
    # At this same detector flight time, optimizing sigma_q for an initially
    # unchirped minimum-uncertainty Gaussian gives sigma(t)^2 >= t/m.
    min_detector_sigma = math.sqrt(detector_flight/mass)
    return dict(mass=mass, mean_kinetic=Emean, std_kinetic=math.sqrt(Evar),
                sigma_source=st_source, sigma_detector=st_detector,
                min_detector_sigma_unchirped=min_detector_sigma,
                gate_fourier_eta=eta, gate_weyl_commutator_norm=weyl_commutator_norm)


def test():
    H, g, A, B = cp09_two_event()
    U = expm(-1j*H*.4)
    At = U.conj().T @ A @ U
    ideal = float((1j*np.vdot(g,(At@B-B@At)@g)).real)
    assert abs(ideal-0.05258851353874493)<1e-11
    widths = [0.2, 0.15, 0.1, 0.05, 0.025]
    chiral = {str(w): chiral_square_width_response(H,g,A,B,w) for w in widths}
    assert abs(chiral['0.15']-0.05234919260)<2e-9
    assert abs(chiral['0.2']-0.052159928)<2e-7
    bounded = 0.0539630583 # published 320/2048/dt=.00125 run, not rerun here
    gap_pulse = chiral['0.15']-ideal
    gap_clock = bounded-chiral['0.15']
    assert abs((bounded-ideal)-(gap_pulse+gap_clock)) < 1e-14
    rows=[free_gaussian_parameters(m) for m in (20,80,320)]
    for r in rows:
        mass=r['mass']
        assert r['sigma_detector']>=r['min_detector_sigma_unchirped']-1e-12
        # Both propagating-clock unitaries are noncommuting Weyl controls,
        # except in the infinite-mass/classical limit.
        assert 0 < r['gate_weyl_commutator_norm'] < 2
        assert abs(r['mean_kinetic']-(mass/2+19.53125/mass)) < 1e-11
    assert rows[0]['gate_weyl_commutator_norm'] > rows[1]['gate_weyl_commutator_norm'] > rows[2]['gate_weyl_commutator_norm']
    result={'ideal':ideal,'chiral_finite_width_same_eps_zeta_as_bounded':chiral,
        'bounded_m320_reported':bounded,'finite_width_shift':gap_pulse,
        'additional_finite_mass_shift_vs_chiral':gap_clock,
        'net_error':bounded-ideal,'clock_gaussian_controls':rows,
        'limitations':'No reexecution of bounded-clock grid simulation; analytic free-clock packet and exact chiral finite-width comparator only.'}
    print(json.dumps(result,indent=2))
    print('CP09_CLOCK_ERROR_BUDGET_INDEPENDENT_CHECK_PASS')
    return result


if __name__=='__main__':
    test()
