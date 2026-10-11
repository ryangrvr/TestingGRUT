"""CP09 bounded clock: classical velocity-mixture control, independent of grid code.

Rebuilds frozen N=2 CP09 Hamiltonian and probes independently; evaluates:
(1) fixed-speed exact finite source/pointer window protocol,
(2) exact weak-pulse spectral sinc^2 convolution,
(3) free Gaussian momentum distribution integrated by Gauss-Hermite,
(4) reported corrected bounded-clock number as a declared external datum.
Requires NumPy/SciPy; no modification to the frozen Hamiltonian.
"""
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import eigh, expm

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
K = np.kron
WIDTH = 0.15
TAU = 0.4
EPS = 0.001
ZETA = 0.001
MASS = 320
SIGMA_Q = 0.08
CORRECTED_BOUNDED_RESPONSE = 0.05368309506521327  # from 986394a, NOT rerun


def cp09_two_event():
    gam = math.pi / 4
    omega = .25
    u01 = math.cos(gam)*np.eye(4) - 1j*math.sin(gam)*K(Z, X)
    u10 = math.cos(gam)*np.eye(4) - 1j*math.sin(gam)*K(X, Z)
    H = np.zeros((12, 12), complex)
    onsite = omega*(2*np.eye(4) - K(X, I) - K(I, X))
    for j, u in [(1, u01), (2, u10)]:
        H[:4, :4] += np.eye(4)
        H[j*4:(j+1)*4, j*4:(j+1)*4] += np.eye(4)
        H[j*4:(j+1)*4, :4] -= u
        H[:4, j*4:(j+1)*4] -= u.conj().T
    for j in range(3):
        H[j*4:(j+1)*4, j*4:(j+1)*4] += onsite
    A, B = np.zeros_like(H), np.zeros_like(H)
    A[4:8, 4:8], A[8:12, 8:12] = K(I, Z), K(Z, I)
    B[4:8, 4:8], B[8:12, 8:12] = K(X, I), K(I, X)
    Es, V = eigh(H)
    assert abs(Es[0] - 0.09475528369770218) < 1e-10
    return H, V[:, 0], A, B, Es, V


def chiral_finite_pulse_response(v, H, g, A, B):
    """Classical constant speed v>0: source/detector exposure lasts width/v."""
    assert v > 0
    H0, source, detector = K(H, I), K(B, I), K(A, Y)
    pointer_x = K(np.eye(12), X)
    dur, between = WIDTH/v, (TAU-WIDTH)/v
    vals = []
    for sign in (+1, -1):
        state = K(g, np.array([1., 0.]))
        state = expm(-1j*(H0-sign*EPS/WIDTH*source)*dur) @ state
        state = expm(-1j*H0*between) @ state
        state = expm(-1j*(H0-ZETA/WIDTH*detector)*dur) @ state
        vals.append(float(np.vdot(state, pointer_x @ state).real))
    return -(vals[0]-vals[1])/(4*EPS*ZETA)


def spectral_weak_pulse_response(v, Es, V, A, B):
    """Independent exact leading-order-in-(EPS,ZETA) Kubo convolution.

    Clock traverses each spatial window in time width/v; each integrated
    strength is EPS/v or ZETA/v. Ground-state spectral kernel is multiplied
    by sinc^2[(En-E0)*width/(2v)] and clock-time separation tau/v.
    """
    energies = Es - Es[0]
    Aev = V.conj().T @ A @ V
    Bev = V.conj().T @ B @ V
    coefficients = Aev[0, :] * Bev[:, 0]
    arg = energies * WIDTH/(2*v)
    factor = np.sinc(arg/np.pi)**2
    correlation = np.sum(coefficients * np.exp(-1j*energies*TAU/v) * factor)
    return float(-2*correlation.imag/v**2)


def gaussian_momentum_average(fn, mass=MASS, n=16):
    """Free initial minimum-uncertainty packet: p~N(m,1/(4*sigma_q^2)).

    Only classical positive-speed clock trajectories are included. At m=320
    negative-speed Gaussian tails are negligible at machine precision.
    """
    x, w = hermgauss(n)
    sigma_v = 1/(2*SIGMA_Q*mass)
    v = 1 + np.sqrt(2)*sigma_v*x
    assert min(v) > 0
    return float(np.dot(w, [fn(float(vel)) for vel in v])/np.sqrt(np.pi))


def main():
    H, g, A, B, Es, V = cp09_two_event()
    U = expm(-1j*H*TAU)
    A_t = U.conj().T @ A @ U
    ideal = float((1j*np.vdot(g, (A_t@B-B@A_t)@g)).real)
    assert abs(ideal-.052588513538745)<1e-10
    chiral = lambda v: chiral_finite_pulse_response(v,H,g,A,B)
    spectral = lambda v: spectral_weak_pulse_response(v,Es,V,A,B)
    assert abs(chiral(1)-.052349192601293)<1e-10
    assert abs(spectral(1)-chiral(1)) < 6e-8
    for v in (0.96,0.98,1,1.02,1.04):
        assert abs(spectral(v)-chiral(v)) < 7e-8
    averages = [gaussian_momentum_average(chiral,n=n) for n in (8,12,16)]
    assert max(averages)-min(averages) < 1e-10
    classical = averages[-1]
    classical_spectral = gaussian_momentum_average(spectral,n=16)
    assert abs(classical_spectral-classical)<6e-8
    assert abs(classical-.05248912523974)<1e-9
    # Second-order classical velocity-spread control, error O(m^-4) for
    # a symmetric centered Gaussian under smooth response.
    h=1e-3
    curvature=(chiral(1+h)-2*chiral(1)+chiral(1-h))/h**2
    sigma_v=1/(2*SIGMA_Q*MASS)
    second_order_estimate=curvature*sigma_v**2/2
    assert abs((classical-chiral(1))-second_order_estimate)<4e-7
    corrected=CORRECTED_BOUNDED_RESPONSE
    budget=dict(
        ideal_instantaneous=ideal,
        matched_finite_width_chiral_at_v1=chiral(1),
        classical_gaussian_velocity_mixture=classical,
        classical_weak_pulse_spectral=classical_spectral,
        bounded_grid_published_corrected=corrected,
        finite_window_shift=chiral(1)-ideal,
        classical_velocity_mixture_shift=classical-chiral(1),
        residual_bounded_minus_classical=corrected-classical,
        net_bounded_minus_ideal=corrected-ideal,
        velocity_dispersion_std=sigma_v,
        response_velocity_curvature_at_v1=curvature,
        leading_small_dispersion_shift=second_order_estimate,
        important_boundary=("The bounded-grid result is imported; residual includes "
            "quantum packet coherences, recoil, finite-grid effects, circle "
            "boundary and incomplete clock passage; no unique backaction decomposition."),
    )
    assert abs(sum((budget['finite_window_shift'],
                    budget['classical_velocity_mixture_shift'],
                    budget['residual_bounded_minus_classical']))
               -budget['net_bounded_minus_ideal']) < 1e-12
    print(json.dumps(budget, indent=2))
    Path(__file__).with_name('results.json').write_text(
        json.dumps(budget, indent=2)+'\n')
    print('CP09_VELOCITY_ENSEMBLE_CONTROL_PASS')


if __name__=='__main__':
    main()
