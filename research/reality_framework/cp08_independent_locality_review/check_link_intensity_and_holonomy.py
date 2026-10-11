"""Independent CP-08 locality boundary controls (no new physical law).

Checks an exact Poisson-sprinkling link expectation in unbounded 4D
Minkowski space, plus the CP-08 flat-layer/nonflat-chain transport algebra.
Requires numpy and scipy. Does not edit CP-07/08 or Stage-3.
"""
import math
import numpy as np
from scipy.integrate import quad


def expected_4d_links(rapid_max: float, density: float) -> float:
    """Palm conditional expectation given event at origin, all tau > 0.

    Future y=(tau cosh chi, tau sinh chi n). Its Alexandrov volume
    from the origin is pi*tau**4/24. Link probability = exp(-rho V).
    dVol = tau**3 sinh(chi)**2 d tau d chi dOmega_2.
    """
    c4 = math.pi / 24
    radial, re = quad(lambda tau: density * tau**3 *
                      math.exp(-density*c4*tau**4), 0, math.inf,
                      epsabs=2e-12, epsrel=2e-12)
    angular, ae = quad(lambda chi: math.sinh(chi)**2,
                       0, rapid_max, epsabs=2e-12, epsrel=2e-12)
    assert re < 1e-9 and ae < 1e-9
    return 4*math.pi*radial*angular


def transport_algebra():
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], complex)
    Z = np.diag([1., -1.]).astype(complex)
    def op(a, b, c):
        return np.kron(np.kron(a, b), c)
    p = op(Z, X, I)  # 0 -> 1
    q = op(I, Z, X)  # 1 -> 2; composable with p
    r = op(Z, I, X)  # 0 -> 2; height-two pair with p
    assert np.linalg.norm(p@q + q@p) < 1e-12
    assert np.linalg.norm(p@r - r@p) < 1e-12
    gamma = math.pi / 4
    def U(P):
        return math.cos(gamma)*np.eye(8)-1j*math.sin(gamma)*P
    A, B, R = U(p), U(q), U(r)
    nonflat = np.trace(A.conj().T@B.conj().T@A@B).real/8
    flat = np.trace(A.conj().T@R.conj().T@A@R).real/8
    assert abs(nonflat-(1-2*math.sin(gamma)**4)) < 1e-12
    assert abs(flat-1.) < 1e-12
    return nonflat, flat


def run():
    print('4D future links, rapidity truncated but proper time unbounded:')
    for rho in (0.1, 1., 10.):
        for cutoff in (0.5, 1., 2., 3.):
            numerical = expected_4d_links(cutoff, rho)
            exact = 6*math.sinh(2*cutoff)-12*cutoff
            assert abs(numerical-exact) < 2e-9*max(1, exact)
            print(f'rho={rho:>4g} rapidity<={cutoff:g}: '
                  f'numerical={numerical:.10f} exact={exact:.10f}')
    assert 6*math.sinh(20)-120 > 1e9
    nonflat, flat = transport_algebra()
    print('CP08 nontrivial 2-step transport holonomy:', nonflat)
    print('CP08 height-two transport holonomy:', flat)
    print('PASS: density-independent 4D link divergence; holonomy local to composable arrows')


if __name__ == '__main__':
    run()
