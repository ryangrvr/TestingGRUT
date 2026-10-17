"""Small algebra and exact constants for the receiver-dependence lemma."""
from fractions import Fraction as Q
import numpy as np
from charge_controls import Blocks


def run():
    assert Q(160, 69) * Q(33, 80)**128 < Q(1, 2**120)
    assert (4 * Q(9, 10)**2 - 1) > Q(6, 5)**2 * Q(11, 10)**4
    # exp(486/125)>4 suffices; a positive rational partial sum proves it.
    assert 1 + Q(486, 125) > 4
    rates = []
    for a in [.9, 1., 1.1]:
        for omega in [-.25, 0., .25]:
            for delta in [-.75, -.17, .41, .75]:
                gamma = np.sqrt(4 * a*a - (omega + delta)**2) / a**2
                assert gamma > 1.2
                rates.append(float(gamma))
    flat = []
    for N in [3, 4, 5, 6]:
        b = Blocks(N - 1, 1, np.linspace(-.75, .75, N - 1), 1.1)
        states, _ = b.basis(0)
        K = np.array([len(w) for w, n in states])
        diagonal = np.array([sum(b.delta[v] + 2 for v in w) for w, n in states])
        H = b.H(0) + np.diag(2 * K)
        V = H - np.diag(diagonal)
        vnorm = float(np.linalg.norm(V, 2))
        assert vnorm <= 1.1 + 1e-14
        assert np.all(diagonal >= 1.25 * K - 1e-14)
        assert H[0, 0] == 0
        flat.append(dict(N=N, dimension=len(H), interaction_norm=vnorm,
                         norm_bound=1.1, diagonal_lower_per_attachment=1.25,
                         initial_mean_energy=0, evolution='not run; energy inequality is analytic'))
    return dict(status='PASS', finite_tail_constant_exact=True,
                uniform_rate_bound_exact='6/5', sampled_minimum_rate=min(rates),
                flat_receiver_operator_checks=flat,
                scope='sufficient-criterion diagnostics; no physical-L trajectory or reduced receiver claim')
