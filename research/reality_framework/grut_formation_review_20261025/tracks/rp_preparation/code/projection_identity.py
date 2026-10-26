"""Exact projection identity behind the torque law (track rp_preparation, own code).

For a coherent state |alpha> and the charge projector P0 = (1/2pi) int dtheta exp(i theta Q),
any monomial O of annihilation operators with total charge zero satisfies O|alpha_theta> =
o(alpha)|alpha_theta> with the same eigenvalue for every theta, so
<alpha|P0 O P0|alpha> / <alpha|P0|alpha> = o(alpha) exactly. For O = c b^2 this gives
<c b^2>_projected = alpha_c alpha_b^2 = 2 sqrt2 e^{i delta}, hence
E(delta) = E_0 + 2 G2 Re(alpha_c alpha_b^2) = E_0 + 0.226274 cos(delta) and
T(delta) = i<[H, L]> = -12 G2 Im(alpha_c alpha_b^2) = -1.357645 sin(delta) = 6 dE/d delta,
independent of the projection (truncation apart). Mixed monomials such as n_k = a_k^dag a_k
are NOT preserved, so mean occupations do change (6, 2, 2 -> 5.83, 1.79, 1.86).
Usage: python3 projection_identity.py OUT.json
"""
import sys, json, math
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from neutral_boundary import Model, G2

out = []
for cuts in [(24, 12, 9), (36, 18, 15)]:
    m = Model(cuts)
    # c b^2 restricted to the neutral block: matrix elements from occupation numbers
    d = m.d
    look = {tuple(n): i for i, n in enumerate(m.nocc)}
    O = np.zeros((d, d))
    for j, (na, nb, nc) in enumerate(m.nocc):
        tgt = (na, nb - 2, nc - 1)
        if nb >= 2 and nc >= 1 and tgt in look:
            O[look[tgt], j] = math.sqrt(nc * nb * (nb - 1))
    rows = []
    for delta in (0.3, math.pi / 2, 2.0):
        psi, _, _ = m.state(delta)
        ev = np.vdot(psi, O @ psi)
        exact = 2 * math.sqrt(2) * np.exp(1j * delta)
        rows.append(dict(delta=delta, cb2_projected=[ev.real, ev.imag], cb2_coherent=[exact.real, exact.imag],
                         abs_error=float(abs(ev - exact)),
                         torque=m.torque(psi), torque_law=-12 * G2 * exact.imag,
                         n_a_projected=float(np.vdot(psi, np.diag(m.nocc[:, 0].astype(float)) @ psi).real), n_a_coherent=6.0))
    out.append(dict(cutoffs=list(cuts), neutral_dimension=d, rows=rows))
    print(cuts, [(round(r['delta'], 3), r['abs_error'], r['torque'], r['torque_law'], r['n_a_projected']) for r in rows])
json.dump(dict(scope=__doc__, models=out), open(sys.argv[1], 'w'), indent=1)
