"""st_physics skeptic, check k5 (own code, [C] on the declared couplings, [E] for its consequence).

The declared light sector (WORKING_THEORY.md: V >= g_X^2 phi^2 |X|^2 + lambda_X |X|^4 + lambda_phi (phi^2 - v^2)^2 / 4,
g_X = 0.1, lambda_X = 1e-3, v_phi = 0.924 GeV << H_i) couples the neutral field phi to X.  Codex's six-field mode system
(X, Y, Sigma) does not contain phi fluctuations; phi = 0 is only pinned homogeneously.  On the matched conformal
X orbit, the conformal phi mode (real field, conformal variables, R = 0 in radiation, a v_phi negligible) obeys
    dphi'' + (k^2 + g^2 r(eta)^2) dphi = 0,          r = |chi_X|,  m_phi^2 = 2 g^2 |X|^2 -> g^2 r^2 conformally,
a Hill equation with the radial period T.  We compute its Floquet exponent mu_phi(k) per conformal time on the same
matched orbits as k1 and compare with the X self-resonance.  Also: the growth needed for the phi-induced X mass
g^2 <dphi^2> to reach the quartic lambda r^2 at pericentre, starting from a vacuum-like dphi.

usage: python3 -I k5_phi_resonance.py <k1_corot_floquet.json> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
sys.path.insert(0, __import__('os').path.dirname(__file__))
from k1_corot_floquet import CoRot  # noqa: E402

G2 = 0.1 ** 2
k1f, out = sys.argv[1:3]
k1 = json.load(open(k1f))
res = dict(model='V >= g^2 phi^2 |X|^2, g = 0.1 (declared); conformal phi mass^2 = g^2 r^2; Hill equation on the matched orbit')


def mu_phi(orb, k, g2=G2):
    def rhs(t, y):
        r, rp = y[0], y[1]
        m2 = k * k + g2 * r * r
        Y = y[2:].reshape(2, 2)
        return np.r_[rp, orb.L ** 2 / r ** 3 - orb.lam * r ** 3, Y[1], -m2 * Y[0]]
    y0 = np.r_[orb.rmin, 0.0, np.eye(2).ravel()]
    s = solve_ivp(rhs, (0, orb.T), y0, method='DOP853', rtol=1e-12, atol=1e-13)
    M = s.y[2:, -1].reshape(2, 2)
    tr = np.trace(M)
    if abs(tr) <= 2:
        return 0.0, float(tr)
    m = (abs(tr) + np.sqrt(tr * tr - 4)) / 2
    return float(np.log(m) / orb.T), float(tr)


for name in ('quiet3_late', 'quiet5_late', 'quiet3_tau10'):
    E, L = k1[name]['E'], k1[name]['L']
    orb = CoRot(E, L)
    grid = np.r_[0.0, np.arange(0.05, 8.001, 0.05)]
    mus = np.array([mu_phi(orb, k)[0] for k in grid])
    j = int(np.argmax(mus))
    lo, hi = max(grid[j] - 0.05, 0.0), grid[j] + 0.05
    pk = minimize_scalar(lambda k: -mu_phi(orb, k)[0], bounds=(lo, hi), method='bounded', options=dict(xatol=1e-5))
    bands = []
    on = mus > 1e-9
    i = 0
    while i < len(grid):
        if on[i]:
            j0 = i
            while i + 1 < len(grid) and on[i + 1]:
                i += 1
            bands.append(dict(k_from=float(grid[j0]), k_to=float(grid[i]), mu_max=float(mus[j0:i + 1].max())))
        i += 1
    # g^2 sensitivity of the peak
    gs = {}
    for g in (0.05, 0.07, 0.1, 0.14):
        mm = [mu_phi(orb, k, g * g)[0] for k in np.arange(0.0, 8.001, 0.1)]
        gs[str(g)] = float(max(mm))
    X_peak = k1[name]['peak']['mu']
    row = dict(rmin=orb.rmin, rmax=orb.rmax, T_eta=orb.T, phi_mass2_range=[G2 * orb.rmin ** 2, G2 * orb.rmax ** 2],
               mu_phi_k0=float(mus[0]), mu_phi_peak=dict(k=float(pk.x), mu=float(-pk.fun)), bands_mu_gt_0=bands[:12],
               mu_phi_peak_vs_g=gs, X_self_resonance_peak=X_peak,
               ratio_phi_over_X=float(-pk.fun / X_peak))
    # [E] growth needed: induced X mass^2 g^2 <dphi_c^2> vs lam r_min^2 at pericentre.
    # vacuum <dphi_c^2> per ln k ~ k^2 / (4 pi^2) at the resonant k (massless normalization, order-of-magnitude)
    kp = float(pk.x) if pk.x > 0.05 else 1.0
    var0 = kp ** 2 / (4 * np.pi ** 2)
    need = 1e-3 * orb.rmin ** 2 / G2
    row['estimate'] = dict(vacuum_var_per_lnk=var0, var_needed_for_g2phi2_eq_lam_rmin2=need,
                           amplitude_growth_needed=float(np.sqrt(need / var0)),
                           d_eta_needed=float(np.log(np.sqrt(need / var0)) / max(-pk.fun, 1e-9)))
    res[name] = row
    print(name, json.dumps(row)[:1500], flush=True)
json.dump(res, open(out, 'w'), indent=1)
