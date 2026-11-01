"""Track m1_vacuum, item 2 (stability): fixed-charge radial Hessian of the Q = 200000 M1 object.
Independent code (Claude, 1 Nov 2026).  The object is re-solved on each grid with qball_m1 (Numerov-Newton),
then the second-variation operators are built with the 3-point Laplacian on u = rho dF, w = rho dS:
    L_+ = [[-d^2 + S^2 + 3aF^2 - Om^2, 2FS], [2FS, -d^2 + D]],  D = U_phi''(S) + F^2
    K   = 4 Om^2 |rho F><rho F| / int rho^2 F^2        (fixed-charge rank-one term, derived in NOTES.md)
    L_- = -d^2 + S^2 + aF^2 - Om^2                        (phase sector; rho F is its exact zero mode)
    ell=1 blocks: L_+ + 2/rho^2                           (translation zero mode)
Dirichlet at rho = 0 and rho = L.  The point of the box-size scan: the continuum of free X quanta starts at
m_X^2 - omega^2 = 1 - Om^2; in a box of radius L its lowest level sits near 1 - Om^2 + (pi/L)^2.
"""
import json, math, sys
import numpy as np
from scipy.linalg import eigh

sys.path.insert(0, sys.argv[1])
import qball_m1 as q

a, b, c = q.a, q.b, q.c
out = {'boxes': []}
seed = None
for (L, h) in [(40.0, 0.1), (60.0, 0.1), (80.0, 0.1), (80.0, 0.05), (120.0, 0.1), (200.0, 0.1), (300.0, 0.1)]:
    sv, u, w, Om, ob = q.solve_fixed_Q(200000, L, h, seed=seed)
    seed = (sv, u, w, Om)
    r = sv.r
    n = sv.n
    F = u / r
    S = 1 + w / r
    D = q.Uphi_pp(S) + F * F
    lap = (np.diag(2 * np.ones(n)) - np.diag(np.ones(n - 1), 1) - np.diag(np.ones(n - 1), -1)) / h**2
    Lp = np.zeros((2 * n, 2 * n))
    Lp[:n, :n] = lap + np.diag(S * S + 3 * a * F * F - Om * Om)
    Lp[n:, n:] = lap + np.diag(D)
    Lp[:n, n:] = np.diag(2 * F * S)
    Lp[n:, :n] = np.diag(2 * F * S)
    vec = np.r_[u, np.zeros(n)]
    K = 4 * Om * Om * np.outer(vec, vec) / np.dot(u, u)          # h cancels between <u,eta> and <u,u>
    Lm = lap + np.diag(S * S + a * F * F - Om * Om)
    ev_Lp = eigh(Lp, eigvals_only=True, subset_by_index=[0, 4])
    ev_LpK = eigh(Lp + K, eigvals_only=True, subset_by_index=[0, 4])
    ev_Lm = eigh(Lm, eigvals_only=True, subset_by_index=[0, 2])
    ell1 = Lp + np.diag(np.tile(2 / r**2, 2))
    ev_l1 = eigh(ell1, eigvals_only=True, subset_by_index=[0, 2])
    # eigenvector of the lowest L_+ + K mode: is it localised on the object or spread over the box?
    evals, evecs = eigh(Lp + K, subset_by_index=[0, 0])
    v0 = evecs[:, 0]
    frac_inside = float(np.sum(v0[:n][r < 15]**2 + v0[n:][r < 15]**2) / np.sum(v0**2))
    # the negative L_+ eigenvector: overlap with the charge direction u
    evn, evcn = eigh(Lp, subset_by_index=[0, 0])
    vneg = evcn[:, 0]
    overlap = float(np.dot(vneg[:n], u) / math.sqrt(np.dot(u, u)))
    row = dict(L=L, h=h, Omega=Om, continuum_edge_1_minus_Om2=1 - Om * Om, box_estimate=1 - Om * Om + (math.pi / L)**2,
               Lplus_lowest=ev_Lp.tolist(), Lplus_plus_K_lowest=ev_LpK.tolist(), Lminus_lowest=ev_Lm.tolist(),
               ell1_lowest=ev_l1.tolist(), lowest_LpK_mode_fraction_inside_rho15=frac_inside,
               Lplus_negative_mode_overlap_with_u=overlap,
               lowest_LpK_minus_edge_times_L2_over_pi2=(ev_LpK[0] - (1 - Om * Om)) * L * L / math.pi**2)
    out['boxes'].append(row)
    print('L=%5.0f h=%.2f  L+ min %.5f | L+ +K lowest %s | edge 1-Om^2 %.5f box est %.5f | L- min %.2e | ell1 min %.2e | inside15 %.3f | (lam-edge)L^2/pi^2 %.3f' % (
        L, h, ev_Lp[0], np.array2string(ev_LpK[:3], precision=5), 1 - Om * Om, 1 - Om * Om + (math.pi / L)**2, ev_Lm[0], ev_l1[0], frac_inside, row['lowest_LpK_minus_edge_times_L2_over_pi2']), flush=True)

out['codex'] = dict(radial_fixed_charge_lowest={'L80': 0.2243396857099204, 'L120': 0.2233198384209463},
                    claim='Four radial grids give positive lowest fixed-charge energy eigenvalues (.2233-.2244 in m_X^2 units)')
json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)
