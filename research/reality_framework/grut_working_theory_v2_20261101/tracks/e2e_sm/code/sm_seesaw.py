#!/usr/bin/env python3
"""e2e_sm: type-I seesaw completion of the HFa/HFb axial singlets (model choices throughout).

Light spectrum: normal ordering, m1 = 0 (choice), NuFIT 6.0 IC24 (with SK-atm) best fit [Q, arXiv:2410.05380 Table 1]:
  dm21^2 = 7.49e-5 eV^2, dm31^2 = 2.513e-3 eV^2.
Heavy spectrum: the four axial singlets of mf_hybrid, M_j/H_i = 250 (N1, sigma), 400, 200, 120 (N2-4, chi).
Casas-Ibarra (hep-ph/0103065): m_D = i U* sqrt(m) R sqrt(M), R (3x4) complex with R R^T = 1_3, built as the first
three rows of a 4x4 complex orthogonal O = prod of six complex plane rotations (12 real parameters).
In the unflavoured regime h^dag h = sqrt(M) R^dag m R sqrt(M) / v^2 does not depend on U [D].

Usage: python3 -B sm_seesaw.py [H_i_GeV]   -> results/seesaw.json
"""
import json, sys, math
import numpy as np
from scipy.optimize import least_squares

V_EW = 174.1                      # GeV (as mf_hybrid / e2e_kinetic)
DM21, DM31 = 7.49e-5, 2.513e-3    # eV^2, NuFIT 6.0 IC24 NO best fit [Q]
M_LIGHT = np.array([0.0, math.sqrt(DM21), math.sqrt(DM31)])   # eV, m1 = 0 (choice)
MSTAR = 1.08e-3                   # eV, equilibrium neutrino mass (Buchmuller-Di Bari-Plumacher hep-ph/0401240) [Q]
M_OVER_H = np.array([250.0, 400.0, 200.0, 120.0])             # N1..N4 in H_i (mf_hybrid bench)
PARENT = ['sigma', 'chi', 'chi', 'chi']
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def ortho4(p):
    """4x4 complex orthogonal matrix from 6 complex angles p = [re..., im...]."""
    z = p[:6] + 1j * p[6:]
    O = np.eye(4, dtype=complex)
    for (a, b), zz in zip(PAIRS, z):
        G = np.eye(4, dtype=complex)
        c, s = np.cos(zz), np.sin(zz)
        G[a, a] = c; G[b, b] = c; G[a, b] = s; G[b, a] = -s
        O = O @ G
    return O


def fV(x):  # vertex
    return math.sqrt(x) * (1.0 - (1.0 + x) * math.log((1.0 + x) / x))


def fS(x):  # self-energy (non-degenerate)
    return math.sqrt(x) / (1.0 - x)


def derived(R, M_GeV):
    """m-tilde_j, h^dag h, Gamma_j (GeV), K_j, eps_j for a CI matrix R (3x4)."""
    m_eV = M_LIGHT
    A = R.conj().T @ np.diag(m_eV) @ R           # 4x4, eV: (R^dag m R)
    mt = np.real(np.diag(A))                      # m-tilde_j (eV)
    sq = np.sqrt(M_GeV)
    hh = (sq[:, None] * A * sq[None, :]) * 1e-9 / V_EW**2      # h^dag h (dimensionless)
    Gam = np.real(np.diag(hh)) * M_GeV / (8 * math.pi)         # rest width, GeV
    K = mt / MSTAR
    eps = np.zeros(4)
    for i in range(4):
        for j in range(4):
            if j == i:
                continue
            x = (M_GeV[j] / M_GeV[i])**2
            eps[i] += -(1 / (8 * math.pi)) * np.imag(hh[i, j]**2) / np.real(hh[i, i]) * (fV(x) + fS(x))
    return mt, hh, Gam, K, eps


def mnu_tree(R, M_GeV):
    """Rebuild the light mass matrix from m_D (U = 1 suffices for eigenvalues) and return its singular values."""
    mD = 1j * np.diag(np.sqrt(M_LIGHT * 1e-9)) @ R @ np.diag(np.sqrt(M_GeV))     # GeV
    mnu = -mD @ np.diag(1.0 / M_GeV) @ mD.T
    return np.sort(np.linalg.svd(mnu, compute_uv=False)) * 1e9, mD


def mnu_loop(mD, M_GeV):
    """One-loop (Z, h) correction, Grimus-Lavoura-type [E]: dm = g^2/(64 pi^2 mW^2) mD diag(f(M)/M) mD^T with
    f = 3 mZ^2 ln(M^2/mZ^2) + mh^2 ln(M^2/mh^2) (heavy-M limit). Returns singular values of tree+loop and of loop
    with the universal log removed (the part that does not track the tree cancellation)."""
    g2, mW, mZ, mh = 0.4225, 80.37, 91.19, 125.2
    f = 3 * mZ**2 * np.log(M_GeV**2 / mZ**2) + mh**2 * np.log(M_GeV**2 / mh**2)
    c = g2 / (64 * math.pi**2 * mW**2)
    dm = c * mD @ np.diag(f / M_GeV) @ mD.T
    fbar = np.mean(f)
    dm_nonuni = c * mD @ np.diag((f - fbar) / M_GeV) @ mD.T
    tree = -mD @ np.diag(1.0 / M_GeV) @ mD.T
    sv = lambda X: np.sort(np.linalg.svd(X, compute_uv=False)) * 1e9
    return dict(universal_fraction=float(c * fbar), sv_tree_plus_loop_eV=sv(tree + dm).tolist(),
                sv_nonuniversal_loop_eV=sv(dm_nonuni).tolist())


def solve(mt_target, eps4_target, M_GeV, seed=0, reg=1e-3, eps_scale=None):
    rng = np.random.default_rng(seed)
    p0 = np.concatenate([rng.uniform(-1.5, 1.5, 6), rng.uniform(-1.0, 1.0, 6)])
    es = eps_scale if eps_scale else 1e-6

    def res(p):
        R = ortho4(p)[:3, :]
        mt, hh, G, K, eps = derived(R, M_GeV)
        r = list(mt / mt_target - 1.0)
        if eps4_target is not None:
            r.append((eps[3] - eps4_target) / es)
        r += list(reg * p)
        return np.array(r)
    sol = least_squares(res, p0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=400)
    return sol.x, sol


def card(p, M_GeV, H):
    R = ortho4(p)[:3, :]
    mt, hh, G, K, eps = derived(R, M_GeV)
    sv, mD = mnu_tree(R, M_GeV)
    loop = mnu_loop(mD, M_GeV)
    return dict(R_re=np.real(R).tolist(), R_im=np.imag(R).tolist(), z_re=p[:6].tolist(), z_im=p[6:].tolist(),
                max_abs_R=float(np.max(np.abs(R[1:, :]))),
                sum_abs_R2_rows23=[float(np.sum(np.abs(R[1]) ** 2)), float(np.sum(np.abs(R[2]) ** 2))],
                mtilde_eV=mt.tolist(), y2_hdagh_jj=np.real(np.diag(hh)).tolist(),
                Gamma_rest_GeV=G.tolist(), Gamma_rest_over_Hi=(G / H).tolist(), K_washout=K.tolist(),
                eps=eps.tolist(), light_masses_rebuilt_eV=sv.tolist(), sum_mnu_eV=float(np.sum(sv)),
                RRT_minus_1=float(np.max(np.abs(R @ R.T - np.eye(3)))), one_loop=loop)


ETA = np.array([1.0, -1.0, 1.0, 1.0])   # CP parities of N1..N4 at delta = 0 (model choice: N2 odd)
CPH = np.where(ETA > 0, 1.0 + 0j, 1j)  # column phases c_j: R0_ij = r_ij c_j
# CP-conserving R (columns real or imaginary) requires sum_j eta_j mt_j = sum_i m_i  [D]
A_ODD = np.array([[0, 0.6, -0.3, 0.8], [-0.6, 0, 0.5, -0.7], [0.3, -0.5, 0, 0.4], [-0.8, 0.7, -0.4, 0]])  # fixed


def cp_even_base(mt_target, seeds=12):
    """Real rows u (m2), w (m3), r0 (m1 = 0) with eta-metric orthonormality and m2 u_j^2 + m3 w_j^2 = mt_j."""
    m2, m3 = M_LIGHT[1], M_LIGHT[2]
    eta = ETA

    def res(x):
        u, w = x[:4], x[4:]
        r = [np.sum(eta * u * u) - 1, np.sum(eta * w * w) - 1, np.sum(eta * u * w)]
        r += list((m2 * u**2 + m3 * w**2) / mt_target - 1)
        return np.array(r + list(1e-6 * x))
    best = None
    for sd in range(seeds):
        x0 = np.random.default_rng(100 + sd).normal(0, 1.5, 8)
        sol = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
        err = np.max(np.abs(res(sol.x)[:7]))
        if err < 1e-8 and (best is None or np.sum(sol.x**2) < np.sum(best**2)):
            best = sol.x
    if best is None:
        raise RuntimeError('no CP-even base')
    u, w = best[:4], best[4:]

    def res0(r0):
        return np.array([np.sum(eta * r0 * r0) - 1, np.sum(eta * r0 * u), np.sum(eta * r0 * w)] + list(1e-6 * r0))
    s0 = least_squares(res0, np.array([0.3, 0.1, 0.9, 0.2]), xtol=1e-15, ftol=1e-15, gtol=1e-15)
    r = np.vstack([s0.x, u, w])
    return r * CPH[None, :]


def expm_small(A):
    from scipy.linalg import expm
    return expm(A)


def family(mt_target, M_GeV, deltas):
    """R(delta) = R0 exp(delta A), A complex antisymmetric (so R R^T = 1 exactly) and CP-odd: P A* P = -A with
    P = diag(eta), i.e. imaginary entries between same-parity columns, real entries between opposite-parity ones.
    eps_j linear in delta at small delta; m-tilde shifts are O(delta^2)."""
    R0 = cp_even_base(mt_target)
    same = (ETA[:, None] * ETA[None, :]) > 0
    A = np.where(same, 1j, 1.0) * A_ODD
    return R0, [(d, R0 @ expm_small(d * A)) for d in deltas]


def card_R(R, M_GeV, H):
    mt, hh, G, K, eps = derived(R, M_GeV)
    sv, mD = mnu_tree(R, M_GeV)
    loop = mnu_loop(mD, M_GeV)
    return dict(R_re=np.real(R).tolist(), R_im=np.imag(R).tolist(),
                max_abs_R=float(np.max(np.abs(R[1:, :]))),
                mtilde_eV=mt.tolist(), y2_hdagh_jj=np.real(np.diag(hh)).tolist(),
                Gamma_rest_GeV=G.tolist(), Gamma_rest_over_Hi=(G / H).tolist(), K_washout=K.tolist(),
                eps=eps.tolist(), light_masses_rebuilt_eV=sv.tolist(), sum_mnu_eV=float(np.sum(sv)),
                RRT_minus_1=float(np.max(np.abs(R @ R.T - np.eye(3)))), one_loop=loop)


def main():
    H = float(sys.argv[1]) if len(sys.argv) > 1 else 1.3e11
    M_GeV = M_OVER_H * H
    out = dict(label='C on model choices (e2e_sm); NO, m1 = 0; unflavoured CP asymmetries (Covi-Roulet-Vissani)',
               H_i_GeV=H, M_GeV=M_GeV.tolist(), M_over_Hi=M_OVER_H.tolist(), parent=PARENT,
               light_masses_eV=M_LIGHT.tolist(), sum_mnu_eV=float(M_LIGHT.sum()), v_EW=V_EW,
               nufit='NuFIT 6.0 IC24 NO: dm21^2 = 7.49e-5, dm31^2 = 2.513e-3 eV^2 (arXiv:2410.05380)',
               family_note='R(delta) = R0 exp(delta A): R0 CP-conserving (columns real, N2 column imaginary: CP parities '
                           'eta = (+,-,+,+)), least-norm solution of m-tilde_j = target; A_jk = a_jk (x i for same-parity pairs) with a '
                           'fixed real antisymmetric a (A_ODD), times i between same-parity columns: the CP-odd direction. delta = 0: eps_j = 0 exactly. '
                           'delta is THE explicit CP phase of the card (supplied); m-tilde shifts are O(delta^2).',
               A_ODD=A_ODD.tolist(), ETA=ETA.tolist())
    out['eps_DI_lightest'] = float(3 / (16 * math.pi) * M_GeV[3] * (M_LIGHT[2] - M_LIGHT[0]) * 1e-9 / V_EW**2)
    deltas = [0.0, 1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0]
    scen = {}
    S = M_LIGHT.sum()
    # eta = (+, -, +, +): S1 keeps every m-tilde >= 0.05 eV (mf_hybrid's rates and N-gas medium mass, or faster)
    for name, mtv in [('S1_mt050', np.array([0.05, 0.15 - S, 0.05, 0.05])),
                      ('S2_mt030', np.array([0.03, 0.03, 0.03, S - 0.03]))]:
        mt_target = mtv
        R0, fam = family(mt_target, M_GeV, deltas)
        rows = []
        for d, R in fam:
            c = card_R(R, M_GeV, H)
            err = float(np.max(np.abs(np.array(c['mtilde_eV']) / mt_target - 1)))
            c.update(delta=d, mt_err=err)
            rows.append(c)
            print(name, 'delta %.0e' % d, 'err %.1e' % err, 'eps', ['%.2e' % e for e in c['eps']],
                  'sum %.4f' % c['sum_mnu_eV'], 'maxR %.2f' % c['max_abs_R'],
                  'loopNU', ['%.1e' % x for x in c['one_loop']['sv_nonuniversal_loop_eV']])
        scen[name] = dict(mtilde_target_eV=mtv.tolist(), tuning_sum_mt_over_sum_m=float(np.sum(mtv) / M_LIGHT.sum()),
                          rows=rows)
    out['scenarios'] = scen
    path = sys.argv[2] if len(sys.argv) > 2 else '../results/seesaw.json'
    json.dump(out, open(path, 'w'), indent=1)
    print('eps_DI_lightest', out['eps_DI_lightest'])


if __name__ == '__main__':
    main()
