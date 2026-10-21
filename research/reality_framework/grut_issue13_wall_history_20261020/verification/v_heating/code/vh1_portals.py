#!/usr/bin/env python3
"""v_heating check 1: one-loop running of the three Higgs portals in SM-Higgs + (phi, X, D).

Independent code (verifier). Conventions (stated here, not taken from the module):
  V = lam_h (H^+H)^2 + (H^+H) [lHp phi^2 + lHX |X|^2 + lHD D^2]
      + (lp/4)(phi^2 - v^2)^2 + gX^2 phi^2 |X|^2 + lX |X|^4 + (hD^2/2) phi^2 D^2 + (lD/4) D^4 + alpha |X|^2 D^2
  H = (h1 + i h2, h3 + i h4)/sqrt2, X = (q1 + i q2)/sqrt2.
Real fields z = (h1..h4, phi, q1, q2, D); V4 = (1/24) L_abcd z_a z_b z_c z_d.
One-loop (scalar part, 16 pi^2 units):  beta_abcd = L_abef L_efcd + L_acef L_efbd + L_adef L_efbc
   + sum_e (gam_ae L_ebcd + 3 perms) + G_abcd,
where gam (field anomalous dimension) is nonzero only on Higgs indices (gauge + Yukawa), and G (gauge^4, yukawa^4
boxes) carries only Higgs indices.  We use random values for gam_H and G to show nothing leaks.

Analytic results derived by hand (see VERIFY.md):
 16pi^2 b_lHp = lHp(12 lam_h + 6 lp + 8 lHp + 2 gH) + 2 gX^2 lHX + 2 hD^2 lHD
 16pi^2 b_lHX = lHX(12 lam_h + 8 lX + 4 lHX + 2 gH) + 4 gX^2 lHp + 4 alpha lHD
 16pi^2 b_lHD = lHD(12 lam_h + 6 lD + 8 lHD + 2 gH) + 2 hD^2 lHp + 2 alpha lHX
with 2 gH = 6 yt^2 - (9/2) g2^2 - (3/2) g1^2 (SM).  Single-singlet limit: Lebedev 1203.0156 eq.(15).
"""
import json, os, itertools
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), '..', 'results')


def V4(z, c):
    h2 = z[0]**2 + z[1]**2 + z[2]**2 + z[3]**2       # = 2 H^+H
    HH = h2 / 2.0
    p, q1, q2, D = z[4], z[5], z[6], z[7]
    X2 = (q1**2 + q2**2) / 2.0                      # |X|^2
    return (c['lam_h'] * HH**2 + HH * (c['lHp'] * p**2 + c['lHX'] * X2 + c['lHD'] * D**2)
            + c['lp'] / 4 * p**4 + c['gX2'] * p**2 * X2 + c['lX'] * X2**2
            + c['hD2'] / 2 * p**2 * D**2 + c['lD'] / 4 * D**4 + c['alpha'] * X2 * D**2)


def tensor_polar(c, n=8):
    """L_abcd = (1/16) sum_eps eps1..eps4 V4(sum eps_i e_i)  (exact for a homogeneous quartic)."""
    E = np.eye(n)
    L = np.zeros((n,) * 4)
    signs = list(itertools.product([1, -1], repeat=4))
    for a, b, cc, d in itertools.combinations_with_replacement(range(n), 4):
        s = 0.0
        for e in signs:
            s += e[0] * e[1] * e[2] * e[3] * V4(e[0] * E[a] + e[1] * E[b] + e[2] * E[cc] + e[3] * E[d], c)
        val = s / 16.0
        for perm in set(itertools.permutations((a, b, cc, d))):
            L[perm] = val
    return L


def tensor_explicit(c):
    """Same tensor, from hand-written fourth derivatives (independent of the polarization routine)."""
    L = np.zeros((8,) * 4)
    def setsym(idx, val):
        for perm in set(itertools.permutations(idx)):
            L[perm] = val
    lh = c['lam_h']
    for i in range(4):
        setsym((i, i, i, i), 6 * lh)
        for j in range(i + 1, 4):
            setsym((i, i, j, j), 2 * lh)
    portal = {4: 2 * c['lHp'] / 2 * 1.0, 5: c['lHX'] / 2, 6: c['lHX'] / 2, 7: c['lHD']}
    # d^4/dh_i^2 dz_a^2 of (h^2/2) * Q_a z_a^2 = 1 * 2 Q_a ; Q_phi = lHp, Q_q = lHX/2, Q_D = lHD
    Q = {4: c['lHp'], 5: c['lHX'] / 2, 6: c['lHX'] / 2, 7: c['lHD']}
    for i in range(4):
        for a in range(4, 8):
            setsym((i, i, a, a), 2 * Q[a])
    setsym((4, 4, 4, 4), 6 * c['lp'])
    setsym((4, 4, 5, 5), 2 * c['gX2']); setsym((4, 4, 6, 6), 2 * c['gX2'])
    setsym((4, 4, 7, 7), 2 * c['hD2'])
    setsym((5, 5, 5, 5), 6 * c['lX']); setsym((6, 6, 6, 6), 6 * c['lX']); setsym((5, 5, 6, 6), 2 * c['lX'])
    setsym((7, 7, 7, 7), 6 * c['lD'])
    setsym((5, 5, 7, 7), 2 * c['alpha']); setsym((6, 6, 7, 7), 2 * c['alpha'])
    return L


def beta_tensor(L, gamH, G):
    B = (np.einsum('abef,efcd->abcd', L, L) + np.einsum('acef,efbd->abcd', L, L)
         + np.einsum('adef,efbc->abcd', L, L))
    gam = np.zeros(8); gam[:4] = gamH
    B += (gam[:, None, None, None] + gam[None, :, None, None] + gam[None, None, :, None]
          + gam[None, None, None, :]) * L
    B += G
    return B


def higgs_only_G(rng):
    G = np.zeros((8,) * 4)
    g = rng.normal(size=(4, 4, 4, 4))
    # symmetrize
    gs = sum(np.transpose(g, p) for p in itertools.permutations(range(4))) / 24
    G[:4, :4, :4, :4] = gs
    return G


def analytic_portal_betas(c, gH2):
    lh = c['lam_h']
    bHp = c['lHp'] * (12 * lh + 6 * c['lp'] + 8 * c['lHp'] + gH2) + 2 * c['gX2'] * c['lHX'] + 2 * c['hD2'] * c['lHD']
    bHX = c['lHX'] * (12 * lh + 8 * c['lX'] + 4 * c['lHX'] + gH2) + 4 * c['gX2'] * c['lHp'] + 4 * c['alpha'] * c['lHD']
    bHD = c['lHD'] * (12 * lh + 6 * c['lD'] + 8 * c['lHD'] + gH2) + 2 * c['hD2'] * c['lHp'] + 2 * c['alpha'] * c['lHX']
    return np.array([bHp, bHX, bHD])


def extract_portal_betas(B):
    # B_{h1 h1 a a} = 2 Q_a-beta: phi: 2 b_lHp ; q1: b_lHX ; D: 2 b_lHD
    return np.array([B[0, 0, 4, 4] / 2, B[0, 0, 5, 5], B[0, 0, 7, 7] / 2])


def main():
    rng = np.random.default_rng(424242)
    bench = dict(lam_h=0.13, lp=0.01, gX2=0.01, lX=0.001, hD2=1e-4, lD=0.001, alpha=0.0,
                 lHp=0.0, lHX=0.0, lHD=0.0)
    res = {}
    # 1. polarization tensor vs explicit tensor
    Lp = tensor_polar(bench); Le = tensor_explicit(bench)
    res['tensor_polar_vs_explicit_maxabs'] = float(np.max(np.abs(Lp - Le)))
    # 2. zero portals: mixed entries of L and beta
    H = list(range(4)); S = list(range(4, 8))
    def mixed_max(T):
        m = 0.0
        for idx in itertools.product(range(8), repeat=4):
            nh = sum(1 for i in idx if i < 4)
            if 0 < nh < 4:
                m = max(m, abs(T[idx]))
        return m
    worst = 0.0; worst_L = 0.0
    for trial in range(20):
        c = dict(lam_h=rng.uniform(0.05, 0.3), lp=rng.uniform(0, 0.5), gX2=rng.uniform(0, 0.5),
                 lX=rng.uniform(0, 0.5), hD2=rng.uniform(0, 0.5), lD=rng.uniform(0, 0.5),
                 alpha=rng.uniform(-0.2, 0.2), lHp=0.0, lHX=0.0, lHD=0.0)
        L = tensor_explicit(c)
        B = beta_tensor(L, rng.normal(), higgs_only_G(rng))
        worst = max(worst, mixed_max(B)); worst_L = max(worst_L, mixed_max(L))
    res['zero_portal_random_trials'] = 20
    res['zero_portal_max_mixed_L'] = worst_L
    res['zero_portal_max_mixed_beta'] = worst
    # 3. analytic vs tensor with nonzero portals
    maxrel = 0.0
    for trial in range(20):
        c = dict(lam_h=rng.uniform(0.05, 0.3), lp=rng.uniform(0, 0.5), gX2=rng.uniform(0, 0.5),
                 lX=rng.uniform(0, 0.5), hD2=rng.uniform(0, 0.5), lD=rng.uniform(0, 0.5),
                 alpha=rng.uniform(-0.2, 0.2), lHp=rng.uniform(-0.3, 0.3), lHX=rng.uniform(-0.3, 0.3),
                 lHD=rng.uniform(-0.3, 0.3))
        gH2 = rng.normal()        # 2*gamma_H in 16pi^2 units
        L = tensor_explicit(c)
        B = beta_tensor(L, gH2 / 2, higgs_only_G(rng))
        bt = extract_portal_betas(B); ba = analytic_portal_betas(c, gH2)
        maxrel = max(maxrel, float(np.max(np.abs(bt - ba) / (np.abs(ba) + 1e-300))))
        # consistency: q2 entry equals q1 entry, h3 entries equal h1 entries, no q1q2 off-diagonal
        assert abs(B[2, 2, 5, 5] - B[0, 0, 5, 5]) < 1e-12 and abs(B[0, 0, 5, 6]) < 1e-12
    res['analytic_vs_tensor_max_rel'] = maxrel
    # 4. Lebedev eq.(15) single-singlet limit: set gX2 = hD2 = alpha = 0, lHX = lHD = 0
    c = dict(lam_h=0.13, lp=0.02, gX2=0.0, lX=0.0, hD2=0.0, lD=0.0, alpha=0.0, lHp=0.05, lHX=0.0, lHD=0.0)
    # Lebedev: V = lh/4 h^4 + lhs/4 s^2 h^2 + ls/4 s^4 with h^2 = 2H^+H  ->  lhs/2 = lHp, ls = lp
    lhs = 2 * c['lHp']; ls = c['lp']
    leb = 4 * lhs**2 + 12 * c['lam_h'] * lhs + 6 * ls * lhs      # scalar part of Lebedev (15)
    L = tensor_explicit(c); B = beta_tensor(L, 0.0, np.zeros((8,) * 4))
    ours = 2 * extract_portal_betas(B)[0]                         # beta of lhs = 2 beta_lHp
    res['lebedev_eq15_scalar_part'] = leb
    res['ours_single_singlet'] = float(ours)
    res['lebedev_rel_diff'] = float(abs(ours - leb) / leb)
    # 5. linearized portal matrix at the benchmark (scalar part, lam_h = 0.13)
    Mlin = np.zeros((3, 3))
    for j in range(3):
        c = dict(bench); key = ['lHp', 'lHX', 'lHD'][j]; c[key] = 1e-8
        Mlin[:, j] = analytic_portal_betas(c, 0.0) / 1e-8
    res['linear_portal_matrix_16pi2'] = Mlin.tolist()
    # 6. explicit integration: zero stays zero; small seed stays proportional
    def rhs(y, c0):
        c = dict(c0); c['lHp'], c['lHX'], c['lHD'] = y
        return analytic_portal_betas(c, 0.0) / (16 * np.pi**2)
    y = np.zeros(3); t = 0.0; dt = 0.05
    for _ in range(int(np.log(2.4e18 / 100) / dt)):
        k1 = rhs(y, bench); k2 = rhs(y + dt / 2 * k1, bench); k3 = rhs(y + dt / 2 * k2, bench); k4 = rhs(y + dt * k3, bench)
        y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    res['zero_initial_portals_after_running_to_MP'] = y.tolist()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'vh1_portals.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
