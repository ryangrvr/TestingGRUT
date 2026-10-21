"""Part 2 (X fluctuations): 2D Fokker-Planck in (phi, rho), rho^2 = 2|X|^2 = q1^2 + q2^2, D = 0.

X is light wherever g_X |phi| < H, so it fluctuates during inflation. With the O(2) phase integrated out
the stochastic process lives on (y, r) = (phi, rho)/H with measure r dr dy and
    P_eq ~ r exp(-s(y,r)),  s = (8 pi^2/3) [ w_phi(y) + c y^2 r^2/2 + lx r^4/4 + b2 r^2/2 ].
SqRA (detailed-balance-exact) discretisation on a cylindrical cell grid; phi-odd sector = half grid y>0
with killing 2D/hy^2 on the cells next to y = 0 (same identity as fp1d).

Outputs:
  (i) Lambda_1(2D) vs Lambda_1(1D, X=0) and vs the entropic Langer estimate Lambda_1(1D) R_X;
  (ii) H >> v: marginal P_eq(phi) with X integrated (shape, <phi^2>), p_- for the fresh start
       (phi = +v, X = 0) and the distribution of the minority fraction for equilibrium draws (phi*, X*).
"""
import time
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, bmat
from scipy.sparse import identity
from scipy.sparse.linalg import eigsh, splu
from common import Spectator, V, DIFF, BETA_FP, PC_SITE_SC, dump
import fp1d
from stochastic_scan import entropic_factors


def grid2d(sp, hy, hr, smax=60.0):
    ym = sp.ymin()
    smin = sp.s(np.array([ym]))[0]
    y = ym + 1.0
    while sp.s(np.array([y]))[0] - smin < smax:
        y *= 1.1
    r = 1.0
    while sp.s(np.array([0.0]), np.array([r]))[0] - sp.s(np.array([0.0]))[0] < smax:
        r *= 1.1
    ys = (np.arange(int(np.ceil(y / hy))) + 0.5) * hy
    rs = (np.arange(int(np.ceil(r / hr))) + 0.5) * hr
    return ys, rs


def build(sp, ys, rs, hy, hr):
    Ny, Nr = len(ys), len(rs)
    Y, R = np.meshgrid(ys, rs, indexing="ij")
    s = sp.s(Y, R)
    idx = np.arange(Ny * Nr).reshape(Ny, Nr)
    rows, cols, vals = [], [], []        # generator A (row = from), off-diagonal rates
    out = np.zeros((Ny, Nr))
    srows, scols, svals = [], [], []     # symmetrised off-diagonals
    # phi direction
    kf = DIFF / hy ** 2 * np.exp(-(s[1:, :] - s[:-1, :]) / 2)
    kb = DIFF / hy ** 2 * np.exp(-(s[:-1, :] - s[1:, :]) / 2)
    a, b = idx[:-1, :].ravel(), idx[1:, :].ravel()
    rows += [a, b]; cols += [b, a]; vals += [kf.ravel(), kb.ravel()]
    out[:-1, :] += kf; out[1:, :] += kb
    sym = -DIFF / hy ** 2 * np.ones(a.size)
    srows += [a, b]; scols += [b, a]; svals += [sym, sym]
    # rho direction
    rf = (np.arange(Nr - 1) + 1.0) * hr        # r_{j+1/2}
    kup = DIFF / hr ** 2 * (rf / rs[:-1])[None, :] * np.exp(-(s[:, 1:] - s[:, :-1]) / 2)
    kdn = DIFF / hr ** 2 * (rf / rs[1:])[None, :] * np.exp(-(s[:, :-1] - s[:, 1:]) / 2)
    a, b = idx[:, :-1].ravel(), idx[:, 1:].ravel()
    rows += [a, b]; cols += [b, a]; vals += [kup.ravel(), kdn.ravel()]
    out[:, :-1] += kup; out[:, 1:] += kdn
    sym = (-DIFF / hr ** 2 * rf[None, :] / np.sqrt(rs[:-1] * rs[1:])[None, :] * np.ones((Ny, 1))).ravel()
    srows += [a, b]; scols += [b, a]; svals += [sym, sym]
    kill = np.zeros((Ny, Nr))
    kill[0, :] = 2 * DIFF / hy ** 2
    n = Ny * Nr
    diag = (out + kill).ravel()
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    A = A - csr_matrix((diag, (np.arange(n), np.arange(n))), shape=(n, n))
    M = coo_matrix((np.concatenate(svals), (np.concatenate(srows), np.concatenate(scols))), shape=(n, n)).tocsr()
    M = M + csr_matrix((diag, (np.arange(n), np.arange(n))), shape=(n, n))
    logpi = (np.log(R) - s).ravel()
    return A, M, s, logpi, kill.ravel(), idx


def lambda1_2d(M):
    vals = eigsh(M, k=2, sigma=0.0, which="LM", return_eigenvectors=False)
    return float(np.min(vals))


def bdf2(A, v0, Ns, dt=0.05):
    """L-stable BDF2 propagation of dv/dN = A v (first step backward Euler); returns {N: v(N)}.
    (scipy's expm_multiply was ~100x slower here; dt convergence is checked in main().)"""
    n = A.shape[0]
    I = identity(n, format="csc")
    lu1 = splu((I - dt * A).tocsc())
    lu2 = splu((I - (2.0 / 3.0) * dt * A).tocsc())
    steps = {int(round(N / dt)): N for N in Ns}
    vprev, v = v0, lu1.solve(v0)
    out = {}
    if 1 in steps:
        out[steps[1]] = v.copy()
    for k in range(2, max(steps) + 1):
        vprev, v = v, lu2.solve((4.0 / 3.0) * v - (1.0 / 3.0) * vprev)
        if k in steps:
            out[steps[k]] = v.copy()
    return out


def survival(A, Ns, dt=0.05):
    res = bdf2(A, np.ones(A.shape[0]), Ns, dt)
    return {N: np.clip(v, 0.0, 1.0) for N, v in res.items()}


def main():
    t0 = time.time()
    out = {"module": "inflation_spectator", "part": "2 X fluctuations (2D phi-|X| Fokker-Planck)",
           "label": "computed; D set to zero (its entropic factor R_D is computed in stochastic_1d.json)",
           "rate_checks": [], "high_H": []}
    # (i) rate checks at moderate barriers
    for H in (1000.0, 500.0, 400.0, 330.0, 280.0, 250.0):
        sp = Spectator(H)
        vH = V / H
        sw = 1.0 / np.sqrt(BETA_FP * 2 * sp.lp * max(sp.ymin(), 0.5) ** 2)
        hy = min(0.04, sw / 10)
        sr = 1.0 / np.sqrt(BETA_FP * sp.c * max(sp.ymin(), 0.5) ** 2)
        hr = min(0.05, sr / 8)
        ys, rs = grid2d(sp, hy, hr)
        A, M, s, logpi, kill, idx = build(sp, ys, rs, hy, hr)
        L2 = lambda1_2d(M)
        yy, h = fp1d.make_grid(sp)
        L1, _ = fp1d.lambda1_mfpt(sp, yy, h)
        RX, RD = entropic_factors(sp)
        row = {"H_GeV": H, "B_HM": sp.barrier_B(), "grid": [len(ys), len(rs)], "hy": hy, "hr": hr,
               "Lambda1_2D": L2, "Lambda1_1D_X0": L1, "ratio_2D_over_1D": L2 / L1, "R_X_entropic_estimate": RX,
               "ratio_over_RX": L2 / (L1 * RX)}
        out["rate_checks"].append(row)
        print(row, flush=True)
    # grid refinement of one rate check
    H = 330.0
    sp = Spectator(H)
    sw = 1.0 / np.sqrt(BETA_FP * 2 * sp.lp * sp.ymin() ** 2)
    sr = 1.0 / np.sqrt(BETA_FP * sp.c * sp.ymin() ** 2)
    res = []
    for f in (1.0, 0.7):
        hy, hr = min(0.04, sw / 10) * f, min(0.05, sr / 8) * f
        ys, rs = grid2d(sp, hy, hr)
        A, M, s, logpi, kill, idx = build(sp, ys, rs, hy, hr)
        res.append(lambda1_2d(M))
    out["rate_check_refinement_H330"] = {"Lambda1_2D_coarse_fine": res, "rel_change": res[1] / res[0] - 1}
    # fresh 2D start (phi=+v, X=0) at H = 330 GeV, N = 50, 60
    hy, hr = min(0.04, sw / 10), min(0.05, sr / 8)
    ys, rs = grid2d(sp, hy, hr)
    A, M, s, logpi, kill, idx = build(sp, ys, rs, hy, hr)
    i0 = int(np.argmin(np.abs(ys - V / H)))
    m0 = idx[i0, 0]
    centre = idx[0, :]
    # forward absorbed flux sum_c kill_c p_c(t), with p_c = (pi_c/pi_m0) [e^{At} e_m0]_c (reversibility)
    w = kill[centre] * np.exp(logpi[centre] - logpi[m0])
    n = A.shape[0]
    e = csr_matrix((w, (np.zeros(centre.size, dtype=int), centre)), shape=(1, n))
    Aug = bmat([[A, None], [e, csr_matrix((1, 1))]], format="csr")
    v0 = np.zeros(n + 1)
    v0[m0] = 1.0
    res = bdf2(Aug, v0, (50.0, 60.0), dt=0.05)
    res2 = bdf2(Aug, v0, (50.0,), dt=0.025)
    I50, I60, I50b = res[50.0][-1], res[60.0][-1], res2[50.0][-1]
    yy, h = fp1d.make_grid(sp)
    m01 = int(np.argmin(np.abs(yy - V / H)))
    p1 = fp1d.p_minus_delta(sp, yy, h, m01, [50.0, 60.0])
    out["fresh_start_2D_H330"] = {"p_minus_2D_N50": 0.5 * I50, "p_minus_2D_N60": 0.5 * I60,
                                  "p_minus_1D_N50": p1[50.0]["p_minus"], "p_minus_1D_N60": p1[60.0]["p_minus"],
                                  "p_minus_2D_N50_dt_half": 0.5 * I50b}
    print(out["fresh_start_2D_H330"], flush=True)
    # (ii) H >> v
    for H in (1e4, 1e13):
        sp = Spectator(H)
        hy, hr = 0.05, 0.05
        ys, rs = grid2d(sp, hy, hr)
        A, M, s, logpi, kill, idx = build(sp, ys, rs, hy, hr)
        L2 = lambda1_2d(M)
        pi = np.exp(logpi - logpi.max())
        pi /= pi.sum()
        Ny, Nr = len(ys), len(rs)
        Pphi = pi.reshape(Ny, Nr).sum(axis=1)
        Pphi /= Pphi.sum()
        phi2 = float(np.sum(Pphi * ys ** 2))
        ipk = int(np.argmax(Pphi))
        Peq_1d = np.exp(-(sp.s(ys) - sp.s(ys).min()))
        Peq_1d /= Peq_1d.sum()
        S = survival(A, (50.0, 60.0))
        i0 = int(np.argmin(np.abs(ys - V / H)))
        row = {"H_GeV": H, "grid": [Ny, Nr], "Lambda1_2D": L2,
               "marginal_phi2_over_H2_with_X": phi2, "phi2_over_H2_X0_slice": float(np.sum(Peq_1d * ys ** 2)),
               "marginal_peak_y": float(ys[ipk]), "X0_slice_peak_y": float(ys[int(np.argmax(Peq_1d))]),
               "rho2_over_H2_mean": float(np.sum(pi.reshape(Ny, Nr) * (rs ** 2)[None, :])),
               "symmetry": "exact: P(phi) = P(-phi) by Z2 (the odd sector is solved explicitly)"}
        for N in (50.0, 60.0):
            m = 0.5 * (1 - S[N])
            order = np.argsort(m)
            cdf = np.cumsum(pi[order])
            q = lambda x: float(m[order][np.searchsorted(cdf, x)])
            y0 = V / H
            pf = float(m[idx[i0, 0]]) if y0 >= 2 * hy else float(0.5 - (0.5 - m[idx[0, 0]]) * y0 / (hy / 2))
            row["N%d" % N] = {"p_minus_fresh_start_phi_v_X0": pf,
                              "eq_draws_mean_minority_fraction": float(np.sum(pi * m)),
                              "eq_draws_minority_quantiles_10_50_90": [q(0.1), q(0.5), q(0.9)],
                              "eq_draws_prob_minority_below_pc": float(np.sum(pi[m < PC_SITE_SC]))}
        if H == 1e13:
            Sb = survival(A, (50.0,), dt=0.025)
            row["dt_check_max_abs_dS_N50"] = float(np.max(np.abs(Sb[50.0] - S[50.0])))
        out["high_H"].append(row)
        print(row, flush=True)
    out["runtime_s"] = time.time() - t0
    dump("stochastic_2d.json", out)
    return out


if __name__ == "__main__":
    main()
