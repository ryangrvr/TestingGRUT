"""Per-e-fold production of minority (phi < 0) volume during the last N_* e-folds, from the inflation_spectator
Fokker-Planck chain (1D, X = D = 0 along the path, couplings and xi run to mu = H). [P] inflation is preparation data.

Identity used [derived, NOTES §1]. The phi-odd sector is the half chain y > 0 with killing rate kappa_0 = 2D/h^2
in cell 0 (inflation_spectator fp1d.odd_chain). For a start in cell m0 the minority fraction is p_-(n) = F(n)/2 with
F the absorbed probability, so

    dp_-/dn (fresh, delta at m0) = (kappa_0/2) P(X_n = 0 | X_0 = m0)
                                 = (kappa_0/2) (pi_0/pi_m0) P(X_n = m0 | X_0 = 0)        (reversibility)
    dp_-/dn (relaxed, start ~ pi restricted to y > 0) = (kappa_0/2) (pi_0/Z_+) sum_m P(X_n = m | X_0 = 0)

ONE forward solve from the barrier cell 0 gives both. P(X_n = m0 | X_0 = 0) is the downhill probability, O(1e-2)
after the roll time, so the exponentially small e^{-B} sits in the closed-form Boltzmann factor and no small
number is formed by subtraction. The forward ODE dp/dn = Q^T p is solved with BDF and component-wise error control
(atol 1e-120; 1e-80 and 1e-140 agree to 4 decimals in log10), so the deep-lag tail (fresh start) is resolved relative to itself.

Verifier corrections applied (v_inflation): the fresh start is placed EXACTLY at y = v/H (grid spacing adjusted so
that a cell centre sits at v/H; the module's start-cell offset biased p_- up by <= 2 %), Lambda_1 is reported from
the eigenvalue where resolvable (MFPT only at large barrier where both agree).
"""
import math
import time
import numpy as np
from scipy.sparse import diags
from scipy.integrate import solve_ivp
from scipy.special import logsumexp
from scipy.linalg import eigh_tridiagonal
from mr_common import Spectator, spec_fp1d, V, dump

N_MAX = 62.0
DN = 0.05


def build(H):
    sp = Spectator(H)
    yy0, h0 = spec_fp1d.make_grid(sp)
    y0 = V / H
    m0 = max(2, int(round(y0 / h0 - 0.5)))
    h = y0 / (m0 + 0.5)                       # cell centre exactly at v/H (verifier correction)
    yy, h = spec_fp1d.make_grid(sp, h=h)
    assert abs(yy[m0] - y0) < 1e-12 * max(1.0, y0)
    s, up, dn, kill, k = spec_fp1d.odd_chain(sp, yy, h)
    return sp, yy, h, m0, s, up, dn, kill


def lambda1(sp, yy, h):
    s, up, dn, kill, k = spec_fp1d.odd_chain(sp, yy, h)
    lam = eigh_tridiagonal(up + dn + kill, -k * np.ones(len(s) - 1), select="i", select_range=(0, 0))[0][0]
    L1m, _ = spec_fp1d.lambda1_mfpt(sp, yy, h)
    return float(lam), float(L1m)


def production(H, nmax=N_MAX, dn_=DN, atol=1e-120, rtol=1e-9):
    t0 = time.time()
    sp, yy, h, m0, s, up, dn, kill = build(H)
    M = len(s)
    # forward generator transpose: dp_j/dn = up[j-1] p_{j-1} + dn[j+1] p_{j+1} - (up_j + dn_j + kill_j) p_j
    QT = diags([up[:-1], -(up + dn + kill), dn[1:]], [-1, 0, 1], format="csc")
    p0 = np.zeros(M)
    p0[0] = 1.0
    ns = np.round(np.arange(0.0, nmax + 1e-9, dn_), 10)
    sol = solve_ivp(lambda t, p: QT @ p, (0.0, nmax), p0, method="BDF", jac=QT, rtol=rtol, atol=atol, t_eval=ns)
    if not sol.success:
        raise RuntimeError(sol.message)
    P = sol.y
    pm0 = np.clip(P[m0], 0.0, None)
    Ssum = np.clip(P.sum(axis=0), 0.0, None)
    log_pref_fresh = math.log(0.5 * kill[0]) - (s[0] - s[m0])
    logZp = logsumexp(-s)
    log_pref_rel = math.log(0.5 * kill[0]) - s[0] - logZp
    with np.errstate(divide="ignore"):
        ldf = log_pref_fresh + np.log(pm0)
        ldr = log_pref_rel + np.log(Ssum)
    # integrated p_-(N) (trapezoid in n, relative-precision because integrand is positive)
    def cum(ld):
        f = np.exp(ld)
        return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(ns))])
    pf = cum(ldf)
    pr = cum(ldr)
    lam_eig, lam_mfpt = lambda1(sp, yy, h)
    i50, i60 = int(round(50 / dn_)), int(round(60 / dn_))
    return {"H_GeV": H, "B_HM": sp.barrier_B(), "v_over_H": V / H, "grid_h": h, "cells": M, "m0": m0,
            "Lambda1_eig": lam_eig, "Lambda1_mfpt": lam_mfpt,
            "Lambda1_adopted": lam_eig if lam_eig > 1e-12 else lam_mfpt,
            "log10_dpdn_fresh": np.round(ldf / math.log(10), 5).tolist(),
            "log10_dpdn_relaxed": np.round(ldr / math.log(10), 5).tolist(),
            "p_minus_fresh_N50": pf[i50], "p_minus_fresh_N60": pf[i60],
            "p_minus_relaxed_N50": pr[i50], "p_minus_relaxed_N60": pr[i60],
            "runtime_s": time.time() - t0, "n_grid": {"dn": dn_, "nmax": nmax}}


def H_grid():
    return np.unique(np.round(np.concatenate([np.logspace(np.log10(130.0), np.log10(500.0), 49),
                                              np.logspace(np.log10(550.0), np.log10(2000.0), 8)]), 4))


def main():
    t0 = time.time()
    from multiprocessing import Pool
    with Pool(2) as pool:                      # 2 single-threaded workers (BRIEF: <= 2 threads)
        rows = pool.map(production, [float(H) for H in H_grid()])
    out = {"module": "minority_relic", "part": "per-e-fold minority production dp_-/dn, fresh and relaxed starts",
           "labels": "computed with inflation_spectator's SqRA chain (imported read-only) [Q-code]; identities derived "
                     "here; inflation/H_inf/N_*/start condition are outside-action preparation data [P]",
           "n_axis": "n = e-folds after the start of the last-N_* window (fresh start at phi = +v exactly)",
           "rows": rows, "runtime_s": time.time() - t0}
    dump("production.json", out)
    return out


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        for H in (250.0, 335.0, 1000.0, 2000.0, 200.0, 150.0):
            r = production(H)
            print("%7.1f B=%7.2f M=%d  L1e=%.4e L1m=%.4e pf50=%.5e pf60=%.5e pr50=%.5e pr60=%.5e  t=%.1fs" % (
                H, r["B_HM"], r["cells"], r["Lambda1_eig"], r["Lambda1_mfpt"], r["p_minus_fresh_N50"],
                r["p_minus_fresh_N60"], r["p_minus_relaxed_N50"], r["p_minus_relaxed_N60"], r["runtime_s"]))
    else:
        main()
