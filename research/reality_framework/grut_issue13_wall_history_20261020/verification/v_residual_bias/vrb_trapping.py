"""Verifier v_residual_bias, part 3: route B (inflationary prior + post-inflation trapping map).

Independent method (no module code imported):
  * RG: one-loop betas re-typed from grut_candidate_1_quantum_20261018/rg_beta_functions.txt [Q]; the xi betas
    derived here as the mirror of the mass betas, 16 pi^2 beta_xi_ab = sum_cd Lambda_abcd (xi_cd - delta_cd/6).
  * FP: Starobinsky-Yokoyama, dP/dN = D d/dy (dP/dy + s' P), D = 1/(8 pi^2), s = (8 pi^2/3)(a2 y^2/2 + lam y^4/4).
    Solved by the SPECTRAL method: P = e^{-s/2} psi, -psi'' + (s'^2/4 - s''/2) psi = (Lambda/D) psi on a fine
    symmetric node grid (Dirichlet), eigh_tridiagonal; P_end = e^{-(s-s0)/2} sum_n u_n u_n(y0) e^{-Lambda_n N}/h.
    The start y0 = v/H sits exactly on a node.
  * Trapping map: vectorized fixed-step RK4 in conformal time for chi = a phi (a = tau, RD from H(1) = 1),
        chi'' = -lam chi^3 + lam u^2 tau^2 chi,  chi(1) = y, chi'(1) = y,  u = v/H,
    trapped when the physical energy phidot^2/2 + lam/4 (phi^2 - u^2)^2 < lam u^4/4; final sign = sign(phi) then.
    Thresholds by vectorized bisection.  Step halving checks the thresholds.
  * m0 = Int_0^inf P_odd(a) s(a) da with P_odd = P(a) - P(-a), integrated with a cubic spline between exact
    thresholds.
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)
V = 1000.0
LAMPOST = 0.01
BETA = 8 * math.pi ** 2 / 3
DIFF = 1 / (8 * math.pi ** 2)


# ------------------------------------------------------------------ RG
def rg_rhs(t, p):
    lp, lx, ld, c, d, al, mp, mx, md, xp, xx, xd = p
    k = 1 / (16 * math.pi ** 2)
    s6 = 1 / 6
    return k * np.array([
        18 * lp ** 2 + 4 * c ** 2 + 2 * d ** 2,
        20 * lx ** 2 + 2 * c ** 2 + 2 * al ** 2,
        18 * ld ** 2 + 2 * d ** 2 + 4 * al ** 2,
        6 * lp * c + 8 * lx * c + 8 * c ** 2 + 2 * d * al,
        6 * lp * d + 6 * ld * d + 4 * c * al + 8 * d ** 2,
        8 * lx * al + 6 * ld * al + 2 * c * d + 8 * al ** 2,
        6 * mp * lp + 4 * mx * c + 2 * md * d,
        2 * mp * c + 8 * mx * lx + 2 * md * al,
        2 * mp * d + 4 * mx * al + 6 * md * ld,
        6 * lp * (xp - s6) + 4 * c * (xx - s6) + 2 * d * (xd - s6),
        2 * c * (xp - s6) + 8 * lx * (xx - s6) + 2 * al * (xd - s6),
        2 * d * (xp - s6) + 4 * al * (xx - s6) + 6 * ld * (xd - s6)])


P0 = [0.01, 0.001, 0.001, 0.01, 1e-4, 0.0, -0.01 * V ** 2, 0.0, 1e4, 0.0, 0.0, 0.0]


def couplings(H):
    if H <= 100.0:
        p = P0
    else:
        sol = solve_ivp(rg_rhs, (0, math.log(H / 100.0)), P0, method="LSODA", rtol=1e-11, atol=1e-15)
        p = sol.y[:, -1]
    lam, mp, xp = p[0], p[6], p[9]
    return dict(lam=float(lam), a2=float(mp / H ** 2 + 12 * xp), xi=float(xp))


# ------------------------------------------------------------------ FP spectral
class FP:
    def __init__(self, H, h_target=0.002, Y=7.0, Emax=400.0):
        cp = couplings(H)
        self.cp = cp
        self.u = V / H
        k = max(1, int(round(self.u / h_target)))
        self.h = self.u / k                       # y0 = u is node k
        n = int(math.ceil(Y / self.h))
        self.y = np.arange(-n, n + 1) * self.h
        self.i0 = n + k
        self.iz = n
        lam, a2 = cp["lam"], cp["a2"]
        y = self.y
        self.s = BETA * (0.5 * a2 * y ** 2 + 0.25 * lam * y ** 4)
        sp = BETA * (a2 * y + lam * y ** 3)
        spp = BETA * (a2 + 3 * lam * y ** 2)
        Vq = sp ** 2 / 4 - spp / 2
        dg = 2 / self.h ** 2 + Vq
        off = -np.ones(len(y) - 1) / self.h ** 2
        E, U = eigh_tridiagonal(dg, off, select="v", select_range=(-1e3, Emax))
        self.Lam = E * DIFF
        self.U = U

    def P_fresh(self, N):
        c = self.U[self.i0, :] * np.exp(-self.Lam * N)
        return np.exp(-(self.s - self.s[self.i0]) / 2) * (self.U @ c) / self.h

    def P_relaxed(self, N):
        w = np.where(self.y > 0, 1.0, np.where(self.y == 0, 0.5, 0.0))
        smin = self.s.min()
        P0 = w * np.exp(-(self.s - smin))
        P0 /= P0.sum() * self.h
        phi0 = np.exp(self.s / 2 - smin / 2) * P0 * np.exp(smin / 2)   # e^{s/2} P0
        c = self.U.T @ phi0 * np.exp(-self.Lam * N)
        return np.exp(-self.s / 2) * (self.U @ c)

    def mean_slope(self, N):
        """d<y_N>/dy0 at y0 = 0 (odd linear response), and <y_N^2> at y0 = 0."""
        e = np.exp(-self.Lam * N)
        dpsi0 = (self.U[self.iz + 1, :] - self.U[self.iz - 1, :]) / (2 * self.h)
        # E[f(y_N)|y0] = e^{s(y0)/2} sum_n psi_n(y0) e^{-Lam N} Int f e^{-s/2} psi_n  (continuum psi = u/sqrt h)
        Iy = (self.U.T @ (self.y * np.exp(-self.s / 2)))
        Iy2 = (self.U.T @ (self.y ** 2 * np.exp(-self.s / 2)))
        slope = float(np.sum(dpsi0 * e * Iy))        # s'(0) = 0 so d/dy0 e^{s/2} term vanishes at 0
        m2 = float(np.sum(self.U[self.iz, :] * e * Iy2))
        return slope, m2


# ------------------------------------------------------------------ trapping map (vectorized RK4)
def final_sign(a, u, lam=LAMPOST, dtau=None, tmax_fac=4.0):
    a = np.atleast_1d(np.asarray(a, float))
    f0 = a * (1 + 2 / (lam * a ** 2)) ** 0.25
    if dtau is None:
        dtau = 0.02 * 7.4163 / (math.sqrt(lam) * f0.max())
    tau = 1.0
    chi = a.copy()
    pc = a.copy()
    sign = np.zeros_like(a)
    alive = np.ones(a.shape, bool)
    tend = max(50.0, tmax_fac * f0.max() / u)

    def acc(c, t):
        return -lam * c ** 3 + lam * u * u * t * t * c
    while alive.any() and tau < tend:
        c, p = chi[alive], pc[alive]
        k1c, k1p = p, acc(c, tau)
        k2c, k2p = p + 0.5 * dtau * k1p, acc(c + 0.5 * dtau * k1c, tau + 0.5 * dtau)
        k3c, k3p = p + 0.5 * dtau * k2p, acc(c + 0.5 * dtau * k2c, tau + 0.5 * dtau)
        k4c, k4p = p + dtau * k3p, acc(c + dtau * k3c, tau + dtau)
        c = c + dtau / 6 * (k1c + 2 * k2c + 2 * k3c + k4c)
        p = p + dtau / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        tau += dtau
        chi[alive], pc[alive] = c, p
        phi = c / tau
        phid = (p / tau - c / tau ** 2) / tau
        E = 0.5 * phid ** 2 + 0.25 * lam * (phi ** 2 - u * u) ** 2
        tr = E < 0.25 * lam * u ** 4
        if tr.any():
            idx = np.where(alive)[0][tr]
            sign[idx] = np.sign(phi[tr])
            alive[idx] = False
    if alive.any():
        raise RuntimeError("untrapped trajectories")
    return sign


def thresholds(u, amax, dfac=1.0, nbis=36):
    da = (u / 0.313) / 10
    ag = np.arange(0.3 * u, amax + da, da)
    sg = np.concatenate([final_sign(ag[i:i + 400], u, dtau=None if dfac == 1 else None) for i in range(0, len(ag), 400)])
    if dfac != 1.0:
        sg = np.concatenate([final_sign(ag[i:i + 400], u, dtau=0.02 * 7.4163 / (0.1 * (ag[i:i + 400] * (1 + 200 / ag[i:i + 400] ** 2) ** 0.25).max()) / dfac)
                             for i in range(0, len(ag), 400)])
    flips = np.where(sg[1:] != sg[:-1])[0]
    lo, hi = ag[flips].copy(), ag[flips + 1].copy()
    slo = sg[flips].copy()
    for _ in range(nbis):
        mid = 0.5 * (lo + hi)
        if dfac == 1.0:
            sm = final_sign(mid, u)
        else:
            f0 = mid * (1 + 2 / (LAMPOST * mid ** 2)) ** 0.25
            sm = final_sign(mid, u, dtau=0.02 * 7.4163 / (0.1 * f0.max()) / dfac)
        same = sm == slo
        lo = np.where(same, mid, lo)
        hi = np.where(same, hi, mid)
    return 0.5 * (lo + hi), (sg[0] if len(sg) else 1.0)


def m0_from(fp, P, th):
    y = fp.y
    pos = y >= 0
    a = y[pos]
    Pp = P[pos]
    Pm = P[::-1][pos]          # P(-a) on the same nodes (grid symmetric)
    Podd = Pp - Pm
    cs = CubicSpline(a, Podd)
    pts = np.concatenate([[0.0], th[th < a[-1]], [a[-1]]])
    tot = 0.0
    for k in range(len(pts) - 1):
        tot += (-1) ** k * cs.integrate(pts[k], pts[k + 1])
    end = cs.integrate(0.0, a[-1])
    return float(tot), float(end)


def emulate_cellcentre(H, N, th, h=0.02, smax=40.0, y0cell=None):
    """inflation_spectator-like: SqRA chain on cell centres (k+1/2)h, delta start in a chosen cell, cell-centre signs."""
    cp = couplings(H)
    lam, a2 = cp["lam"], cp["a2"]
    s_of = lambda y: BETA * (0.5 * a2 * y ** 2 + 0.25 * lam * y ** 4)
    ym = 1.0
    s0 = s_of(np.linspace(0, 5, 5001)).min()
    while s_of(ym) - s0 < smax:
        ym *= 1.15
    M = int(math.ceil(ym / h))
    y = (np.arange(-M, M) + 0.5) * h
    s = s_of(y)
    k = DIFF / h ** 2
    up = np.zeros_like(s)
    dn = np.zeros_like(s)
    up[:-1] = k * np.exp(-(s[1:] - s[:-1]) / 2)
    dn[1:] = k * np.exp(-(s[:-1] - s[1:]) / 2)
    # symmetrize: detailed balance pi_i up_i = pi_{i+1} dn_{i+1}; similarity with sqrt(pi)
    # use dense expm on the small chain via eigen-decomposition of the symmetrized generator
    pi = np.exp(-(s - s.min()))
    sq = np.sqrt(pi)
    d = -(up + dn)
    off = up[:-1] * sq[:-1] / sq[1:]
    E, W = eigh_tridiagonal(d, off)
    u0 = V / H
    i0 = int(np.argmin(np.abs(y - u0))) if y0cell is None else int(np.argmin(np.abs(y - y0cell)))
    p0 = np.zeros_like(y)
    p0[i0] = 1.0
    # forward: p(N) = S^{1/2} W e^{EN} W^T S^{-1/2} p0
    p = sq * (W @ (np.exp(E * N) * (W.T @ (p0 / sq))))
    p = np.clip(p, 0, None)
    a = np.abs(y)
    kk = np.searchsorted(th, a)
    sgn = np.sign(y) * (-1.0) ** kk
    return dict(start_cell=float(y[i0]), p_minus_cellcentre=float(p[sgn < 0].sum() / p.sum()),
                p_minus_end=float(p[y < 0].sum() / p.sum()))


def run_H(H, Ns=(50.0, 60.0), h_target=0.002, check_steps=False):
    t0 = time.time()
    fp = FP(H, h_target=h_target)
    th, s_first = thresholds(fp.u, 6.0)
    row = dict(H=H, u=fp.u, lam=fp.cp["lam"], a2=fp.cp["a2"], xi_phi=fp.cp["xi"],
               y_min=(math.sqrt(-fp.cp["a2"] / fp.cp["lam"]) if fp.cp["a2"] < 0 else 0.0),
               h=fp.h, n_eig=int(len(fp.Lam)), Lambda1=float(fp.Lam[1]), n_thresholds=int(len(th)),
               first_threshold_over_u=float(th[0] / fp.u) if len(th) else None,
               spacing_over_u=(float(np.mean(np.diff(th[-6:])) / fp.u) if len(th) > 6 else None),
               thresholds_first5=th[:5].tolist(), sign_below_first=float(s_first))
    for N in Ns:
        mf, endf = m0_from(fp, fp.P_fresh(N), th)
        mr, endr = m0_from(fp, fp.P_relaxed(N), th)
        row["N%d" % N] = dict(m0_fresh=mf, p_minus_fresh=0.5 * (1 - mf), p_minus_end_fresh=0.5 * (1 - endf),
                              m0_relaxed=mr, m0_end_relaxed=endr)
    if check_steps:
        th2, _ = thresholds(fp.u, 6.0, dfac=2.0)
        n = min(len(th), len(th2))
        row["threshold_maxdiff_dtau_half"] = float(np.max(np.abs(th[:n] - th2[:n]))) if n else None
        mf2, _ = m0_from(fp, fp.P_fresh(50.0), th2)
        row["m0_fresh_N50_dtau_half"] = mf2
        fp2 = FP(H, h_target=2 * h_target)
        row["m0_fresh_N50_h_double"] = m0_from(fp2, fp2.P_fresh(50.0), th)[0]
        row["m0_relaxed_N50_h_double"] = m0_from(fp2, fp2.P_relaxed(50.0), th)[0]
    row["wall_s"] = time.time() - t0
    return row, th


def main(which):
    out = {}
    if which in ("rows", "all"):
        rows = []
        for H in (1400.0, 2000.0, 3000.0, 5000.0, 1e4, 2e4):
            r, th = run_H(H, check_steps=H in (5000.0, 1e4))
            if H in (5000.0, 1e4):
                r["emulation_h0.02_cellcentre"] = {
                    "N50_nearest_cell": emulate_cellcentre(H, 50.0, th),
                    "N50_start_0.11" if H == 1e4 else "N50_start_0.21": emulate_cellcentre(
                        H, 50.0, th, y0cell=(0.11 if H == 1e4 else 0.21)),
                    "N50_start_0.09" if H == 1e4 else "N50_start_0.19": emulate_cellcentre(
                        H, 50.0, th, y0cell=(0.09 if H == 1e4 else 0.19)),
                    "N60_start_0.11" if H == 1e4 else "N60_start_0.21": emulate_cellcentre(
                        H, 60.0, th, y0cell=(0.11 if H == 1e4 else 0.21)),
                    "N50_h0.005_start_exact": emulate_cellcentre(H, 50.0, th, h=0.005,
                                                                 y0cell=(0.1025 if H == 1e4 else 0.2025)),
                }
            rows.append(r)
            print(json.dumps({k: v for k, v in r.items() if k not in ("thresholds_first5",)}), flush=True)
        out["rows"] = rows
        with open(os.path.join(OUT, "trapping.json"), "w") as fh:
            json.dump(out, fh, indent=1)
    if which in ("mean", "all"):
        cm = []
        for H in (1e4, 1e5, 1e8):
            fp = FP(H, h_target=0.002)
            sl, m2 = fp.mean_slope(50.0)
            u = V / H
            cm.append(dict(H=H, mean_phi_over_v=sl, mean_over_rms=sl * u / math.sqrt(m2)))
            print(cm[-1], flush=True)
        with open(os.path.join(OUT, "coherent_mean.json"), "w") as fh:
            json.dump(cm, fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
