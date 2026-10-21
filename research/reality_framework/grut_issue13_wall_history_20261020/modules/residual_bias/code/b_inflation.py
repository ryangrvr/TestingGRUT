"""Route B: randomized prior from inflation with H_inf >> v [P].  Systematic bias left by the post-inflation
trapping map, and the coherent ensemble mean that would enter a later thermal restoration (route A).

Model (same as inflation_spectator, re-implemented here, not imported):
  * 1D Starobinsky-Yokoyama Fokker-Planck for y = phi/H with X = D = 0, couplings and xi_phi run to mu = H
    (one-loop betas of rg_beta_functions.txt as transcribed in inflation_spectator/code/common.py, re-typed here),
    w = a2 y^2/2 + lam y^4/4, a2 = m_p^2(mu)/H^2 + 12 xi_phi(mu); P_eq ~ exp(-s), s = (8 pi^2/3) w.
  * SqRA reversible chain on a symmetric cell-centred grid, rates (D/h^2) exp(-(s_j - s_i)/2), D = 1/(8 pi^2).
  * Instant radiation era after inflation, phi frozen at the end of inflation, chi = a phi:
        chi'' + lam chi^3 - lam (v/H)^2 tau^2 chi = 0, chi(1) = chi'(1) = y_end   (units H_inf = 1, a = tau),
    integrated until the physical energy falls below the barrier lam v^4/4 ("trapped"); final sign F(y_end).
Key identity [D]: F is odd, F(-y) = -F(y) (Z2 of the ODE).  Hence
    m0 = E[F(y_end)] = Int_0^inf [P_end(a) - P_end(-a)] s(a) da,   s(a) = F(a) for a > 0,
so a symmetric end distribution gives m0 = 0 exactly, and any nonzero m0 is the overlap of the ODD part of P_end
(present for a fresh +v start or a relaxed + start) with the oscillating trapping map s(a).
Computation: backward Kolmogorov equation m(y, n) = [exp(n Q) F](y) on the grid, with F as exact cell averages
of the threshold map (thresholds a_k located by bisection), so no cell-centre sign assignment error.
"""
import math
import time
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from common import V, dump, log10_m0

warnings.filterwarnings("ignore", category=RuntimeWarning)

MU0 = 100.0
K16 = 16.0 * math.pi ** 2
DIFF = 1.0 / (8.0 * math.pi ** 2)
BETA_FP = 8.0 * math.pi ** 2 / 3.0
LAM_POST = 0.01       # post-inflation quartic (as in inflation_spectator post_inflation.py)
P0 = np.array([0.01, 0.001, 0.001, 0.01, 1e-4, 0.0, -0.01 * V ** 2, 0.0, (0.1 * V) ** 2, 0.0, 0.0, 0.0])


def betas(p):
    """16 pi^2 beta for (lp, lx, ld, c, d, alpha, mp, mx, md, xp, xx, xd) [Q rg_beta_functions.txt via
    inflation_spectator/code/common.py B_all; re-typed]."""
    lp, lx, ld, c, d, a, mp, mx, md, xp, xx, xd = p
    s = 1.0 / 6.0
    return np.array([
        18 * lp ** 2 + 4 * c ** 2 + 2 * d ** 2,
        20 * lx ** 2 + 2 * c ** 2 + 2 * a ** 2,
        18 * ld ** 2 + 2 * d ** 2 + 4 * a ** 2,
        6 * lp * c + 8 * lx * c + 8 * c ** 2 + 2 * d * a,
        6 * lp * d + 6 * ld * d + 4 * c * a + 8 * d ** 2,
        8 * lx * a + 6 * ld * a + 2 * c * d + 8 * a ** 2,
        6 * mp * lp + 4 * mx * c + 2 * md * d,
        2 * mp * c + 8 * mx * lx + 2 * md * a,
        2 * mp * d + 4 * mx * a + 6 * md * ld,
        6 * lp * (xp - s) + 4 * c * (xx - s) + 2 * d * (xd - s),
        2 * c * (xp - s) + 8 * lx * (xx - s) + 2 * a * (xd - s),
        2 * d * (xp - s) + 4 * a * (xx - s) + 6 * ld * (xd - s),
    ])


_RUN = solve_ivp(lambda t, y: betas(y) / K16, (0.0, math.log(1e15 / MU0)), P0, method="DOP853", rtol=1e-12,
                 atol=1e-16, dense_output=True)


def spectator(H):
    p = _RUN.sol(math.log(H / MU0)) if H >= MU0 else P0
    lp, mp, xp = p[0], p[6], p[9]
    a2 = mp / H ** 2 + 12.0 * xp
    return dict(H=H, lam=lp, a2=a2, xi_phi=xp, vH=V / H)


def s_of(sp, y):
    return BETA_FP * (0.5 * sp["a2"] * y ** 2 + 0.25 * sp["lam"] * y ** 4)


# ------------------------------------------------------------------------------------------------ trapping map
def n_cross(y_end, vH, lam=LAM_POST):
    """Number of zero crossings before trapping (final sign = (-1)^N for y_end > 0)."""
    f0 = y_end * (1.0 + 2.0 / (lam * y_end ** 2)) ** 0.25
    tmax = max(50.0, 4.0 * f0 / vH)
    period = 7.4163 / (math.sqrt(lam) * f0)

    def rhs(t, u):
        return [u[1], -lam * u[0] ** 3 + lam * vH ** 2 * t * t * u[0]]

    def zero(t, u):
        return u[0]

    def trapped(t, u):
        phi = u[0] / t
        phidot = (u[1] - u[0] / t) / t ** 2
        E = 0.5 * phidot ** 2 + 0.25 * lam * (phi ** 2 - vH ** 2) ** 2
        return E - 0.25 * lam * vH ** 4 * (1 - 1e-12)
    trapped.terminal = True
    trapped.direction = -1
    sol = solve_ivp(rhs, (1.0, tmax), [y_end, y_end], method="DOP853", rtol=1e-10, atol=1e-12,
                    events=(zero, trapped), max_step=period / 20)
    return len(sol.t_events[0])


def thresholds(vH, amax, lam=LAM_POST):
    """Positions a_k (in units of H) where the final sign of a positive start flips (parity of N_cross changes)."""
    da = (vH / 0.313) / 8.0
    a_grid = np.arange(0.5 * vH, amax + da, da)
    par = np.array([n_cross(a, vH, lam) % 2 for a in a_grid])
    th = []
    for i in np.where(par[1:] != par[:-1])[0]:
        lo, hi = a_grid[i], a_grid[i + 1]
        plo = par[i]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if n_cross(mid, vH, lam) % 2 == plo:
                lo = mid
            else:
                hi = mid
        th.append(0.5 * (lo + hi))
    return np.array(th), a_grid, par


def cell_average_F(y, h, th):
    """Exact cell averages of F(y) = sgn(y) * (-1)^{#thresholds below |y|}."""
    def G(x):   # integral of s(a) from 0 to x >= 0
        x = np.asarray(x, dtype=float)
        pts = np.concatenate([[0.0], th])
        out = np.zeros_like(x)
        for k in range(len(pts)):
            left = pts[k]
            right = pts[k + 1] if k + 1 < len(pts) else np.inf
            sign = (-1) ** k
            out += sign * np.clip(np.minimum(x, right) - left, 0.0, None)
        return out
    lo = y - h / 2
    hi = y + h / 2
    # F odd: Int_lo^hi F = sgn-aware via G(|.|)
    def I(x):          # Int_0^x F(y) dy  (even function of x, since F odd)
        return G(np.abs(x))
    return (I(hi) - I(lo)) / h


def sign_cell_center(y, th):
    a = np.abs(y)
    k = np.searchsorted(th, a)
    return np.sign(y) * (-1.0) ** k


# ------------------------------------------------------------------------------------------------ FP chain
def make_grid(sp, h, smax=60.0):
    y = 1.0
    s0 = min(s_of(sp, np.linspace(0, 5, 5001)))
    while s_of(sp, y) - s0 < smax:
        y *= 1.1
    M = int(math.ceil(y / h))
    yy = (np.arange(-M, M) + 0.5) * h
    return yy


def generator(sp, yy, h):
    s = s_of(sp, yy)
    k = DIFF / h ** 2
    up = np.zeros_like(s)
    dn = np.zeros_like(s)
    up[:-1] = k * np.exp(-(s[1:] - s[:-1]) / 2)
    dn[1:] = k * np.exp(-(s[:-1] - s[1:]) / 2)
    Q = diags([dn[1:], -(up + dn), up[:-1]], [-1, 0, 1], format="csr")   # backward: (Qf)_i = sum_j k_ij (f_j-f_i)
    return Q, s


def odd_interp(yy, m, y0):
    return float(np.interp(y0, yy, m))


def analyse(H, h, Ns=(50.0, 60.0), th=None):
    sp = spectator(H)
    yy = make_grid(sp, h)
    Q, s = generator(sp, yy, h)
    amax = yy[-1]
    if th is None:
        th, _, _ = thresholds(sp["vH"], amax)
    F_avg = cell_average_F(yy, h, th)
    F_ctr = sign_cell_center(yy, th)
    F_end = np.sign(yy)
    pi = np.exp(-(s - s.min()))
    pos = yy > 0
    out = dict(H_GeV=H, h=h, vH=sp["vH"], lam=sp["lam"], xi_phi=sp["xi_phi"], a2=sp["a2"],
               n_thresholds=len(th), first_threshold_over_vH=(th[0] / sp["vH"] if len(th) else None))
    for N in Ns:
        mA = expm_multiply(Q * N, F_avg)
        mC = expm_multiply(Q * N, F_ctr)
        mE = expm_multiply(Q * N, F_end)
        y0 = sp["vH"]
        rel = lambda m: float((pi[pos] * m[pos]).sum() / pi[pos].sum())
        eq_rms = float(np.sqrt((pi * mA ** 2).sum() / pi.sum()))
        # realized-patch spread [E]: one e-fold after the fresh start, sub-patch values y1 ~ G(y0 -> ., 1)
        mA1 = expm_multiply(Q * (N - 1.0), F_avg)
        d0 = np.zeros_like(yy)
        d0[int(np.argmin(np.abs(yy - y0)))] = 1.0
        P1 = expm_multiply(Q.T * 1.0, d0)
        P1 = np.clip(P1, 0, None)
        P1 /= P1.sum()
        mean1 = float((P1 * mA1).sum())
        sd1 = float(np.sqrt(max((P1 * mA1 ** 2).sum() - mean1 ** 2, 0.0)))
        out["N%d" % N] = dict(
            m0_fresh_after_trapping=odd_interp(yy, mA, y0),
            m0_fresh_cellcentre_sign=odd_interp(yy, mC, y0),
            m0_fresh_end_of_inflation=odd_interp(yy, mE, y0),
            p_minus_fresh_after_trapping=0.5 * (1 - odd_interp(yy, mA, y0)),
            m0_relaxed_after_trapping=rel(mA),
            m0_relaxed_end_of_inflation=rel(mE),
            m_patch_rms_equilibrium_draw=eq_rms,
            fresh_one_efold_subpatch_mean=mean1, fresh_one_efold_subpatch_sd=sd1)
    return out, th


def coherent_mean(H, h=0.01, N=50.0):
    """Ensemble mean <phi_end>/v and <phi_end>/phi_rms for a fresh +v start (input to route A)."""
    sp = spectator(H)
    yy = make_grid(sp, h)
    Q, s = generator(sp, yy, h)
    my = expm_multiply(Q * N, yy.copy())
    my2 = expm_multiply(Q * N, yy ** 2)
    y0 = sp["vH"]
    m1 = odd_interp(yy, my, y0)             # odd, linear near 0 -> fine even for y0 << h
    m2 = float(np.interp(y0, yy, my2))
    return dict(H_GeV=H, vH=y0, mean_phi_over_v=m1 / y0, mean_over_rms=m1 / math.sqrt(m2),
                rms_phi_over_H=math.sqrt(m2))


def main():
    t0 = time.time()
    out = {"module": "residual_bias", "part": "B inflationary prior, trapping map, propagation to route A",
           "inputs [P]": "H_inf, N_* = 50/60, start condition (fresh +v / relaxed + / equilibrium draw), instant RD, "
                         "no preheating or restoration (trapping-only regime H_inf <~ 15-20 v)"}
    rows = []
    conv = []
    for H in (700.0, 1000.0, 1400.0, 2000.0, 3000.0, 5000.0, 7000.0, 1e4, 1.4e4, 2e4):
        sp = spectator(H)
        h0 = min(0.01, sp["vH"] / 10)
        r, th = analyse(H, h0)
        r["thresholds_first5"] = th[:5].tolist()
        rows.append(r)
        print("H %.0f done %.1fs  m0 fresh N50 %.3e relaxed %.3e" % (H, time.time() - t0,
              r["N50"]["m0_fresh_after_trapping"], r["N50"]["m0_relaxed_after_trapping"]), flush=True)
        if H in (3000.0, 5000.0, 1e4):
            for fac in (2.0, 0.5):
                rr, _ = analyse(H, h0 / fac, Ns=(50.0,), th=th)
                conv.append(dict(H_GeV=H, h=h0 / fac, m0_fresh=rr["N50"]["m0_fresh_after_trapping"],
                                 m0_relaxed=rr["N50"]["m0_relaxed_after_trapping"],
                                 m0_fresh_cellcentre=rr["N50"]["m0_fresh_cellcentre_sign"]))
    out["rows [C]"] = rows
    out["grid_convergence [C]"] = conv
    # module emulation: their grid h = 0.02 with cell-centre signs, start cell nearest v/H
    emu = []
    for H in (3000.0, 5000.0, 1e4):
        r, _ = analyse(H, 0.02, Ns=(50.0, 60.0))
        emu.append(dict(H_GeV=H, h=0.02, p_minus_cellcentre_N50=0.5 * (1 - r["N50"]["m0_fresh_cellcentre_sign"]),
                        p_minus_cellavg_N50=0.5 * (1 - r["N50"]["m0_fresh_after_trapping"]),
                        p_minus_cellcentre_N60=0.5 * (1 - r["N60"]["m0_fresh_cellcentre_sign"])))
    out["module_emulation_h0.02 [C]"] = emu
    out["inflation_spectator_quoted [Q]"] = {"p_minus_after_trapping_N50": {"3000": 0.5139226569730427,
                                                                            "5000": 0.4999299890642369,
                                                                            "10000": 0.49969144387635966},
                                              "p_minus_after_trapping_N60": {"5000": 0.499696313731129,
                                                                             "10000": 0.49967334917889267}}
    # coherent ensemble mean entering a later restoration (route A input)
    out["coherent_mean_fresh_start [C]"] = [coherent_mean(H) for H in (1e4, 1e5, 1e6, 1e8, 1e10, 1e13)]
    out["runtime_s"] = time.time() - t0
    dump("b_inflation.json", out)
    return out


if __name__ == "__main__":
    o = main()
    for r in o["rows [C]"]:
        print(r["H_GeV"], {k: (round(v, 8) if isinstance(v, float) else v) for k, v in r["N50"].items()})
    print(o["grid_convergence [C]"])
    print(o["module_emulation_h0.02 [C]"])
    print(o["coherent_mean_fresh_start [C]"])
