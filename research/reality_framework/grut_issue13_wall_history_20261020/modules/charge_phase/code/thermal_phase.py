"""Part 3: leading-order finite-temperature phase structure at nonzero X charge.

Units: v = 1 (T/v, n/v^3, f/v^4).  Benchmark couplings.

Free energy density at fixed hidden temperature T and net X charge density n (canonical):
  f(phi; T, n) = U_T(phi) + f_gas(n_g; T, M) + E_cond(n - n_g; M),
  U_T(phi) = (lam_phi/4)(phi^2-1)^2 + 1/2 (lam_phi/4 + hD^2/12) T^2 phi^2      [phi and D loops, O(T^2) only]
  M^2(phi) = gX^2 phi^2 + Pi_X T^2,  Pi_X = gX^2/12 + lam_X/3                     [leading thermal X mass]
  f_gas    = -P(T,mu,M) + mu n_g : exact one-loop ideal Bose gas of X and anti-X with chemical potential mu<=M
             (contains the X-loop thermal mass of phi, gX^2 T^2/6 at mu=0, and its mu dependence)
  E_cond(q;M) = min_c [ q^2/(4c) + M^2 c + lam_X c^2 ]   (homogeneous condensate c=|X_c|^2 carrying charge q)
  n_g = min(n, n_max(T,M)), n_max = n(T, mu=M, M)  (Bose-Einstein condensation when mu reaches M).
The T -> 0 limit is exactly the part-1 problem.  Zero-temperature loop (Coleman-Weinberg) pieces are dropped,
phi/D cubic (daisy) terms are dropped, gas-condensate Hartree terms beyond Pi_X are dropped: LEADING ORDER ONLY.
Run: python3 thermal_phase.py   (about 2-4 minutes, 1 core)
"""
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicHermiteSpline
from scipy.integrate import quad
from common import G_X, LAM_PHI, LAM_X, H_D, V_GEV, save_json

g2 = G_X**2
PI_X = g2 / 12 + LAM_X / 3
PI_P = LAM_PHI / 4 + H_D**2 / 12        # phi self + D loops (X loop is inside the gas)
PI_PHI_TOTAL = LAM_PHI / 4 + g2 / 6 + H_D**2 / 12

# ---------------- ideal Bose gas (particle + antiparticle) ----------------
_bp = [0, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.6, 1.5, 3.0, 6.0, 12.0, 25.0, 60.0]
_gx, _gw = np.polynomial.legendre.leggauss(24)
_Y, _W = [], []
for lo, hi in zip(_bp[:-1], _bp[1:]):
    _Y.append(0.5 * (hi - lo) * _gx + 0.5 * (hi + lo)); _W.append(0.5 * (hi - lo) * _gw)
YB = np.concatenate(_Y); WB = np.concatenate(_W)


def bose(mt, mut):
    """n/T^3, P/T^4, <|X|^2>/T^2 for arrays mt=M/T, mut=mu/T (|mu|<=M). Broadcast over leading dims."""
    mt = np.asarray(mt, float)[..., None]
    mut = np.asarray(mut, float)[..., None]
    sc = np.maximum(1.0, np.sqrt(mt))
    y = sc * YB
    w = sc * WB
    eps = np.sqrt(y * y + mt * mt)
    zm = np.maximum((y * y + (mt - mut) * (mt + mut)) / (eps + mut), 1e-300)   # = eps - mu, cancellation-free
    zp = eps + mut
    with np.errstate(over="ignore", divide="ignore"):
        nbm = 1.0 / np.expm1(zm)
        nbp = 1.0 / np.expm1(zp)
        lm = np.log(-np.expm1(-zm))
        lp = np.log(-np.expm1(-zp))
    pref = 1.0 / (2 * np.pi**2)
    n = pref * np.sum(w * y * y * (nbm - nbp), axis=-1)
    P = -pref * np.sum(w * y * y * (lm + lp), axis=-1)
    X2 = pref * np.sum(w * y * y * (nbm + nbp) / (2 * eps), axis=-1)
    return n, P, X2


def bose_quad(mt, mut):
    E = lambda y: np.sqrt(y * y + mt * mt)
    zm = lambda y: (y * y + (mt - mut) * (mt + mut)) / (E(y) + mut)
    f_n = lambda y: y * y * (1 / np.expm1(zm(y)) - 1 / np.expm1(E(y) + mut))
    f_P = lambda y: -y * y * (np.log(-np.expm1(-zm(y))) + np.log(-np.expm1(-(E(y) + mut))))
    f_X = lambda y: y * y * (1 / np.expm1(zm(y)) + 1 / np.expm1(E(y) + mut)) / (2 * E(y))
    top = 60 * max(1.0, np.sqrt(mt))
    pts = [1e-3, 0.1, 1, 10]
    res = []
    for f in (f_n, f_P, f_X):
        res.append(quad(f, 1e-9, top, points=pts, limit=400, epsabs=1e-14, epsrel=1e-12)[0] / (2 * np.pi**2))
    return res


def quadrature_check():
    rows = []
    for mt, r in [(0.03, 0.0), (0.03, 1.0), (0.3, 0.5), (1.0, 1.0), (5.0, 0.9), (30.0, 1.0)]:
        a = bose(mt, r * mt)
        q = bose_quad(mt, r * mt)
        rows.append(dict(mt=mt, mu_over_M=r, err_rel_or_abs=[float(abs(x - y) / max(1e-12, abs(y))) for x, y in zip(a, q)]))
    return rows


# ---------------- condensate ----------------
def cond(q, M2):
    """c, E_cond, mu_c for charge q>=0 in condensate with X mass^2 M2 (arrays)."""
    q = np.asarray(q, float); M2 = np.asarray(M2, float) * np.ones_like(q)
    ca = np.where(M2 > 0, q / (2 * np.sqrt(np.maximum(M2, 1e-300))), np.inf)
    cb = (q * q / (8 * LAM_X)) ** (1 / 3)
    c = np.where(np.isfinite(ca), np.maximum(ca, cb), cb)
    for _ in range(80):
        F = 8 * LAM_X * c**3 + 4 * M2 * c * c - q * q
        dF = 24 * LAM_X * c * c + 8 * M2 * c
        c = np.where(dF > 0, c - F / np.where(dF > 0, dF, 1.0), c)
    c = np.maximum(c, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        E = np.where(q > 0, q * q / (4 * c) + M2 * c + LAM_X * c * c, 0.0)
        mu = np.where(q > 0, q / (2 * c), np.sqrt(M2))
    return c, E, mu


def UT(phi, T):
    return LAM_PHI / 4 * (phi * phi - 1) ** 2 + 0.5 * PI_P * T * T * phi * phi


# ---------------- per-temperature tables ----------------
PHI = np.concatenate([[0.0], np.linspace(0.004, 1.3, 325)])
NR = 400
_s = np.linspace(0, 1, NR)
RGRID = 1 - (1 - _s) ** 2


class Table:
    def __init__(self, T):
        self.T = T
        self.M2 = g2 * PHI**2 + PI_X * T * T
        self.M = np.sqrt(self.M2)
        if T == 0:
            return
        mt = self.M / T
        mut = mt[:, None] * RGRID[None, :]
        n = np.empty_like(mut); P = np.empty_like(mut)
        for i in range(PHI.size):
            n[i], P[i], _ = bose(np.full(NR, mt[i]), mut[i])
        self.n = n * T**3; self.P = P * T**4; self.mu = mut * T
        self.nmax = self.n[:, -1]

    def f_matrix(self, ngrid):
        """f(phi_i, n_j) for all PHI x ngrid."""
        T = self.T
        out = np.empty((PHI.size, ngrid.size))
        for i in range(PHI.size):
            if T == 0:
                _, E, _ = cond(ngrid, self.M2[i])
                out[i] = UT(PHI[i], 0.0) + E
                continue
            ng = np.minimum(ngrid, self.nmax[i])
            mu = np.interp(ng, self.n[i], self.mu[i])
            Pspl = CubicHermiteSpline(self.mu[i], self.P[i], self.n[i])
            fg = -Pspl(mu) + mu * ng
            _, E, _ = cond(ngrid - ng, self.M2[i])
            out[i] = UT(PHI[i], T) + fg + E
        return out


def refine_min(fcol):
    """global min over PHI with parabolic refinement; returns (phi_min, f_min, index)."""
    i = int(np.argmin(fcol))
    if 0 < i < PHI.size - 1:
        x0, x1, x2 = PHI[i - 1:i + 2]; y0, y1, y2 = fcol[i - 1:i + 2]
        den = (x0 - x1) * (x0 - x2) * (x1 - x2)
        A_ = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den
        B_ = (x2**2 * (y0 - y1) + x1**2 * (y2 - y0) + x0**2 * (y1 - y2)) / den
        C_ = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / den
        if A_ > 0:
            xm = -B_ / (2 * A_)
            return xm, C_ - B_**2 / (4 * A_), i
    return PHI[i], fcol[i], i


# ---------------- local stability of phi = 0 (exact, no interpolation) ----------------
def m2_at_zero(T, n):
    """d^2 f/d phi^2 at phi=0 for fixed (T, n)."""
    M2 = PI_X * T * T
    M = np.sqrt(M2)
    if T == 0:
        c, _, _ = cond(np.array([n]), 0.0)
        return -LAM_PHI + 2 * g2 * c[0]
    nmax = bose(M / T, M / T)[0] * T**3
    if n <= nmax:
        if n == 0:
            mu = 0.0
        else:
            mu = brentq(lambda m: bose(M / T, m / T)[0] * T**3 - n, 0.0, M, xtol=1e-15, rtol=1e-14)
        X2 = bose(M / T, mu / T)[2] * T * T
        c = 0.0; extra = 0.0
    else:
        X2 = bose(M / T, 1.0 * M / T)[2] * T * T
        c, _, muc = cond(np.array([n - nmax]), M2)
        c = c[0]; muc = muc[0]
        dM2 = 1e-6 * M2
        Mp, Mm = np.sqrt(M2 + dM2), np.sqrt(M2 - dM2)
        dn = (bose(Mp / T, Mp / T)[0] - bose(Mm / T, Mm / T)[0]) * T**3 / (2 * dM2)
        extra = (M - muc) * dn
    return -LAM_PHI + PI_P * T * T + 2 * g2 * (X2 + c + extra)


def Tc_zero_charge():
    return brentq(lambda T: m2_at_zero(T, 0.0), 1.0, 2.0, xtol=1e-12)


def local_line(T):
    if m2_at_zero(T, 0.0) > 0:
        return 0.0
    return brentq(lambda n: m2_at_zero(T, n), 0.0, 1.0, xtol=1e-12)


# ---------------- per-T analysis ----------------
NGRID = np.linspace(0.0, 0.40, 4001)


def analyse_T(T):
    tab = Table(T)
    fm = tab.f_matrix(NGRID)
    phimin = np.empty(NGRID.size); fmin = np.empty(NGRID.size); has_broken_local = np.zeros(NGRID.size, bool)
    for j in range(NGRID.size):
        col = fm[:, j]
        pm, fv, i = refine_min(col)
        phimin[j] = pm if i > 0 else 0.0
        fmin[j] = fv if i > 0 else col[0]
        interior = (col[1:-1] < col[:-2]) & (col[1:-1] < col[2:])
        has_broken_local[j] = bool(np.any(interior))
    sym_global = phimin < 2e-3
    # global line: first n above which phi=0 is the global minimum (for all larger n on grid)
    idx = np.where(~sym_global)[0]
    n_glob = 0.0 if idx.size == 0 else (NGRID[idx[-1]] + NGRID[min(idx[-1] + 1, NGRID.size - 1)]) / 2
    idx2 = np.where(has_broken_local)[0]
    n_spin_broken = 0.0 if idx2.size == 0 else NGRID[idx2[-1]]
    # convex hull of F(n) = fmin
    pts = np.column_stack([NGRID, fmin])
    hull = [0]
    for k in range(1, len(pts)):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = pts[hull[-2]], pts[hull[-1]]
            x3, y3 = pts[k]
            if (y2 - y1) * (x3 - x1) >= (y3 - y1) * (x2 - x1) - 1e-18:
                hull.pop()
            else:
                break
        hull.append(k)
    hull = np.array(hull)
    gaps = np.where(np.diff(hull) > 1)[0]
    coex = []
    for gi in gaps:
        lo, hi = NGRID[hull[gi]], NGRID[hull[gi + 1]]
        if hi - lo > 5 * (NGRID[1] - NGRID[0]):
            slope = (fmin[hull[gi + 1]] - fmin[hull[gi]]) / (hi - lo)
            coex.append(dict(n_low=lo, n_high=hi, mu_coex=slope,
                             phi_low=phimin[hull[gi]], phi_high=phimin[hull[gi + 1]]))
    n_loc = local_line(T)
    nmax0 = 0.0 if T == 0 else float(bose(np.sqrt(PI_X), np.sqrt(PI_X))[0] * T**3)
    return dict(T=T, n_local=n_loc, n_global_homogeneous=n_glob, n_broken_spinodal=n_spin_broken,
                coexistence=coex, gas_nmax_at_phi0=nmax0,
                phi_min_at_n0=float(phimin[0]))


if __name__ == "__main__":
    qc = quadrature_check()
    Tc0 = Tc_zero_charge()
    # order of the n=0 transition: scan phi_min(T) around Tc0
    near = []
    for T in np.linspace(Tc0 - 0.02, Tc0 + 0.005, 11):
        tab = Table(T)
        col = tab.f_matrix(np.array([0.0]))[:, 0]
        pm, fv, i = refine_min(col)
        interior = (col[1:-1] < col[:-2]) & (col[1:-1] < col[2:])
        near.append(dict(T=T, phi_min=pm if i > 0 else 0.0, n_local_minima_phi_gt0=int(interior.sum()),
                         m2_phi0=m2_at_zero(T, 0.0)))
    Ts = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.53,
          1.535, 1.54, 1.545, 1.55]
    rows = []
    for T in Ts:
        r = analyse_T(T)
        rows.append(r)
        print(T, r["n_local"], r["n_global_homogeneous"], r["n_broken_spinodal"], r["coexistence"][:1], flush=True)
    # dense local line for interpolation by cosmo_map
    Td = np.linspace(0.0, Tc0, 121)
    nloc_d = [local_line(T) for T in Td]
    out = dict(
        status="computed at leading order (see module docstring); T,n in units of v, v^3",
        Pi_phi_total_over_T2=PI_PHI_TOTAL, Pi_X_over_T2=PI_X, Pi_phi_without_X_loop=PI_P,
        Tc_zero_charge_over_v=Tc0, Tc_leading_T2_only_over_v=float(np.sqrt(LAM_PHI / PI_PHI_TOTAL)),
        near_Tc_scan=near, quadrature_check=qc, rows=rows,
        local_line_dense=dict(T_over_v=Td, n_local_over_v3=nloc_d),
        GeV_conversion="T_GeV = T*v, n_GeV3 = n*v^3, v=%g GeV" % V_GEV,
    )
    save_json("thermal_phase.json", out)
    print("Tc0", Tc0, "quad", qc)
