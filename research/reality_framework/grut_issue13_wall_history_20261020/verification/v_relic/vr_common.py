"""v_relic independent verifier: shared constants and hand-derived Feynman rules.

Written from scratch for the adversarial check; no code imported from modules/relic_yield.
Units: GeV. Potential (BRIEF.md):
  U = lam/4 (phi^2-v^2)^2 + g^2 phi^2 |X|^2 + lamX |X|^4 + 1/2 (mD0^2 + hD^2 phi^2) D^2 + lamD/4 D^4
Broken phase phi = v + h, X = (q1 + i q2)/sqrt2.
"""
import numpy as np

v = 1000.0
gX = 0.1
lam = 0.01
lamX = 0.001
hD = 0.01
lamD = 0.001
mD0 = 100.0
MP = 2.435323e18  # reduced Planck mass, GeV

mX = gX * v
mh = np.sqrt(2.0 * lam) * v
mD = np.sqrt(mD0**2 + hD**2 * v**2)

# --- Feynman vertices (i*L_int convention: vertex = -i * coefficient * symmetry factor) ---
# returned as the real number V such that vertex factor = -i V
V_hhh = 6.0 * lam * v          # lam v h^3  -> 3! lam v
V_hXX = 2.0 * gX**2 * v        # 2 g^2 v h X X*  (complex X: one X, one X* leg, no extra factor)
V_hhXX = 2.0 * gX**2           # g^2 h^2 X X* -> 2! g^2
V_hDD = 2.0 * hD**2 * v        # hD^2 v h D^2 -> 2! hD^2 v
V_hhDD = 2.0 * hD**2           # 1/2 hD^2 h^2 D^2 -> 2!2!/2 hD^2


def U(phi, X2, D):
    """potential with X2 = |X|^2"""
    return (lam / 4) * (phi**2 - v**2)**2 + gX**2 * phi**2 * X2 + lamX * X2**2 \
        + 0.5 * (mD0**2 + hD**2 * phi**2) * D**2 + (lamD / 4) * D**4


def check_vertices_fd():
    """finite-difference check of vertex coefficients from U (independent of the hand expansion).
    Use real basis X = q1/sqrt2 (q2=0): |X|^2 = q1^2/2. Coefficient of h q1 q1/2 ... we compare
    d^3U/dh dq1 dq1 = V_hXX (for real component q1, the h q1 q1 vertex equals 2g^2 v which coincides
    with the complex-field h X X* rule)."""
    hstep = 1.0
    def Ur(h, q1, d):
        return U(v + h, q1**2 / 2.0, d)
    # third derivatives by central differences of a polynomial (exact up to roundoff for step 1)
    def d3(f, i, j, k):
        e = np.eye(3) * hstep
        tot = 0.0
        for si in (+1, -1):
            for sj in (+1, -1):
                for sk in (+1, -1):
                    x = si * e[i] + sj * e[j] + sk * e[k]
                    tot += si * sj * sk * f(*x)
        return tot / (8 * hstep**3)
    def d4(f, i, j, k, l):
        e = np.eye(3) * hstep
        tot = 0.0
        for si in (+1, -1):
            for sj in (+1, -1):
                for sk in (+1, -1):
                    for sl in (+1, -1):
                        x = si * e[i] + sj * e[j] + sk * e[k] + sl * e[l]
                        tot += si * sj * sk * sl * f(*x)
        return tot / (16 * hstep**4)
    out = {
        "d3U_hhh": d3(Ur, 0, 0, 0), "V_hhh": V_hhh,
        "d3U_hq1q1": d3(Ur, 0, 1, 1), "V_hXX": V_hXX,
        "d4U_hhq1q1": d4(Ur, 0, 0, 1, 1), "V_hhXX": V_hhXX,
        "d3U_hDD": d3(Ur, 0, 2, 2), "V_hDD": V_hDD,
        "d4U_hhDD": d4(Ur, 0, 0, 2, 2), "V_hhDD": V_hhDD,
    }
    # second derivatives -> masses
    def d2(f, i, j):
        e = np.eye(3) * hstep
        tot = 0.0
        for si in (+1, -1):
            for sj in (+1, -1):
                tot += si * sj * f(*(si * e[i] + sj * e[j]))
        return tot / (4 * hstep**2)
    out["m2_h_fd"] = d2(Ur, 0, 0); out["m2_h"] = mh**2
    out["m2_q1_fd"] = d2(Ur, 1, 1); out["m2_X"] = mX**2
    out["m2_D_fd"] = d2(Ur, 2, 2); out["m2_D"] = mD**2
    return out


def M_hh_XX(s, cth, g2=None, lm=None):
    """tree amplitude hh -> X Xbar (sign convention: iM = sum of diagrams, vertex -iV, prop i/(p^2-m^2))
    mostly-minus metric. Optional overrides of g_X^2 and lambda_phi with v fixed (masses follow)."""
    g2 = gX**2 if g2 is None else g2
    lm = lam if lm is None else lm
    mh2 = 2 * lm * v**2
    mX2 = g2 * v**2
    Vhhh = 6 * lm * v
    VhXX = 2 * g2 * v
    VhhXX = 2 * g2
    p = np.sqrt(np.maximum(s / 4 - mh2, 0.0))
    k = np.sqrt(np.maximum(s / 4 - mX2, 0.0))
    t = mh2 + mX2 - 2 * (s / 4 - p * k * cth)
    u = mh2 + mX2 - 2 * (s / 4 + p * k * cth)
    # iM = (-iV)  ->  M = -V ; exchange: (-iV1)(-iV2) i/(q^2-m^2) = -i V1 V2/(q^2-m^2) -> M = -V1V2/(q^2-m^2)
    M = -VhhXX - Vhhh * VhXX / (s - mh2) - VhXX**2 / (t - mX2) - VhXX**2 / (u - mX2)
    return M


def M_hh_DD(s, cth):
    p = np.sqrt(np.maximum(s / 4 - mh**2, 0.0))
    k = np.sqrt(np.maximum(s / 4 - mD**2, 0.0))
    t = mh**2 + mD**2 - 2 * (s / 4 - p * k * cth)
    u = mh**2 + mD**2 - 2 * (s / 4 + p * k * cth)
    return -V_hhDD - V_hhh * V_hDD / (s - mh**2) - V_hDD**2 / (t - mD**2) - V_hDD**2 / (u - mD**2)


def M_DD_XX(s, cth):
    return -V_hDD * V_hXX / (s - mh**2) + 0.0 * cth


_xg, _wg = np.polynomial.legendre.leggauss(64)


def sigma(s, Mfun, mi, mf, Sf):
    """2->2 cross section for equal-mass initial pair mi, final pair mf; Sf final symmetry factor."""
    s = np.atleast_1d(s)
    out = np.zeros_like(s)
    ok = (s > 4 * mi**2) & (s > 4 * mf**2)
    ss = s[ok]
    p = np.sqrt(ss / 4 - mi**2)
    k = np.sqrt(ss / 4 - mf**2)
    M2 = np.array([np.sum(_wg * Mfun(si, _xg)**2) for si in ss])
    out[ok] = Sf * k / (32 * np.pi * ss * p) * M2
    return out
