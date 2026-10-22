"""Track `dissipation`: right-hand sides of the three Klinkhamer-Volovik dissipative q-models and the
QV1-with-dissipation transcription (QV2-diss), plus a common integrator wrapper.

Models (dimensionless; units are those of each paper, stated per model):

KV16a  arXiv:1601.00601 (constant G).  Units E_* = sqrt(3) M (see s1_symbolic.py, check B), r = rho/E_*^4,
       h = H/E_*, tau = E_* t.  eps(f) = lambda + f^2 + 1/f^2 (eq 6.2b), r_V = lambda - f^2 + 3/f^2.
       Paper's system (6.1): f fdot eps'' = gamma |h| k^2,  k = hdot + 2 h^2 = 2 r_V  (radiation only).
       Our extension adds dust r_d (w=0): Friedmann h^2 = r_V + r_r + r_d,  hdot = -(3/2) r_d - 2 r_r,
       k = 2 r_V + r_d/2, r_V-dot = -gamma|h| k^2, r_r receives +gamma|h| k^2.  State: [h, f, r_d, ln a];
       r_r = h^2 - r_V - r_d (Friedmann used as definition; its consistency is exact, s1 check B/E).

KSV16  arXiv:1601.04676 (G^-1 = (8 pi/3) f q0, sigma = 1 units of 0806.2805: rho = r/chi0, H = h/sqrt(chi0 q0)).
       eps(f) = lambda + (1/2)(-f^2 + f^4/3)  (eq 18c; lambda = 0 in the paper, added here for the basin map).
       Regular first-order form (our rewriting; equivalent to (17a-c) on the constraint, s1 check C):
         fdd  = -(2 hdot + 3 h^2) f - 3 w r_M + 3 (eps - u f) - 2 h fdot      (space-space eq.)
         hdot = eps'(f) - u - 2 h^2                                          (q-equation, u = u_eff)
         udot = s = gamma_f fdot^2 + gamma_h hdot^2                           (eqs 18b, 21)
         r_M' = -3 h (1+w) r_M + f s                                          (eq 17c)
       Constraint C = f h^2 + h fdot - (eps - u f) - r_M obeys Cdot = -3 h C (s1 check C).
       Optional source 'pp': f s = gamma |h| (hdot + 2h^2)^2  (the 1601.00601 drain with nonconstant G;
       our transcription, KV do not print it).  State: [f, fdot, h, u, r_M, ln a].

K22    arXiv:2207.03453 (brane; E_P units, G = 1).  r_V(xi) = lambda + xi^2 - mu xi (eqs 13, 15e, 16a).
         r_M' + 3 h (1+w) r_M = Gamma,   xi-dot r_V'(xi) = -Gamma,   2 hdot + 3 h^2 = -8 pi (w r_M - r_V),
         3 h^2 = 8 pi (r_M + r_V) (constraint).  Gamma = g~ |h| h^4 (eq 17, 'adhoc') or
         gamma |h| (hdot + 2h^2)^2 (eq 8, 'pp').  State: [h, xi, r_M, ln a].
"""
import numpy as np
from scipy.integrate import solve_ivp


# ----------------------------------------------------------------------------------------------
# KV16a
# ----------------------------------------------------------------------------------------------
def kv16a_rV(f, lam):
    return lam - f**2 + 3.0 / f**2


def kv16a_epspp(f):
    return 2.0 + 6.0 / f**4


def kv16a_rhs(tau, y, lam, gam):
    h, f, rd, lna = y
    rV = kv16a_rV(f, lam)
    k = 2.0 * rV + 0.5 * rd
    rr = h * h - rV - rd
    hdot = -1.5 * rd - 2.0 * rr
    drain = gam * abs(h) * k * k          # = -rV-dot
    fdot = drain / (f * kv16a_epspp(f))   # since rV-dot = -f eps'' fdot
    return [hdot, fdot, -3.0 * h * rd, h]


def kv16a_f_of_rV(rV, lam):
    """invert r_V(f) = lam - x + 3/x, x = f^2 > 0 (monotone decreasing in x)."""
    x = 0.5 * ((lam - rV) + np.sqrt((lam - rV)**2 + 12.0))
    return np.sqrt(x)


def kv16a_fbar(lam):
    return kv16a_f_of_rV(0.0, lam)


# ----------------------------------------------------------------------------------------------
# KSV16 (and the nonconstant-G 'pp' variant)
# ----------------------------------------------------------------------------------------------
def ksv_eps(f, lam=0.0):
    return lam + 0.5 * (-f**2 + f**4 / 3.0)


def ksv_epsp(f):
    return -f + (2.0 / 3.0) * f**3


def ksv_epspp(f):
    return -1.0 + 2.0 * f**2


def ksv_rhs(tau, y, gf, gh, w, lam, source='ksv', gam=0.0):
    f, fd, h, u, rM, lna = y
    eps = ksv_eps(f, lam)
    hdot = ksv_epsp(f) - u - 2.0 * h * h
    fdd = -(2.0 * hdot + 3.0 * h * h) * f - 3.0 * w * rM + 3.0 * (eps - u * f) - 2.0 * h * fd
    if source == 'ksv':
        s = gf * fd * fd + gh * hdot * hdot
    elif source == 'pp':
        k = hdot + 2.0 * h * h
        s = gam * abs(h) * k * k / f
    else:
        s = 0.0
    rMd = -3.0 * h * (1.0 + w) * rM + f * s
    return [fd, fdd, hdot, s, rMd, h]


def ksv_initial(h, hdot, f, rM, lam=0.0):
    """KSV boundary data {h, hdot, f, r_M} -> state; u_eff from (19), fdot from constraint (17b)."""
    u = ksv_epsp(f) - hdot - 2.0 * h * h
    fd = (ksv_eps(f, lam) - u * f + rM - f * h * h) / h
    return np.array([f, fd, h, u, rM, 0.0])


def ksv_constraint(y, lam=0.0):
    f, fd, h, u, rM, lna = y
    return f * h * h + h * fd - (ksv_eps(f, lam) - u * f) - rM


def ksv_rVtilde(y, lam=0.0):
    f, fd, h, u, rM, lna = y
    return ksv_eps(f, lam) - u * f


def ksv_u0(lam=0.0):
    f0 = np.sqrt(0.5 * (1.0 + np.sqrt(1.0 + 8.0 * lam)))
    return ksv_epsp(f0), f0


# ----------------------------------------------------------------------------------------------
# K22 brane
# ----------------------------------------------------------------------------------------------
def k22_rV(xi, lam, mu):
    return lam + xi * xi - mu * xi


def k22_rhs(tau, y, lam, mu, w, gt, source='adhoc'):
    h, xi, rM, lna = y
    rV = k22_rV(xi, lam, mu)
    hdot = -0.5 * (3.0 * h * h + 8.0 * np.pi * (w * rM - rV))
    if source == 'adhoc':
        G = gt * abs(h) * h**4
    else:
        G = gt * abs(h) * (hdot + 2.0 * h * h)**2
    dr = 2.0 * xi - mu
    xid = -G / dr
    return [hdot, xid, -3.0 * h * (1.0 + w) * rM + G, h]


def k22_initial(h, xi, lam, mu):
    rM = 3.0 * h * h / (8.0 * np.pi) - k22_rV(xi, lam, mu)
    return np.array([h, xi, rM, 0.0])


# ----------------------------------------------------------------------------------------------
def integrate(rhs, t0, t1, y0, args, rtol=1e-10, atol=1e-13, events=None, t_eval=None, method='DOP853',
              dense=False, max_step=np.inf):
    return solve_ivp(rhs, (t0, t1), y0, args=args, method=method, rtol=rtol, atol=atol, events=events,
                     t_eval=t_eval, dense_output=dense, max_step=max_step)
