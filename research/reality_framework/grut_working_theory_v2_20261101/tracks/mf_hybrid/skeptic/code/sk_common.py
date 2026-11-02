"""mf_hybrid_skeptic common module (Claude, adversarial pass on tracks/mf_hybrid; GitHub issue #13 context).

Independent re-implementation of the HF potential from the parameters exported in bench_*.json (not imported from the
track's code). Units: H_ref = 1.30e11 GeV (the track's design H), reduced Planck mass mP = M_P/H_ref.
Every potential term, renormalisation condition and parameter is the track's model choice; nothing follows from closure.

  dV(sigma) = c1 x - m2 x^2/2 + n_chi M^4 (ln|M^2|/mu^2 - 3/2)/(64 pi^2) - c^2/(2 M^2),   x = sigma - sigma_c,
  M^2(sigma) = g^2 sigma^2 - lam v^2   (chi mass at chi = 0; the chi loop for sigma > sigma_c, one real d.o.f.)
  V = V0 + dV, V0 = lam v^4/4 = 3 mP^2 (the constant offset subtracted analytically to avoid cancellation).
Two-field tree: (lam/4)(chi^2 - v^2)^2 + g^2 sigma^2 chi^2/2 - c chi + c1 x - m2 x^2/2 (+ loop for sigma > sigma_c) + Vc.
"""
import json, math, os
import numpy as np

MP_GEV = 2.435e18
H_REF = 1.30e11
mP = MP_GEV / H_REF
AS = 2.100e-9
HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, '..', '..', 'results')
OUT = os.path.join(HERE, '..', 'results')


def load_bench(name):
    d = json.load(open(os.path.join(BENCH, 'bench_%s.json' % name)))
    pr = d['params']
    P = dict(lam=pr['lam'], v=pr['v_MP'] * mP, g=pr['g'], sc=pr['sigma_c_MP'] * mP, c1=pr['c1_Hi3'], m2=pr['m2_Hi2'],
             mu2=pr['mu_R2_Hi2'], cb=pr['c_bias_Hi3'], m_chi=pr['m_chi_over_Hi'], m_sig=pr['m_sig_over_Hi'],
             s_star=pr['s_star_MP'] * mP, g_s=pr['g_s'], smech=pr['smech'])
    P['V0'] = 0.25 * P['lam'] * P['v'] ** 4
    return d, P


def hist(d):
    return {k: np.array([np.nan if x is None else x for x in v], dtype=float) for k, v in d['history'].items()}


def M2(P, s):
    return P['g'] ** 2 * s * s - P['lam'] * P['v'] ** 2


def dV1(P, s, loop=True, bias=True):
    """V - V0 along the inflationary valley (single-field effective potential), sigma > sigma_c."""
    x = s - P['sc']
    out = P['c1'] * x - 0.5 * P['m2'] * x * x
    m2 = M2(P, s)
    if loop:
        out += m2 * m2 * (math.log(abs(m2) / P['mu2']) - 1.5) / (64 * math.pi ** 2)
    if bias:
        out -= 0.5 * P['cb'] ** 2 / m2
    return out


def dV1_d(P, s, n=1, loop=True, bias=True):
    """Derivatives by Richardson-extrapolated central differences (an independent route from the track's analytic ones)."""
    def cd(h):
        if n == 1:
            return (dV1(P, s + h, loop, bias) - dV1(P, s - h, loop, bias)) / (2 * h)
        if n == 2:
            return (dV1(P, s + h, loop, bias) - 2 * dV1(P, s, loop, bias) + dV1(P, s - h, loop, bias)) / h ** 2
        return (dV1(P, s + 2 * h, loop, bias) - 2 * dV1(P, s + h, loop, bias) + 2 * dV1(P, s - h, loop, bias)
                - dV1(P, s - 2 * h, loop, bias)) / (2 * h ** 3)
    h = 1e-3 * (s - P['sc']) if n == 1 else 3e-3 * (s - P['sc'])
    return (4 * cd(h / 2) - cd(h)) / 3


def V2(P, s, x, Vc=0.0, loop=True):
    d = s - P['sc']
    out = (0.25 * P['lam'] * (x * x - P['v'] ** 2) ** 2 + 0.5 * P['g'] ** 2 * s * s * x * x - P['cb'] * x
           + P['c1'] * d - 0.5 * P['m2'] * d * d + Vc)
    if loop and s > P['sc']:
        m2 = M2(P, s)
        out += m2 * m2 * (math.log(abs(m2) / P['mu2']) - 1.5) / (64 * math.pi ** 2)
    return out


def grad2(P, s, x, loop=True):
    d = s - P['sc']
    Vs = P['g'] ** 2 * s * x * x + P['c1'] - P['m2'] * d
    if loop and s > P['sc']:
        m2 = M2(P, s)
        Vs += m2 * 2 * P['g'] ** 2 * s * (math.log(abs(m2) / P['mu2']) - 1.0) / (32 * math.pi ** 2)
    Vx = P['lam'] * (x * x - P['v'] ** 2) * x + P['g'] ** 2 * s * s * x - P['cb']
    return Vs, Vx


def hess2(P, s, x):
    """tree Hessian (the loop term's sigma curvature is ~1e-2 H^2, irrelevant for the waterfall modes)."""
    Vss = P['g'] ** 2 * x * x - P['m2']
    Vxx = P['lam'] * (3 * x * x - P['v'] ** 2) + P['g'] ** 2 * s * s
    Vsx = 2 * P['g'] ** 2 * s * x
    return Vss, Vxx, Vsx


def vacuum(P):
    """(sigma0, chi0) of the true vacuum by Newton iteration on the tree gradient."""
    s, x = 0.0, P['v']
    for _ in range(60):
        Vs, Vx = grad2(P, s, x, loop=False)
        a, b, c = hess2(P, s, x)
        det = a * b - c * c
        ds = (b * Vs - c * Vx) / det
        dx = (a * Vx - c * Vs) / det
        s, x = s - ds, x - dx
    return s, x


def dump(name, obj):
    path = os.path.join(OUT, name)
    json.dump(obj, open(path, 'w'), indent=1, default=float)
    return path
