"""Shared model for track gp_matching (Claude, independent code).

Fields are written as canonical real components: X=(x1+i x2)/sqrt2, Y=(y1+i y2)/sqrt2,
S=(s1+i s2)/sqrt2 (S = Sigma). Kinetic term |dX|^2+|dY|^2+|dS|^2, so the Hessian of V in
(x1,x2,y1,y2,s1,s2) is directly the mass-squared matrix. Units: H_inf = 1 unless stated.

Potential (7737487 CONSTRUCTION.md, with the review's supplied Hubble mass on X only, plus optional
tagged sensitivity terms c_Y, c_S that are NOT in 7737487):
  V = lam_X|X|^4 - c_H H^2 |X|^2 + lam_Y|Y|^4 + (M^2 + c_Y H^2)|Y|^2 + c_S H^2 |S|^2
      + lam_S(|S|^2 - h^2/2)^2 + [kap conj(Y) X^3 + mu conj(S) Y^2 + c.c.]
phi sits at 0 during the high-field epoch (mass^2 = 2 g_X^2|X|^2 >> H^2) and is omitted.
"""
import numpy as np
import sympy as sp

GEV = dict(M=1e12, h=1e11, mu=5e11, Hinf=4.7e13)
COUP = dict(lam_X=1e-3, lam_Y=1.0, lam_S=1.0, kap=1e-3)


def build(c_H=1.0, c_Y=0.0, c_S=0.0, Hinf=GEV['Hinf'], coup=None, extra_quartic=None):
    """Return dict of lambdified V, grad, hess in canonical real components (H units)."""
    c = dict(COUP) if coup is None else dict(coup)
    M, h, mu = GEV['M'] / Hinf, GEV['h'] / Hinf, GEV['mu'] / Hinf
    q = sp.symbols('x1 x2 y1 y2 s1 s2', real=True)
    X = (q[0] + sp.I * q[1]) / sp.sqrt(2)
    Y = (q[2] + sp.I * q[3]) / sp.sqrt(2)
    S = (q[4] + sp.I * q[5]) / sp.sqrt(2)
    aX = sp.expand(X * sp.conjugate(X))
    aY = sp.expand(Y * sp.conjugate(Y))
    aS = sp.expand(S * sp.conjugate(S))
    inter = c['kap'] * sp.conjugate(Y) * X**3 + mu * sp.conjugate(S) * Y**2
    V = (c['lam_X'] * aX**2 - c_H * aX + c['lam_Y'] * aY**2 + (M**2 + c_Y) * aY + c_S * aS
         + c['lam_S'] * (aS - h**2 / 2)**2 + 2 * sp.re(sp.expand(inter)))
    if extra_quartic:  # dict of portal couplings, e.g. {'XY': l} -> l |X|^2|Y|^2
        amap = dict(X=aX, Y=aY, S=aS)
        for k, val in extra_quartic.items():
            V += val * amap[k[0]] * amap[k[1]]
    V = sp.expand(V)
    g = [sp.diff(V, v) for v in q]
    Hs = [[sp.diff(gi, v) for v in q] for gi in g]
    f = sp.lambdify([q], V, 'numpy')
    fg = sp.lambdify([q], g, 'numpy')
    fh = sp.lambdify([q], Hs, 'numpy')
    return dict(V=lambda z: float(f(z)), grad=lambda z: np.array(fg(z), float),
                hess=lambda z: np.array(fh(z), float), M=M, h=h, mu=mu, coup=c,
                c_H=c_H, c_Y=c_Y, c_S=c_S, Hinf=Hinf, sym=(q, V))


def radial_V(R, m, c_H=None):
    """Potential at aligned phases with moduli R=(|X|,|Y|,|S|) (H units)."""
    c = m['coup']; cH = m['c_H'] if c_H is None else c_H
    x, y, s = R
    return (c['lam_X'] * x**4 - cH * x**2 + c['lam_Y'] * y**4 + (m['M']**2 + m['c_Y']) * y**2
            + m['c_S'] * s**2 + c['lam_S'] * (s**2 - m['h']**2 / 2)**2
            - 2 * c['kap'] * y * x**3 - 2 * m['mu'] * s * y**2)
