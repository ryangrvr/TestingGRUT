#!/usr/bin/env python3
"""kin_rates (Claude, own code): the light-sector model at the heavy-mediator spot point.

Fields: complex X = (x1 + i x2)/sqrt2 and the real radial field phi = v + s.  The D field of the gauge-point packet
is absent from the spot point (parent INPUTS.json light_sector lists only gX, lambdaX, lambdaPhi, vPhi) and is not
included.  Potential [Q, parent INPUTS.json / amp42.py convention]:

    U = lP/4 (phi^2 - v^2)^2 + g^2 phi^2 |X|^2 + lX |X|^4 .

Tree masses: m_X = g v, m_s^2 = 2 lP v^2.  Vertex factors are the exact derivatives of U at the vacuum, in the real
basis (x1, x2, s); the amplitude convention is M = -prod(V) prod 1/(p^2 - m^2) for every tree diagram, so that the
2->2 contact amplitude is -V4.  Useful couplings:
    lam3  = d^3U/ds dx1 dx1 = 2 g^2 v           (s X Xbar trilinear, 0.018482 GeV)  [D]
    lam_s3 = d^3U/ds^3     = 6 lP v             (s^3, 0.16634 GeV)                  [D]
    lam_s4 = d^4U/ds^4     = 6 lP               (s^4)                               [D]
    lam_ssXX = d^4U/ds ds dx1 dx1 = 2 g^2       (s s X Xbar quartic)                [D]
    lam_X4 : d^4U/dx1^4 = 6 lX, d^4U/dx1^2 dx2^2 = 2 lX  (X X Xbar Xbar contact amplitude -4 lX)
Vertex tensors V3[a,b,c], V4[a,b,c,d] are built once (sympy) and cached as numpy arrays; they can be rebuilt for
in-medium values (v_T, couplings unchanged) with `vertices(v)`.
"""
import itertools
import numpy as np
import sympy as sp

# ---- inputs [Q] (parent INPUTS.json at 225028b, light_sector) ----
GX = 0.1
LX = 1e-3
LP = 0.03
V = 0.924089637375362
MX = GX * V                      # 0.0924089637 GeV
MS = np.sqrt(2 * LP) * V         # 0.2263548088 GeV
LAM3 = 2 * GX**2 * V             # 0.018481793 GeV
LAM_S3 = 6 * LP * V              # 0.166336 GeV
LAM_S4 = 6 * LP
LAM_SSXX = 2 * GX**2
GAMMA_S_VAC = LAM3**2 * np.sqrt(1 - 4 * MX**2 / MS**2) / (16 * np.pi * MS)   # 1.7333e-5 GeV [D]

NFIELD = 3          # x1, x2, s
FIELD_NAMES = ('x1', 'x2', 's')
# external "polarisations" in the real basis: incoming X, incoming Xbar, s.  An outgoing X is an incoming Xbar
# with momentum -k (and vice versa); an outgoing s is an incoming s with -k.
E_X_IN = np.array([1, -1j, 0]) / np.sqrt(2)
E_XB_IN = np.array([1, 1j, 0]) / np.sqrt(2)
E_S = np.array([0, 0, 1], dtype=complex)
CHARGE = {'X': 1, 'Xb': -1, 's': 0}
POL_IN = {'X': E_X_IN, 'Xb': E_XB_IN, 's': E_S}
POL_OUT = {'X': E_XB_IN, 'Xb': E_X_IN, 's': E_S}


def vertices(v=V, g=GX, lX=LX, lP=LP):
    """Mass^2 vector and V3, V4 tensors (exact derivatives of U at x=0, phi=v)."""
    x1, x2, s = sp.symbols('x1 x2 s', real=True)
    F = [x1, x2, s]
    phi = v + s
    X2 = (x1**2 + x2**2) / 2
    U = sp.Rational(1, 4) * lP * (phi**2 - v**2)**2 + g**2 * phi**2 * X2 + lX * X2**2
    z = {f: 0 for f in F}
    m2 = np.array([float(sp.diff(U, a, a).subs(z)) for a in F])
    V3 = np.zeros((3, 3, 3))
    V4 = np.zeros((3, 3, 3, 3))
    for i, j, k in itertools.product(range(3), repeat=3):
        V3[i, j, k] = float(sp.diff(U, F[i], F[j], F[k]).subs(z))
    for i, j, k, l in itertools.product(range(3), repeat=4):
        V4[i, j, k, l] = float(sp.diff(U, F[i], F[j], F[k], F[l]).subs(z))
    return m2, V3, V4


_CACHE = {}


def get_vertices(v=V):
    key = round(float(v), 15)
    if key not in _CACHE:
        _CACHE[key] = vertices(v)
    return _CACHE[key]


def masses(v=V):
    return GX * v, np.sqrt(2 * LP) * v


if __name__ == '__main__':
    m2, V3, V4 = get_vertices()
    print('m2 from U:', m2, ' expected', MX**2, MX**2, MS**2)
    print('lam3 (s x1 x1):', V3[2, 0, 0], LAM3, ' s^3:', V3[2, 2, 2], LAM_S3)
    print('s s x1 x1:', V4[2, 2, 0, 0], LAM_SSXX, ' s^4:', V4[2, 2, 2, 2], LAM_S4,
          ' x1^4:', V4[0, 0, 0, 0], 6 * LX, ' x1x1x2x2:', V4[0, 0, 1, 1], 2 * LX)
    print('Gamma_s vacuum:', GAMMA_S_VAC)
