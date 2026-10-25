#!/usr/bin/env python3
"""Tree amplitudes of the light sector by Berends-Giele recursion (own code, track gp_leftovers).

Real fields (x1, x2, s, D) with X = (x1 + i x2)/sqrt2, phi = v + s and
U = lp/4 (phi^2-v^2)^2 + g^2 phi^2 |X|^2 + lX |X|^4 + (mD0^2 + hD^2 phi^2) D^2/2 + lD D^4/4.
Vertices are exact derivatives of U (sympy). Every tree diagram has M = -prod(V) prod(1/(p^2-m^2)).
External polarisations: incoming X (1,-i)/sqrt2, incoming Xbar (1,+i)/sqrt2; an outgoing X is an incoming
Xbar with momentum -k. Outputs threshold amplitudes and the NR rate coefficients
A = |M|^2 Phi_2 / (prod 2m_i) / S_final  for 4 -> 2 processes and checks against analytic 2 -> 2.
"""
import json
import itertools
from pathlib import Path
import numpy as np
import sympy as sp
from common import G, LX, LP, V, HD, MX, MS, MD

OUT = Path(__file__).resolve().parents[1] / 'results'
LD = 1e-3
x1, x2, s, D = sp.symbols('x1 x2 s D', real=True)
F = [x1, x2, s, D]
phi = V + s
X2 = (x1**2 + x2**2) / 2
U = LP / 4 * (phi**2 - V**2)**2 + G**2 * phi**2 * X2 + LX * X2**2 + ((0.1 * V)**2 + HD**2 * phi**2) * D**2 / 2 + LD * D**4 / 4
z = {f: 0 for f in F}
MASS2 = np.array([float(sp.diff(U, a, a).subs(z)) for a in F])
V3 = np.zeros((4, 4, 4)); V4 = np.zeros((4, 4, 4, 4))
for i, j, k in itertools.product(range(4), repeat=3):
    V3[i, j, k] = float(sp.diff(U, F[i], F[j], F[k]).subs(z))
for i, j, k, l in itertools.product(range(4), repeat=4):
    V4[i, j, k, l] = float(sp.diff(U, F[i], F[j], F[k], F[l]).subs(z))
EX_IN = np.array([1, -1j, 0, 0]) / np.sqrt(2)
EXB_IN = np.array([1, 1j, 0, 0]) / np.sqrt(2)
ES = np.array([0, 0, 1, 0], dtype=complex)
ED = np.array([0, 0, 0, 1], dtype=complex)
METRIC = np.array([1, -1, -1, -1])


def amplitude(eps, mom):
    """eps: list of field-space vectors, mom: list of 4-momenta, all incoming, sum = 0."""
    n = len(eps)
    J = {}
    for i in range(n - 1):
        J[1 << i] = np.array(eps[i], dtype=complex)
    full = (1 << (n - 1)) - 1

    def P(S):
        return sum(mom[i] for i in range(n - 1) if S >> i & 1)

    def kernel(S):
        K = np.zeros(4, dtype=complex)
        low = S & -S
        sub = (S - 1) & S
        while sub:
            if sub & low and sub != S:
                S2 = S ^ sub
                K += np.einsum('abc,b,c->a', V3, J[sub], J[S2])
                rest = S2
                low2 = rest & -rest
                sub2 = (rest - 1) & rest
                while sub2:
                    if sub2 & low2 and sub2 != rest:
                        K += np.einsum('abcd,b,c,d->a', V4, J[sub], J[sub2], J[rest ^ sub2])
                    sub2 = (sub2 - 1) & rest
            sub = (sub - 1) & S
        return K

    for size in range(2, n - 1):
        for S in range(1, full + 1):
            if bin(S).count('1') != size:
                continue
            p = P(S)
            p2 = float(np.dot(METRIC, p * p))
            J[S] = kernel(S) / (p2 - MASS2)
    return -np.dot(eps[n - 1], kernel(full))


def rest(m):
    return np.array([m, 0, 0, 0.])


def two_to_two_checks():
    # X(p1) Xbar(p2) -> X(k1) Xbar(k2), generic NR kinematics
    m = MX
    pin = 0.3 * m
    E = np.sqrt(m * m + pin * pin)
    p1 = np.array([E, 0, 0, pin]); p2 = np.array([E, 0, 0, -pin])
    th = 0.7
    k1 = np.array([E, pin * np.sin(th), 0, pin * np.cos(th)]); k2 = np.array([E, -pin * np.sin(th), 0, -pin * np.cos(th)])
    Mnum = amplitude([EX_IN, EXB_IN, EXB_IN, EX_IN], [p1, p2, -k1, -k2])
    sM = float(np.dot(METRIC, (p1 + p2)**2)); tM = float(np.dot(METRIC, (p1 - k1)**2)); uM = float(np.dot(METRIC, (p1 - k2)**2))
    lam3 = 2 * G**2 * V
    Mana = -4 * LX - lam3**2 * (1 / (sM - MS**2) + 1 / (tM - MS**2))
    # X X -> X X (t and u channels)
    Mxx = amplitude([EX_IN, EX_IN, EXB_IN, EXB_IN], [p1, p2, -k1, -k2])
    Mxx_ana = -4 * LX - lam3**2 * (1 / (tM - MS**2) + 1 / (uM - MS**2))
    # charge violating X X -> X Xbar must vanish
    Mcv = amplitude([EX_IN, EX_IN, EXB_IN, EX_IN], [p1, p2, -k1, -k2])
    # D D -> X Xbar via s exchange only
    mD = MD
    ED_ = np.sqrt(mD**2 + (0.1 * mD)**2)
    q1 = np.array([ED_, 0, 0, 0.1 * mD]); q2 = np.array([ED_, 0, 0, -0.1 * mD])
    kk = np.sqrt(ED_**2 - m * m)
    r1 = np.array([ED_, 0, 0, kk]); r2 = np.array([ED_, 0, 0, -kk])
    Mdd = amplitude([ED, ED, EXB_IN, EX_IN], [q1, q2, -r1, -r2])
    sD = 4 * ED_**2
    Mdd_ana = -(2 * HD**2 * V) * lam3 / (sD - MS**2)
    return dict(XXbar=dict(num=complex(Mnum).real, ana=Mana, imag=abs(complex(Mnum).imag)),
                XX=dict(num=complex(Mxx).real, ana=Mxx_ana),
                charge_violating_abs=abs(Mcv),
                DD_to_XXbar=dict(num=complex(Mdd).real, ana=Mdd_ana),
                threshold_XX_to_XX_over4=None)


def threshold_42():
    m = MX
    kf = np.sqrt(3) * m
    p = [rest(m)] * 4
    k1 = np.array([2 * m, 0, 0, kf]); k2 = np.array([2 * m, 0, 0, -kf])
    out = {}
    # a: X X Xbar Xbar -> X Xbar ; b: X X X Xbar -> X X ; cv: X X X X -> X Xbar
    Ma = amplitude([EX_IN, EX_IN, EXB_IN, EXB_IN, EXB_IN, EX_IN], p + [-k1, -k2])
    Mb = amplitude([EX_IN, EX_IN, EX_IN, EXB_IN, EXB_IN, EXB_IN], p + [-k1, -k2])
    Mcv = amplitude([EX_IN, EX_IN, EX_IN, EX_IN, EXB_IN, EX_IN], p + [-k1, -k2])
    # direction independence
    th = 1.1
    k1r = np.array([2 * m, kf * np.sin(th), 0, kf * np.cos(th)]); k2r = k1r * np.array([1, -1, -1, -1])
    Ma_rot = amplitude([EX_IN, EX_IN, EXB_IN, EXB_IN, EXB_IN, EX_IN], p + [-k1r, -k2r])
    # c: X X Xbar Xbar -> D D
    kD = np.sqrt(4 * m * m - MD**2)
    d1 = np.array([2 * m, 0, 0, kD]); d2 = np.array([2 * m, 0, 0, -kD])
    Mc = amplitude([EX_IN, EX_IN, EXB_IN, EXB_IN, ED, ED], p + [-d1, -d2])
    phi2 = lambda k: k / (4 * np.pi * 4 * m)
    flux = (2 * m)**4
    A = dict(a=abs(Ma)**2 * phi2(kf) / flux, b=abs(Mb)**2 * phi2(kf) / flux / 2, c=abs(Mc)**2 * phi2(kD) / flux / 2)
    out = dict(M_a=complex(Ma).real, M_a_imag=abs(complex(Ma).imag), M_a_rotated=complex(Ma_rot).real,
               M_b=complex(Mb).real, M_c=complex(Mc).real, M_charge_violating=abs(Mcv),
               A_GeV_minus8=A,
               lambdaX_only_estimate=None)
    return out


def coupling_scan():
    """How the 4->2 coefficient depends on lambda_X alone (s couplings off) - scale check."""
    global V3, V4
    saved = V3.copy(), V4.copy()
    V3[:] = 0; V4[:] = 0
    V4[0, 0, 0, 0] = V4[1, 1, 1, 1] = 6 * LX
    for (i, j) in ((0, 1),):
        for perm in set(itertools.permutations([i, i, j, j])):
            V4[perm] = 2 * LX
    m = MX; kf = np.sqrt(3) * m
    p = [rest(m)] * 4
    k1 = np.array([2 * m, 0, 0, kf]); k2 = np.array([2 * m, 0, 0, -kf])
    Ma = amplitude([EX_IN, EX_IN, EXB_IN, EXB_IN, EXB_IN, EX_IN], p + [-k1, -k2])
    V3[:], V4[:] = saved
    return complex(Ma).real


def main():
    res = dict(masses_GeV=dict(mX=MX, ms=MS, mD=MD), mass2_from_U=MASS2.tolist(),
               checks=two_to_two_checks(), threshold_4to2=threshold_42())
    res['threshold_4to2']['M_a_lambdaX_only'] = coupling_scan()
    # threshold elastic amplitudes
    lam3 = 2 * G**2 * V
    res['threshold_elastic'] = dict(XXbar=-4 * LX - lam3**2 * (1 / (4 * MX**2 - MS**2) - 1 / MS**2),
                                    XX=-4 * LX + 2 * lam3**2 / MS**2,
                                    attractive_quartic_check=(-4 * LX + 2 * lam3**2 / MS**2) / 4)
    OUT.mkdir(exist_ok=True)
    (OUT / 'amplitudes.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
