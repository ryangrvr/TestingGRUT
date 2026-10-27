#!/usr/bin/env python3
"""cf_scale common RG module (own code): lambdified one-loop betas from beta_derive.py (sympy), the unit-sphere
minimum of V4, the gravity terms and the scale constants.

Couplings vector order: NAMES (11 quartics) + ['e'].
Gravity above M_P (sharp switch, as in de Brito et al. arXiv:2312.02086 eq. 31):
   beta_lambda_i += A_lam * lambda_i   (A_lam > 0: quartic irrelevant; A_lam = -f_lambda in Eichhorn's convention)
   beta_e        += -f_g * e          (f_g > 0: e = 0 relevant, e* = sqrt(16 pi^2 f_g / b_e) irrelevant)
"""
import sys
from pathlib import Path
import numpy as np
import sympy as sp

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import beta_derive as BD

NAMES = BD.NAMES
ALL = NAMES + ['e']
NC = len(NAMES)
MP_RED = 2.435e18          # GeV, reduced Planck mass (gravity switch scale)
B_E = 46.0 / 3.0           # 16 pi^2 beta_e = B_E e^3 (charges 1, 3, 6; complex scalars)
K16 = 16 * np.pi**2

_V4 = sum(BD.C[n] * BD.MONO[n] for n in NAMES)
_B = BD.beta_V(_V4)
_BETAS, _res = BD.project(_B)
assert _res == 0
_syms = [BD.C[n] for n in NAMES] + [BD.e]
_beta_f = sp.lambdify(_syms, [_BETAS[n] for n in NAMES], 'numpy')
_V4f = sp.lambdify([BD.F, _syms[:-1]], _V4, 'numpy')
_V4grad = sp.lambdify([BD.F, _syms[:-1]], [sp.diff(_V4, v) for v in BD.F], 'numpy')


def beta(y, A_lam=0.0, f_g=0.0, b_e=B_E, gauge_active=True):
    lam = y[:NC]; e = y[NC] if gauge_active else 0.0
    b = np.array(_beta_f(*lam, e), float) / K16 + A_lam * lam
    be = (b_e * y[NC]**3 / K16 - f_g * y[NC]) if gauge_active else 0.0
    return np.append(b, be)


def V4_unit_min(lam, n=60, seed=3, fields='all'):
    """min of V4 on the unit sphere in the 7 real field components (or a subset); returns (min, argmin)."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    mask = np.ones(7)
    if fields == 'hidden_charged':
        mask = np.array([1, 1, 1, 1, 1, 1, 0.])
    lam = list(lam)

    def f(v):
        w = v * mask; nn = np.dot(w, w)
        return float(_V4f(w, lam)) / nn**2

    best = None
    for _ in range(n):
        r = minimize(f, rng.normal(size=7), method='BFGS')
        if best is None or r.fun < best.fun:
            best = r
    w = best.x * mask
    return float(best.fun), (w / np.linalg.norm(w)).tolist()


def vec(d):
    return np.array([d.get(n, 0.0) for n in ALL], float)


def asdict(y):
    return {n: float(v) for n, v in zip(ALL, y)}
