#!/usr/bin/env python3
"""Item 1: symbolic re-derivation of the constant-coefficient CFS dispersion polynomial
and its Hurwitz determinants (own derivation, sympy), plus the initial layer charge fraction."""
import json
from pathlib import Path
import numpy as np
import sympy as sp
OUT = Path(__file__).resolve().parents[1] / 'results'
s, k, m, al, sg = sp.symbols('s k m alpha sigma', positive=True)
a = al + sg
Gam = 1 + sg / (s + al)                         # stretch factor
# plane wave exp(i k x): s^2 + m^2 + k^2 / Gam^2 = 0, i.e. (s^2+m^2) Gam^2 + k^2 = 0, times (s+alpha)^2
poly = sp.expand(sp.cancel(((s ** 2 + m ** 2) * Gam ** 2 + k ** 2) * (s + al) ** 2))
P = sp.Poly(poly, s)
c = P.all_coeffs()
claimed = [1, 2 * a, a ** 2 + m ** 2 + k ** 2, 2 * (a * m ** 2 + al * k ** 2), a ** 2 * m ** 2 + al ** 2 * k ** 2]
coeff_ok = all(sp.simplify(x - y) == 0 for x, y in zip(c, claimed))
a0, a1, a2, a3, a4 = c
D2 = sp.factor(a1 * a2 - a0 * a3)
D3 = sp.factor(a3 * (a1 * a2 - a0 * a3) - a1 ** 2 * a4)
D2_ok = sp.simplify(D2 - 2 * (a ** 3 + sg * k ** 2)) == 0
D3_ok = sp.simplify(D3 - 4 * sg * k ** 2 * (a * m ** 2 + al * a ** 2 + al * k ** 2)) == 0
# time-domain elimination: psi_hat = -sf/(s+alpha+sf) G U, q_hat = -sn/(s+alpha+sn) W  =>  C = (s+alpha)/(s+alpha+sigma)
C = sp.simplify(1 - sg / (s + al + sg))
C_ok = sp.simplify(C - 1 / Gam) == 0
# initial layer charge fraction of the prepared Gaussian (h = 0.25, R = 360, L = 450)
h = .25; r = np.arange(1, 1800) * h; f = .18 * np.exp(-r ** 2 / 7200.)
q = (r * f) ** 2
frac = float(np.sum(q[r >= 360]) / np.sum(q))
out = dict(dispersion_coefficients_match=bool(coeff_ok), Delta2=str(D2), Delta2_matches=bool(D2_ok),
           Delta3=str(D3), Delta3_matches=bool(D3_ok), memory_elimination_gives_C_equals_1_over_Gamma=bool(C_ok),
           initial_charge_fraction_in_layer_R360=frac)
(OUT / 'routh_and_memory.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
