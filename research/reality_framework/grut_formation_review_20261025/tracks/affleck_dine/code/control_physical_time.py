#!/usr/bin/env python3
"""Control: integrate the AD equation in physical time and compute n/s directly (own code).

Checks the conformal reduction and the yield formula of window.py for one radiation-era
release: Xddot + 3H Xdot + 2 lX |X|^2 X + 6 lambda6^* X*^5 = 0, H = 1/(2t),
s = (2 pi^2/45) g* T^3 with H = sqrt(pi^2 g*/90) T^2/M.
Units: M = 1 (reduced Planck); H_rel = 1e-5.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ad_condensate import run

lX, gs, M = 1e-3, 106.75, 1.0
H_rel, cH, eta, phi0 = 1e-5, 1.0, 0.1, np.pi / 2
X0 = np.sqrt(cH / (2 * lX)) * H_rel
lam6 = eta * lX / X0 ** 2  # |lambda_6| = c_6/M^2
t_r = 1 / (2 * H_rel)
w0 = np.sqrt(2 * lX) * X0


def rhs(t, y):
    X = y[0] + 1j * y[1]; V = y[2] + 1j * y[3]
    acc = -3 / (2 * t) * V - 2 * lX * abs(X) ** 2 * X - 6 * lam6 * np.conj(X) ** 5
    return [V.real, V.imag, acc.real, acc.imag]


X_init = X0 * np.exp(1j * phi0 / 6)
t_end = t_r * 1e6  # a = (t/t_r)^(1/2) = 1e3
sol = solve_ivp(rhs, (t_r, t_end), [X_init.real, X_init.imag, 0, 0], method='DOP853', rtol=1e-10, atol=1e-14 * X0,
                dense_output=False)
X = sol.y[0] + 1j * sol.y[1]; V = sol.y[2] + 1j * sol.y[3]
a = np.sqrt(sol.t / t_r)
na3 = -2 * np.imag(np.conj(X) * V) * a ** 3
n_tail = float(np.mean(na3[int(0.8 * len(na3)):]))
T_rel = (90 / (np.pi ** 2 * gs)) ** 0.25 * np.sqrt(H_rel * M)
s_rel = (2 * np.pi ** 2 / 45) * gs * T_rel ** 3
Y_direct = abs(n_tail) / s_rel
D = (2 * np.pi ** 2 / 45) * gs * (90 / (np.pi ** 2 * gs)) ** 0.75
r = run(6, eta, phi0, 1 / np.sqrt(cH), 0.0, 'RD')
Y_formula = abs(r['q_inf']) * (H_rel / M) ** 1.5 / (2 * lX * (1 / np.sqrt(cH)) ** 3 * D)
out = dict(Y_direct=Y_direct, Y_formula=Y_formula, rel_diff=abs(Y_direct - Y_formula) / Y_formula,
           q_conformal=r['q_inf'], na3_over_X0sq_w0=n_tail / (X0 ** 2 * w0))
print(out)
(Path(__file__).resolve().parents[1] / 'results' / 'control_physical_time.json').write_text(json.dumps(out, indent=1))
