#!/usr/bin/env python3
"""Item 2: staggered fourth-order gradient/divergence. Own code.
(a) D = -G^T, L_7pt = D G = -G^T G with the odd/even closures (exact SBP);
(b) order of accuracy in the interior and at the closures;
(c) lattice dispersion symbol; (d) lead-root structure (2 oscillatory, 2+2 evanescent);
(e) semi-discrete energy identity.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from lattice import *

OUT = Path(__file__).resolve().parents[1] / 'results'
res = {}

# (a) exact summation by parts
for n, h in [(20, 1.0), (57, 0.3)]:
    G, D, L7 = G_matrix(n, h), D_matrix(n, h), L_matrix(n, h)
    res[f'sbp_n{n}'] = dict(
        max_abs_D_plus_GT=float(np.max(abs(D + G.T))),
        max_abs_L7_minus_DG=float(np.max(abs(L7 - D @ G)) * h * h),
        max_abs_L7_plus_GTG=float(np.max(abs(L7 + G.T @ G)) * h * h),
        L7_symmetric=float(np.max(abs(L7 - L7.T)) * h * h),
        max_eig_L7=float(np.max(np.linalg.eigvalsh(L7))))

# (b) order of accuracy
def errors(func, d1, d2, Lr, hs):
    rows = []
    for h in hs:
        n = int(round(Lr / h))
        rn = np.arange(1, n) * h
        rf = (np.arange(n) + .5) * h
        U = func(rn)
        G, L7 = G_matrix(n, h), L_matrix(n, h)
        eG = G @ U - d1(rf)
        eL = L7 @ U - d2(rn)
        near = rn < 4 * h + 1e-12
        rows.append(dict(h=h, grad_max=float(np.max(abs(eG))), lap_max=float(np.max(abs(eL))),
                         lap_interior_max=float(np.max(abs(eL[(rn > 1) & (rn < Lr - 1)]))),
                         lap_first3_nodes_max=float(np.max(abs(eL[near]))),
                         lap_last3_nodes_max=float(np.max(abs(eL[rn > Lr - 4 * h - 1e-12])))))
    for a, b in zip(rows, rows[1:]):
        for key in ('grad_max', 'lap_max', 'lap_interior_max', 'lap_first3_nodes_max', 'lap_last3_nodes_max'):
            b['order_' + key] = float(np.log(a[key] / b[key]) / np.log(a['h'] / b['h']))
    return rows

hs = [0.4, 0.2, 0.1, 0.05, 0.025]
Lr = 8.0
# regular field: X even in r -> U = r X odd; wall at 8 where it is negligible
f = lambda r: r * np.exp(-(r - 0.0) ** 2 / 2) * (1 + 0.3 * r ** 2)
# derivatives by sympy-free closed forms
def f1(r):
    return np.exp(-r * r / 2) * ((1 + .3 * r * r) + r * (.6 * r) - r * r * (1 + .3 * r * r))
def f2(r):
    g = 1 + .3 * r * r; gp = .6 * r; gpp = .6
    e = np.exp(-r * r / 2)
    # (r g e)'' with e' = -r e, e'' = (r^2-1) e
    return e * (2 * gp + r * gpp + 2 * (g + r * gp) * (-r) + r * g * (r * r - 1))
res['order_regular_odd_U'] = errors(f, f1, f2, Lr, hs)
# irregular field: U = r^2 e^{-r^2/2}-type (U''(0) != 0, i.e. X ~ r, cusp |r|)
g = lambda r: r * r * np.exp(-r * r / 2)
g1 = lambda r: np.exp(-r * r / 2) * (2 * r - r ** 3)
g2 = lambda r: np.exp(-r * r / 2) * (2 - 5 * r * r + r ** 4)
res['order_irregular_U_with_U2_at_0'] = errors(g, g1, g2, Lr, hs)
# Dirichlet wall with U_rr(L) != 0 (odd continuation not exact there)
Lw = 3 * np.pi
w = lambda r: np.sin(r) * (1 + 0.1 * r * r)          # odd at 0, U(3pi)=0, U''(3pi) = -0.6*pi*2 != 0
w1 = lambda r: np.cos(r) * (1 + 0.1 * r * r) + np.sin(r) * 0.2 * r
w2 = lambda r: -np.sin(r) * (1 + 0.1 * r * r) + 2 * np.cos(r) * 0.2 * r + np.sin(r) * 0.2
res['order_wall_with_nonzero_U2'] = errors(w, w1, w2, Lw, [3 * np.pi / 24, 3 * np.pi / 48, 3 * np.pi / 96, 3 * np.pi / 192, 3 * np.pi / 384])

# (c) dispersion symbol from a periodic matrix vs formula
h = 0.5; N = 64
P = np.zeros((N, N))
for j in range(N):
    for d, c in ((0, A0), (1, A1), (-1, A1), (2, A2), (-2, A2), (3, A3), (-3, A3)):
        P[j, (j + d) % N] += c / h ** 2
ks = 2 * np.pi * np.arange(N) / (N * h)
num = -np.real(np.array([np.exp(-1j * k * h * np.arange(N)) @ P @ np.exp(1j * k * h * np.arange(N)) / N for k in ks]))
res['dispersion_symbol_max_rel_error'] = float(np.max(abs(num - lattice_momentum(ks, h) ** 2) / np.maximum(1e-30, lattice_momentum(ks, h) ** 2 + 1e-300) * (ks > 0)))
kk = np.array([.01, .1, .5, 1.0])
for hh in (.5, .25):
    res[f'dispersion_error_k_minus_kh_over_k5h4_h{hh}'] = ((kk - lattice_momentum(kk, hh)) / (kk ** 5 * hh ** 4)).tolist()
res['dispersion_leading_coefficient_expected'] = 3 / 640

# (d) lead roots for omega^2 - m^2 = p^2 with p the lattice momentum of k
roots_rows = []
for hh in (.0625, .125, .25, .5, 2.0):
    for k in np.r_[np.geomspace(1e-4, 1e-2, 5), [.02, .05, .1, .2, .5, 1, 2, 3]]:
        if k * hh >= np.pi:
            continue
        p = lattice_momentum(k, hh)
        z = np.roots([A3, A2, A1, A0 + (p * hh) ** 2, A1, A2, A3])
        mod = np.abs(z)
        on = np.sum(abs(mod - 1) < 1e-8); inside = np.sum(mod < 1 - 1e-8); outside = np.sum(mod > 1 + 1e-8)
        ev = z[mod > 1 + 1e-8]
        osc = z[abs(mod - 1) < 1e-8]
        roots_rows.append(dict(h=hh, k=float(k), n_unit=int(on), n_inside=int(inside), n_outside=int(outside),
                               outside_moduli=sorted(abs(ev).tolist()),
                               outside_are_complex_pair=bool(len(ev) == 2 and abs(ev[0].imag) > 1e-6 and abs(ev[0] - ev[1].conj()) < 1e-8),
                               osc_phase_error=float(min(abs(np.angle(osc) - k * hh).min(), abs(np.angle(osc) + k * hh).min()))))
res['lead_roots'] = roots_rows
res['lead_roots_min_outside_modulus'] = float(min(min(r['outside_moduli']) for r in roots_rows))
res['lead_roots_all_2_2_2'] = bool(all(r['n_unit'] == 2 and r['n_inside'] == 2 and r['n_outside'] == 2 for r in roots_rows))

# (e) semi-discrete energy identity (sigma = 0): dH/dt = V.(LU - m^2 U) + (GU).(GV) + m^2 U.V
rng = np.random.default_rng(1)
n, h, m = 80, .25, 1.0
G = G_matrix(n, h); L = D_matrix(n, h) @ G
U, V = rng.normal(size=n - 1), rng.normal(size=n - 1)
dH = V @ (L @ U - m * m * U) + (G @ U) @ (G @ V) + m * m * U @ V
res['semi_discrete_energy_rate_relative'] = float(abs(dH) / (abs(V @ L @ U) + 1e-300))

(OUT / 'operators_sbp_order_roots.json').write_text(json.dumps(res, indent=1) + '\n')
short = {k: v for k, v in res.items() if not isinstance(v, list)}
print(json.dumps(short, indent=1))
for key in ('order_regular_odd_U', 'order_irregular_U_with_U2_at_0', 'order_wall_with_nonzero_U2'):
    print(key)
    for r in res[key]:
        print({a: (round(b, 3) if 'order' in a else b) for a, b in r.items()})
