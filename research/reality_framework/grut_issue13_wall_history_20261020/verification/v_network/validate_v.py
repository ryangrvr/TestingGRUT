#!/usr/bin/env python3
"""Validation of the verifier's estimators and integrator.
1. Wall-length estimators on exact tanh circles and tilted lines (2D), spheres (3D).
2. Static kink stays static (integrator sanity), energy-like check in flat space (alpha = 0).
"""
import json, sys, time
import numpy as np
from lattice_v import wall_measures, Sim, laplacian, U

out = {}
dx = 0.5
# circles
res = []
for R in (8.0, 16.0, 32.0):
    N = int(4 * R / dx) + 16
    x = (np.arange(N) - N / 2 + 0.37) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2)
    phi = np.tanh(r - R)
    m = wall_measures(phi, dx)
    exact = 2 * np.pi * R
    res.append({"R": R, **{k: m[k] / exact for k in m}})
out["circle_ratio_to_exact"] = res
# tilted periodic lines: phi = tanh-like periodic profile along direction (k1,k2)
res = []
N = 512
L = N * dx
x = np.arange(N) * dx
X, Y = np.meshgrid(x, x, indexing="ij")
for (k1, k2) in [(1, 0), (1, 1), (1, 2), (2, 3)]:
    kk = np.hypot(k1, k2)
    s = (k1 * X + k2 * Y) / kk  # distance along normal
    period = L / kk
    # two walls per period: phi = tanh of distance to nearest of positions 0.25,0.75 period
    u = np.mod(s, period) / period
    phi = np.tanh(np.where(u < 0.5, (u - 0.25), (0.75 - u)) * period)
    m = wall_measures(phi, dx)
    exact = 2 * kk * L  # 2 walls per period, kk periods across, each line of length L*... (total length = 2 * L^2 / period)
    exact = 2 * L * L / period
    res.append({"k": [k1, k2], **{k: m[k] / exact for k in m}})
out["line_ratio_to_exact"] = res
# sphere 3D
res = []
for R in (8.0, 16.0):
    N = int(2.6 * R / dx) + 8
    x = (np.arange(N) - N / 2 + 0.37) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    phi = np.tanh(r - R)
    m = wall_measures(phi, dx)
    exact = 4 * np.pi * R * R
    res.append({"R": R, **{k: m[k] / exact for k in m}})
out["sphere_ratio_to_exact"] = res

# integrator: static kink in flat space (alpha=0) with periodic pair of kinks
N = 256
x = np.arange(N) * dx
L = N * dx
phi0 = np.tanh(np.minimum(x - L / 4, 3 * L / 4 - x))[:, None] * np.ones((1, 8))
s = Sim(phi0, dx, alpha=0.0, beta=0.0, tau0=1.0)
for i in range(2000):
    s.step(0.15)
out["static_kink_maxdev_after_300"] = float(np.abs(s.phi - phi0).max())

# energy conservation in flat space for random field (alpha = 0), dt = 0.15 vs 0.075
def energy(sim, dt):
    v = sim.velocity(dt)
    g2 = 0
    for ax in range(sim.d):
        g2 = g2 + ((np.roll(sim.phi, -1, axis=ax) - sim.phi) / sim.dx) ** 2
    return float((0.5 * v ** 2 + 0.5 * g2 + U(sim.phi)).sum() * sim.dx ** sim.d)

from lattice_v import gaussian_ic
ec = []
for dt in (0.15, 0.075):
    phi0 = gaussian_ic(128, 2, dx, 0.75, 0.1, seed=5)
    s = Sim(phi0, dx, alpha=0.0, beta=0.0, tau0=1.0)
    s.step(dt)
    E0 = energy(s, dt)
    Es = []
    nst = int(round(60 / dt))
    for i in range(nst):
        s.step(dt)
        Es.append(energy(s, dt))
    Es = np.array(Es)
    ec.append({"dt": dt, "max_rel_dev": float(np.abs(Es - E0).max() / abs(E0)), "E0": E0})
out["flat_energy_check"] = ec
json.dump(out, open("results/validation_v.json", "w"), indent=1)
print(json.dumps(out, indent=1))
