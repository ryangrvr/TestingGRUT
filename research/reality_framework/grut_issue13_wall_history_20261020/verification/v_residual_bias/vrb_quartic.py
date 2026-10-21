"""Verifier v_residual_bias, part 5: is a large coherent quartic oscillation stable over ~1e3 oscillations?
(Module's Hubble-only counterfactual: amplitude 99 GeV >> soft rms 2.9 GeV at the KZ point, |m0| = 1 claimed.)

Minkowski / conformal lambda phi^4 with zero mass (RD conformal frame of the quartic regime):
    f'' = lap f - f^3,  units sqrt(lambda) A = 1, initial f = 1 + noise, f' = 0.
2D lattice 192^2, dx = 0.5, velocity Verlet dt = 0.05 (the Lame resonance band is per mode, so 2D suffices).  Two noise levels (white per cell, rms 1e-3 and 1e-6).
Output: coherent amplitude envelope max|<f>| per 7.4 time units (one period = 7.4163) and the rms of f - <f>.
"""
import json
import math
import os
import sys
import time

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def lap3(f, dx):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / (dx * dx)


def run(noise, N=192, dx=0.5, dt=0.05, tmax=3000.0, seed=5):
    rng = np.random.default_rng(seed)
    f = 1.0 + noise * rng.normal(size=(N, N))
    v = np.zeros_like(f)
    a = lap3(f, dx) - f ** 3
    per = 7.4163
    env = []
    cur, t, nxt = 0.0, 0.0, per
    nst = int(tmax / dt)
    for it in range(nst):
        v += 0.5 * dt * a
        f += dt * v
        a = lap3(f, dx) - f ** 3
        v += 0.5 * dt * a
        t += dt
        if it % 5 == 0:
            cur = max(cur, abs(float(f.mean())))
        if t >= nxt:
            env.append([round(t, 2), cur, float(f.std())])
            cur = 0.0
            nxt += per
    return env


if __name__ == "__main__":
    t0 = time.time()
    res = {}
    for nz in (1e-3, 1e-6):
        env = run(nz)
        res["noise_%g" % nz] = env
        e = np.array(env)
        print(nz, [(round(x[0]), round(x[1], 4), round(x[2], 4)) for x in env[::25]], flush=True)
    res["wall_s"] = time.time() - t0
    with open(os.path.join(OUT, "quartic.json"), "w") as fh:
        json.dump(res, fh)
