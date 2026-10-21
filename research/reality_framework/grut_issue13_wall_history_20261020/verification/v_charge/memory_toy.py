#!/usr/bin/env python3
"""Independent version of the spherical sign-memory toy (verifier v_charge). Own construction:
Codex units (length ell = 1/(sqrt(lp) v)), a = 1.  S_tt = S_rr + (2/r) S_r - (S^2 - 1 + a x(r,t)) S.
Pre-quench x = x_before * w(r), w = (1 - tanh((r-R)/0.6))/2, S = static positive solution (gradient flow).
At t = 0: x -> x_after * w(r) (sudden), seed dS = -eps * w(r) (coherent, my choice), zero velocity.
Decision = sign of S(0,t) when |S(0,t)| first exceeds 0.5.  R_first_flip = smallest R on the grid with
negative decision.  Grid R step 1.  Also reports the Dirichlet threshold pi/sqrt(1 - a x_after).
"""
import json, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_memory.json")
a = 1.0
dr = 0.05


def lap(S, r):
    L = np.empty_like(S)
    L[1:-1] = (S[2:] - 2 * S[1:-1] + S[:-2]) / dr**2 + (S[2:] - S[:-2]) / (dr * r[1:-1])
    L[0] = 6 * (S[1] - S[0]) / dr**2
    L[-1] = 0.0
    return L


_cache = {}


def static_profile(R, xb, rmax):
    key = (R, xb, rmax)
    if key in _cache:
        return _cache[key]
    r = np.arange(0, rmax + dr / 2, dr)
    w = 0.5 * (1 - np.tanh((r - R) / 0.6))
    S = np.ones_like(r)
    dt = 0.2 * dr**2
    for it in range(int(60 / dt)):
        S += dt * (lap(S, r) - (S**2 - 1 + a * xb * w) * S)
        S[-1] = 1.0
    _cache[key] = (r, w, S.copy())
    return r, w, S.copy()


def run(R, eps, xb=3.4, xa=0.0, tmax=70.0):
    rmax = R + 80.0
    r, w, S = static_profile(R, xb, rmax)
    S0c = S[0]
    S = S - eps * w
    Sd = np.zeros_like(S)
    dt = 0.25 * dr
    x = xa * w
    acc = lap(S, r) - (S**2 - 1 + a * x) * S
    t = 0.0
    while t < tmax:
        Sd += 0.5 * dt * acc
        S += dt * Sd
        S[-1] = 1.0
        acc = lap(S, r) - (S**2 - 1 + a * x) * S
        Sd += 0.5 * dt * acc
        t += dt
        if abs(S[0]) > 0.5:
            return int(np.sign(S[0])), t, S0c
    return 0, t, S0c


res = {}
for eps in [1e-2, 1e-4, 1e-6]:
    for xa in [0.0, 0.5]:
        rows = []
        first = None
        for R in range(3, 31):
            sgn, td, s0 = run(float(R), eps, xa=xa)
            rows.append(dict(R=R, decision=sgn, t_decision=td, S0_static=s0))
            if sgn < 0 and first is None:
                first = R
        res[f"xa={xa},eps={eps}"] = dict(R_first_flip=first, R_D=np.pi / np.sqrt(1 - a * xa), rows=rows)
        print(xa, eps, first, flush=True)
# heuristic of the module, evaluated
heur = {f"eps={e}": 2 + np.log(0.5 / e) / np.sqrt(a * 3.4 - 1) + np.log(1 / e) / 1.0 for e in [1e-2, 1e-4, 1e-6]}
eps_quantum = np.sqrt(0.01) / (2 * np.sqrt(2) * np.pi)
with open(OUT, "w") as f:
    json.dump(dict(results=res, heuristic_xa0=heur, eps_quantum_estimate=eps_quantum), f, indent=1, default=float)
print(heur, eps_quantum)
