#!/usr/bin/env python3
"""Verifier v_charge: archived Q = 39,944 soliton S-profile (read-only CSV), charge removed instantly (x = 0),
coherent negative seed -eps for r < 3.3, zero velocity; sign of S(0,t) when |S(0)| first exceeds 0.5.
Own spherical leapfrog (same toy class as memory_toy.py). Units ell = 1/(sqrt(lp) v) = 1/m_X at the benchmark.
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROF = os.path.join(HERE, "../../../grut_candidate_1_20261018/code/results_profile.csv")
d = np.loadtxt(PROF, delimiter=",", skiprows=1)
r0, F0, S0 = d[:, 0], d[:, 1], d[:, 3]
dr = 0.05
r = np.arange(0, 120 + dr / 2, dr)
S_init = np.interp(r, r0, S0, right=1.0)


def lap(S):
    L = np.empty_like(S)
    L[1:-1] = (S[2:] - 2 * S[1:-1] + S[:-2]) / dr**2 + (S[2:] - S[:-2]) / (dr * r[1:-1])
    L[0] = 6 * (S[1] - S[0]) / dr**2
    L[-1] = 0.0
    return L


rows = []
for eps in [0.0, 0.01, 0.05, 0.1, 0.3, 0.5]:
    S = S_init - eps * (r < 3.3)
    Sd = np.zeros_like(S)
    dt = 0.25 * dr
    acc = lap(S) - (S**2 - 1) * S
    t = 0.0
    dec = None
    while t < 60:
        Sd += 0.5 * dt * acc; S += dt * Sd; S[-1] = 1.0
        acc = lap(S) - (S**2 - 1) * S
        Sd += 0.5 * dt * acc; t += dt
        if abs(S[0]) > 0.5:
            dec = (int(np.sign(S[0])), t)
            break
    rows.append(dict(eps=eps, S0_initial=float(S_init[0] - eps), decision=dec))
    print(rows[-1], flush=True)
json.dump(dict(profile_S0=float(S0[0]), profile_F2_0=float(F0[0] ** 2), rows=rows),
          open(os.path.join(HERE, "results_archived_quench.json"), "w"), indent=1)
