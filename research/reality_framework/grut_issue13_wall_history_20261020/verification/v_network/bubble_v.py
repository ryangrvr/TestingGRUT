#!/usr/bin/env python3
"""(d) Isolated circular bubble collapse in 2D, field theory with the verifier's lattice code.
phi = tanh((r - R0)/w) (negative inside), released from rest at tau_i = R0/x, PRS alpha=2, beta=0
(D_h = 2) and also flat space (alpha = 0).  Compare with NG: Delta tau_c = (pi/2) R0, Delta tau_10 = arccos(0.1) R0
(flat), and the radiation-era NG ODE (ng_v.py) for the given x."""
import json, sys
import numpy as np
from lattice_v import Sim
sys.path.insert(0, ".")
from ng_v import collapse

out = []
CASES = [(20.0, 0.01, 2.0, 0.5), (20.0, 0.01, 2.0, 0.25), (20.0, 0.01, 2.0, 0.125),
         (40.0, 0.01, 2.0, 0.5), (40.0, 0.01, 2.0, 0.25), (40.0, 0.01, 2.0, 0.125),
         (80.0, 0.01, 2.0, 0.5), (80.0, 0.01, 2.0, 0.25),
         (40.0, 0.0, 0.0, 0.5), (40.0, 0.0, 0.0, 0.25), (40.0, 0.0, 0.0, 0.125),
         (40.0, 0.5, 2.0, 0.5), (40.0, 0.5, 2.0, 0.25)]
for (R0, x, alpha, dx) in CASES:
    L = 2.6 * R0 + 24
    N = int(np.ceil(L / dx / 16) * 16)
    xs = (np.arange(N) - N / 2 + 0.5) * dx
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    r = np.hypot(X, Y)
    phi0 = np.tanh(r - R0)
    tau_i = R0 / x if x > 0 else 1.0
    sim = Sim(phi0, dx, alpha, 0.0, tau_i)
    dt = 0.3 * dx
    rows = []
    t10 = None
    t50 = None
    t0 = None
    A0 = (phi0 < 0).sum() * dx * dx
    while sim.tau < tau_i + 3.0 * R0:
        sim.step(dt)
        Aneg = (sim.phi < 0).sum() * dx * dx
        Reff = np.sqrt(Aneg / np.pi)
        rows.append((sim.tau - tau_i, Reff))
        if t10 is None and Reff <= 0.1 * R0:
            # linear interpolation in R_eff between last two points
            (ta, Ra), (tb, Rb) = rows[-2], rows[-1]
            t10 = ta + (Ra - 0.1 * R0) * (tb - ta) / (Ra - Rb)
        if t50 is None and Reff <= 0.5 * R0:
            (ta, Ra), (tb, Rb) = rows[-2], rows[-1]
            t50 = ta + (Ra - 0.5 * R0) * (tb - ta) / (Ra - Rb)
        if t0 is None and Aneg == 0:
            t0 = sim.tau - tau_i
            break
    rows = np.array(rows)
    if x > 0:
        tc_ng, t10_ng, ti = collapse(1, 2.0 if alpha == 2.0 else 0.0, x, R0)
        tc_ng -= ti
        t10_ng -= ti
        t50_ng = None
    else:
        tc_ng, t10_ng = np.pi / 2 * R0, np.arccos(0.1) * R0
        t50_ng = np.arccos(0.5) * R0
    rec = {"R0_over_w": R0, "x": x, "alpha": alpha, "N": N, "dx": dx, "dt": dt,
           "Reff0_lattice": float(np.sqrt(A0 / np.pi)),
           "dtau_10_field": t10, "dtau_10_NG": float(t10_ng),
           "dtau_10_field_over_R0": t10 / R0, "dtau_10_NG_over_R0": float(t10_ng / R0),
           "dtau_c_field(first no negative cell)": t0, "dtau_c_field_over_R0": (t0 / R0) if t0 else None,
           "dtau_c_NG_over_R0": float(tc_ng / R0),
           "rel_dev_tau10": (t10 - t10_ng) / t10_ng,
           "dtau_50_field_over_R0": t50 / R0,
           "dtau_50_NG_over_R0_flat": (t50_ng / R0) if t50_ng else None,
           "rel_dev_tauc": ((t0 - tc_ng) / tc_ng) if t0 else None,
           "R_of_tau_sample": rows[:: max(1, len(rows) // 40)].tolist()}
    out.append(rec)
    print({k: v for k, v in rec.items() if k != "R_of_tau_sample"})
json.dump(out, open("results/bubble_v.json", "w"), indent=1)
