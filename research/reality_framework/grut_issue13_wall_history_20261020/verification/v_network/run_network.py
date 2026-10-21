#!/usr/bin/env python3
"""Run one PRS network simulation with the verifier's own code and write a JSON time series.

A definition used here (stated): rho_w = A sigma / t with t cosmic time, a proportional to tau
(radiation era, t = a tau / 2).  Physical wall energy per physical volume (3D) or per physical area
(2D): rho_w = sigma (a^(d-1) S_c) / (a^d V_c) = sigma S_c / (a V_c), S_c = comoving wall area (3D)
or length (2D), V_c = comoving volume (area).  Hence
    A = rho_w t / sigma = (S_c / V_c) * tau / 2.
"""
import argparse, json, time
import numpy as np
from lattice_v import Sim, gaussian_ic, wall_measures, wall_v2

ap = argparse.ArgumentParser()
ap.add_argument("--d", type=int, default=2)
ap.add_argument("--N", type=int, default=512)
ap.add_argument("--dx", type=float, default=0.5)
ap.add_argument("--alpha", type=float, default=2.0)
ap.add_argument("--beta", type=float, default=0.0)
ap.add_argument("--tau0", type=float, default=1.0)
ap.add_argument("--tau_end", type=float, default=None)
ap.add_argument("--s_c", type=float, default=0.75)
ap.add_argument("--std", type=float, default=0.1)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--p", type=float, default=None)
ap.add_argument("--dt_frac", type=float, default=0.3)
ap.add_argument("--nout", type=int, default=48)
ap.add_argument("--out", type=str, required=True)
ap.add_argument("--maxwall", type=float, default=1440.0)
args = ap.parse_args()

d, N, dx = args.d, args.N, args.dx
L = N * dx
tau_end = args.tau_end if args.tau_end else L / 2
dt = args.dt_frac * dx
phi0 = gaussian_ic(N, d, dx, args.s_c, args.std, args.seed, mean_shift_frac=args.p)
sim = Sim(phi0, dx, args.alpha, args.beta, args.tau0)
V = L ** d
touts = list(np.geomspace(max(args.tau0 * 2, 2.0), tau_end, args.nout))
rec = []
t_start = time.time()
k = 0
while sim.tau < tau_end - 1e-9:
    sim.step(dt)
    while k < len(touts) and sim.tau >= touts[k] - 1e-9:
        m = wall_measures(sim.phi, dx)
        vel = sim.velocity(dt)
        row = {"tau": sim.tau,
               "fplus": float((sim.phi > 0).mean()),
               "core_frac": float((np.abs(sim.phi) < 0.5).mean()),
               "v2": wall_v2(sim.phi, vel, dx)}
        for key, val in m.items():
            row["S_" + key] = val
            row["A_" + key] = val / V * sim.tau / 2.0
        rec.append(row)
        k += 1
    if time.time() - t_start > args.maxwall:
        print("wall-clock cap reached at tau", sim.tau)
        break
res = {"config": vars(args), "L": L, "dt": dt, "tau_end_reached": sim.tau,
       "wall_seconds": time.time() - t_start, "series": rec}
json.dump(res, open(args.out, "w"), indent=0)
print("done", args.out, "tau", sim.tau, "secs", round(time.time() - t_start, 1))
