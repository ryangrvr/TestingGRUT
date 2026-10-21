#!/usr/bin/env python3
"""Track wbound: 1D (planar-symmetric) flat-space lattice tests.  Units v = 1, lam = 0.01 (m_phi = 0.1414,
wall width delta = sqrt(2/lam) = 14.14), g^2 = 0.01 (m_X = 0.1), lam_X = 0.001.  Periodic boxes.

Families (all with phi and complex X):
  A  bounded one-sign random waves about +v (non-relativistic and relativistic spectra, with X waves)
  B  kink-antikink pairs at speed u (0, 0.5, 0.9, 0.99): wall w = u^2 - 2/3 and post-collision debris
  C  quench to a domain-wall network (random small phi, pi = 0) with X noise
  D  Codex-type kinetic lumps: phi = v, pi = A bump(z; R0) in a sparse set of slabs (instantaneous w = 1 at t = 0)
     with small (linear) and large (over-barrier) A and several R0, plus 1e-3 noise.
For each: time-averaged w, the exact-identity check <3p - rho> = lam <phi^2 - 1> + boundary, extrema of phi,
and the running cumulative wbar_T.  Writes results/sims_1d.json.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lattice as L  # noqa: E402

OUT = {}


def bump(z, R0):
    r = np.abs(z)
    return np.where(r < R0, np.cos(0.5 * np.pi * r / R0) ** 2, 0.0)


def series(rows, nmax=200):
    t, w = L.running_wbar(rows)
    idx = np.unique(np.linspace(0, len(t) - 1, nmax).astype(int))
    inst = np.array([r["p"] / r["rho"] for r in rows])
    return dict(t=t[idx].round(3).tolist(), wbar_T=w[idx].round(5).tolist(), w_inst=inst[idx].round(5).tolist(),
                frac_abs_gt1=[round(rows[i]["frac_abs_gt1"], 5) for i in idx],
                Phi0_over_rho=[round(rows[i]["Phi0"] / rows[i]["rho"], 5) for i in idx])


def family_A(rng):
    res = {}
    N, dx, dt, T = 8192, 0.5, 0.1, 3000.0
    for name, kc, amp, xamp in (("NR_small", 0.05, 0.05, 0.05), ("NR_large", 0.05, 0.3, 0.3),
                                ("rel_small", 1.0, 0.02, 0.02), ("rel_large", 1.0, 0.1, 0.1),
                                ("mixed_large_X", 0.3, 0.2, 0.8)):
        lat = L.Lattice(1, N, dx)
        lat.phi = 1 + L.gaussian_field(rng, 1, N, dx, kc, amp)
        omega = math.sqrt(2 * L.LAM + kc**2)
        lat.Pphi = L.gaussian_field(rng, 1, N, dx, kc, amp * omega)
        lat.x1 = L.gaussian_field(rng, 1, N, dx, kc, xamp)
        lat.x2 = L.gaussian_field(rng, 1, N, dx, kc, xamp)
        omx = math.sqrt(G2_ + kc**2)
        lat.P1 = L.gaussian_field(rng, 1, N, dx, kc, xamp * omx)
        lat.P2 = L.gaussian_field(rng, 1, N, dx, kc, xamp * omx)
        rows = L.run(lat, dt, T, every=10)
        s = L.summarize_flat(rows)
        s["late_half"] = L.summarize_flat(rows, t_from=T / 2)
        res[name] = s
        print("A", name, json.dumps({k: s[k] for k in ("wbar", "one_sign", "identity_residual_rel", "w_inst_first", "phi_min", "phi_max")}), flush=True)
    return res


G2_ = L.G2


def kink(z, u):
    g = 1 / math.sqrt(1 - u * u)
    d = math.sqrt(2 / L.LAM)
    return np.tanh(g * z / d), g / d / np.cosh(g * z / d) ** 2      # profile and d/dz (times -u gives phidot)


def family_B(rng):
    res = {}
    for u in (0.0, 0.5, 0.9, 0.99):
        g = 1 / math.sqrt(1 - u * u)
        dx = min(0.5, 14.14 / g / 12)
        Lbox = 3000.0
        N = int(Lbox / dx)
        dt = 0.4 * dx
        z = (np.arange(N) - N / 2) * dx
        z1, z2 = -300.0, 300.0
        f1, d1 = kink(z - z1, u)       # kink at z1 moving right (+u)
        f2, d2 = kink(z - z2, u)       # kink profile at z2, used as antikink moving left
        lat = L.Lattice(1, N, dx)
        lat.phi = f1 - f2 - 1.0
        lat.Pphi = -u * d1 - u * d2        # d/dt tanh(g(z-z1-ut)/d) = -u * d/dz ; antikink -tanh(g(z-z2+ut)/d): -(+u) d2
        lat.x1 = 1e-4 * rng.normal(size=N)
        lat.x2 = 1e-4 * rng.normal(size=N)
        tcol = (z2 - z1) / (2 * u) if u > 0 else 1e9
        T = 2400.0 if u > 0 else 600.0
        rows = L.run(lat, dt, T, every=max(1, int(1.0 / dt)))
        s = dict(gamma=g, dx=dx, t_collision=tcol if u > 0 else None,
                 NG_w=u * u - 2 / 3,
                 pre_collision=L.summarize_flat(rows, t_to=min(0.8 * tcol, T)),
                 whole=L.summarize_flat(rows))
        if u > 0:
            s["post_collision"] = L.summarize_flat(rows, t_from=tcol + 50)
            post = [r for r in rows if tcol - 5 <= r["t"] <= tcol + 200]
            s["max_abs_phi_near_collision"] = max(max(abs(r["phi_min"]), r["phi_max"]) for r in post)
            s["max_frac_abs_gt1_near_collision"] = max(r["frac_abs_gt1"] for r in post)
            s["max_Phi0_over_rho_after_collision"] = max(r["Phi0"] / r["rho"] for r in rows if r["t"] > tcol)
        s["series"] = series(rows)
        res[f"u={u}"] = s
        print("B", u, json.dumps({"NG_w": s["NG_w"], "pre": s["pre_collision"]["wbar"], "whole": s["whole"]["wbar"],
                                  "post": s.get("post_collision", {}).get("wbar"),
                                  "idres": s["whole"]["identity_residual_rel"]}), flush=True)
    return res


def family_E(rng):
    """Ultra-relativistic kink-antikink collisions (gamma = 22.4, 70.7) in a smaller box."""
    res = {}
    us = tuple(float(x) for x in os.environ.get("WB_E_SPEEDS", "0.999,0.9999").split(","))
    for u in us:     # u = 0.9999 (gamma 70.7) exceeded the 30-min cost cap and was not executed
        g = 1 / math.sqrt(1 - u * u)
        dx = 14.14 / g / 12
        Lbox = 1200.0
        N = int(Lbox / dx)
        dt = 0.4 * dx
        z = (np.arange(N) - N / 2) * dx
        z1, z2 = -100.0, 100.0
        f1, d1 = kink(z - z1, u)
        f2, d2 = kink(z - z2, u)
        lat = L.Lattice(1, N, dx)
        lat.phi = f1 - f2 - 1.0
        lat.Pphi = -u * d1 - u * d2
        lat.x1 = 1e-4 * rng.normal(size=N)
        lat.x2 = 1e-4 * rng.normal(size=N)
        tcol = (z2 - z1) / (2 * u)
        T = 560.0
        rows = L.run(lat, dt, T, every=max(1, int(1.0 / dt)))
        post = [r for r in rows if r["t"] >= tcol - 5]
        s = dict(gamma=g, dx=dx, N=N, t_collision=tcol, NG_w=u * u - 2 / 3,
                 pre_collision=L.summarize_flat(rows, t_to=0.8 * tcol),
                 post_collision=L.summarize_flat(rows, t_from=tcol + 20),
                 post_late=L.summarize_flat(rows, t_from=tcol + 250),
                 whole=L.summarize_flat(rows),
                 max_abs_phi_after_collision=max(max(abs(r["phi_min"]), r["phi_max"]) for r in post),
                 max_Phi0_over_rho_after_collision=max(r["Phi0"] / r["rho"] for r in post),
                 series=series(rows))
        res[f"u={u}"] = s
        print("E", u, json.dumps({"gamma": g, "pre": s["pre_collision"]["wbar"], "post": s["post_collision"]["wbar"],
                                  "post_late": s["post_late"]["wbar"], "maxphi": s["max_abs_phi_after_collision"],
                                  "maxPhi0/rho": s["max_Phi0_over_rho_after_collision"],
                                  "idres": s["whole"]["identity_residual_rel"]}), flush=True)
    return res


def family_C(rng):
    res = {}
    N, dx, dt, T = 16384, 0.5, 0.1, 4000.0
    for name, kc, amp in (("quench_kc0.3", 0.3, 0.1), ("quench_kc1_hot", 1.0, 0.5)):
        lat = L.Lattice(1, N, dx)
        lat.phi = L.gaussian_field(rng, 1, N, dx, kc, amp)
        lat.Pphi = L.gaussian_field(rng, 1, N, dx, kc, 0.3 * amp * kc)
        lat.x1 = L.gaussian_field(rng, 1, N, dx, kc, 0.3)
        lat.x2 = L.gaussian_field(rng, 1, N, dx, kc, 0.3)
        rows = L.run(lat, dt, T, every=10)
        s = dict(whole=L.summarize_flat(rows), late=L.summarize_flat(rows, t_from=T / 2), series=series(rows))
        nw = int(np.sum(np.diff(np.sign(lat.phi)) != 0))
        s["n_sign_changes_end"] = nw
        res[name] = s
        print("C", name, s["whole"]["wbar"], s["late"]["wbar"], nw, s["whole"]["identity_residual_rel"], flush=True)
    return res


def family_D(rng):
    res = {}
    for name, A, R0, nslab, Lbox in (("lin_A0.02_R5", 0.02, 5.0, 4, 4000.0), ("lin_A0.02_R50", 0.02, 50.0, 4, 4000.0),
                                     ("over_A0.2_R20", 0.2, 20.0, 4, 4000.0), ("over_A0.2_R100", 0.2, 100.0, 4, 4000.0),
                                     ("over_A0.2_R400", 0.2, 400.0, 2, 4000.0), ("over_A0.12_R100", 0.12, 100.0, 4, 4000.0)):
        dx, dt = 0.5, 0.1
        N = int(Lbox / dx)
        z = np.arange(N) * dx
        centers = (np.arange(nslab) + 0.5) * Lbox / nslab
        lat = L.Lattice(1, N, dx)
        prof = sum(bump(z - c, R0) for c in centers)
        # seed noise confined to the slabs (smooth, kc = 0.5) so the lumps carry essentially all the energy
        lat.phi = 1 + 1e-3 * prof * L.gaussian_field(rng, 1, N, dx, 0.5, 1.0)
        lat.Pphi = A * prof
        lat.x1 = 1e-3 * prof * L.gaussian_field(rng, 1, N, dx, 0.5, 1.0)
        lat.x2 = 1e-3 * prof * L.gaussian_field(rng, 1, N, dx, 0.5, 1.0)
        T = 3000.0
        rows = L.run(lat, dt, T, every=10)
        t, w = L.running_wbar(rows)
        above = t[w > 1 / 3]
        s = dict(A=A, R0=R0, nslab=nslab, filling=float(prof.mean()),
                 local_eps_over_lam=A * A / 2 / L.LAM, mean_rho=rows[0]["rho"],
                 whole=L.summarize_flat(rows), late=L.summarize_flat(rows, t_from=1500.0),
                 first_window=L.summarize_flat(rows, t_to=min(2 * R0, T)),
                 last_t_wbarT_above_third=float(above.max()) if above.size else None,
                 max_wbarT_after_t10=float(w[t > 10].max()),
                 series=series(rows))
        res[name] = s
        print("D", name, json.dumps({"w0": s["whole"]["w_inst_first"], "whole": s["whole"]["wbar"], "late": s["late"]["wbar"],
                                     "first": s["first_window"]["wbar"], "last_above": s["last_t_wbarT_above_third"],
                                     "idres": s["whole"]["identity_residual_rel"]}), flush=True)
    return res


if __name__ == "__main__":
    rng = np.random.default_rng(13)
    which = sys.argv[1:] or ["A", "B", "C", "D"]
    for fam in which:
        OUT[fam] = {"A": family_A, "B": family_B, "C": family_C, "D": family_D, "E": family_E}[fam](rng)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "results", "sims_1d_%s.json" % "".join(which)), "w") as fh:
        json.dump(OUT, fh, indent=1)
