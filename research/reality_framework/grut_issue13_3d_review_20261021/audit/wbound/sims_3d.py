#!/usr/bin/env python3
"""Track wbound: 3D lattice tests (reduced grids, numpy).  Units v = 1, lam = 0.01, g^2 = 0.01, lam_X = 0.001.

  net      quench to a domain-wall network with X noise (flat, periodic, N = 96, dx = 1.5)
  lump     one isolated over-barrier kinetic lump (phi = v, pi = A cos^2 bump of radius R0) in N = 128, dx = 1.5;
           records the dilatation moment Dm = int x . T^{0i}; checks dDm/dt = 3 P_tot and the bound
           wbar_T <= 1/3 + (Dm(T) - Dm(0) - E T ... ) i.e. (1/T) int 3P dt = (Dm(T) - Dm(0))/T <= E (R0 + T)/T,
           then continues in the periodic box (a cubic lattice of lumps = sparse gas)
  small    small-amplitude kinetic lumps (linear), N = 64
  charged  homogeneous rotating X condensate at phi ~ 0 (restored, |X| = 1.5/sqrt2) plus noise, N = 64
  frw      prescribed radiation a = (t/t0)^(1/2), H0/m_phi = 0.1; network and one-sign wave states, N = 64;
           checks the exact budget  Delta(rho a^4) = -int H a^4 <Phi_0> dt - [H a^4 g0] + int a'' a^3 g0 dt
Writes results/sims_3d_<test>.json.
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lattice as L  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def save(name, obj):
    with open(os.path.join(HERE, "results", f"sims_3d_{name}.json"), "w") as fh:
        json.dump(obj, fh, indent=1)


def series(rows, nmax=200, keys=("frac_abs_gt1", "frac_neg")):
    t, w = L.running_wbar(rows)
    idx = np.unique(np.linspace(0, len(t) - 1, nmax).astype(int))
    out = dict(t=t[idx].round(3).tolist(), wbar_T=w[idx].round(5).tolist(),
               w_inst=[round(rows[i]["p"] / rows[i]["rho"], 5) for i in idx],
               Phi0_over_rho=[round(rows[i]["Phi0"] / rows[i]["rho"], 5) for i in idx])
    for k in keys:
        out[k] = [round(rows[i][k], 5) for i in idx]
    return out


def bump3(lat, center, R0):
    N, dx = lat.N, lat.dx
    r2 = 0
    for ax in range(3):
        c = (np.arange(N) - center[ax] + N / 2) % N - N / 2
        shp = [1, 1, 1]
        shp[ax] = N
        r2 = r2 + (c.reshape(shp) * dx) ** 2
    r = np.sqrt(r2)
    return np.where(r < R0, np.cos(0.5 * np.pi * r / R0) ** 2, 0.0)


def t_net(rng):
    N, dx, dt, T = 96, 1.5, 0.3, 900.0
    lat = L.Lattice(3, N, dx)
    lat.phi = L.gaussian_field(rng, 3, N, dx, 0.15, 0.1)
    lat.x1 = L.gaussian_field(rng, 3, N, dx, 0.3, 0.2)
    lat.x2 = L.gaussian_field(rng, 3, N, dx, 0.3, 0.2)
    rows = L.run(lat, dt, T, every=5)
    out = dict(whole=L.summarize_flat(rows), early=L.summarize_flat(rows, t_to=150.0),
               late=L.summarize_flat(rows, t_from=450.0),
               max_inst_Phi0_over_rho=max(r["Phi0"] / r["rho"] for r in rows),
               max_inst_w=max(r["p"] / r["rho"] for r in rows),
               series=series(rows))
    save("net", out)
    return {k: out[k] for k in ("whole", "early", "late", "max_inst_Phi0_over_rho", "max_inst_w")}


def t_lump(rng, A=0.2, R0=30.0, N=128, dx=1.5, dt=0.3, T=450.0):
    lat = L.Lattice(3, N, dx)
    c = (N / 2, N / 2, N / 2)
    b = bump3(lat, c, R0)
    lat.phi = 1 + 1e-3 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)      # seed noise confined to the lump
    lat.Pphi = A * b
    lat.x1 = 1e-3 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)
    lat.x2 = 1e-3 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)
    t_iso = (N * dx / 2 - R0)          # light-cone of the bump edge reaches the box face
    rows = []

    def cb(lat_, rws):
        pass

    rows = L.run(lat, dt, T, every=3, center=None)
    # second pass for the dilatation moment would double cost; instead recompute isolated phase separately
    out = dict(A=A, R0=R0, N=N, dx=dx, Lbox=N * dx, t_isolated=t_iso,
               E_tot=rows[0]["rho"] * (N * dx) ** 3, mean_rho=rows[0]["rho"],
               local_eps_over_lam=A * A / 2 / L.LAM,
               whole=L.summarize_flat(rows), late=L.summarize_flat(rows, t_from=T / 2),
               series=series(rows))
    t, w = L.running_wbar(rows)
    above = t[w > 1 / 3]
    out["last_t_wbarT_above_third"] = float(above.max()) if above.size else None
    out["max_wbarT_after_t5"] = float(w[t > 5].max())
    save(f"lump_A{A}_R{R0}", out)
    return {k: out[k] for k in ("whole", "late", "last_t_wbarT_above_third", "max_wbarT_after_t5", "t_isolated")}


def tophat3(lat, center, R0, edge=8.0):
    N, dx = lat.N, lat.dx
    r2 = 0
    for ax in range(3):
        c = (np.arange(N) - center[ax] + N / 2) % N - N / 2
        shp = [1, 1, 1]
        shp[ax] = N
        r2 = r2 + (c.reshape(shp) * dx) ** 2
    r = np.sqrt(r2)
    x = np.clip((r - (R0 - edge)) / edge, 0, 1)
    return np.cos(0.5 * np.pi * x) ** 2


def t_lump_dilatation(rng, A=0.2, R0=30.0, N=128, dx=1.5, dt=0.3, profile="cos2"):
    """Isolated phase only, with the dilatation moment Dm(t) recorded."""
    lat = L.Lattice(3, N, dx)
    c = (N / 2, N / 2, N / 2)
    b = bump3(lat, c, R0) if profile == "cos2" else tophat3(lat, c, R0)
    lat.phi = 1 + 1e-3 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)      # seed noise confined to the lump
    lat.Pphi = A * b
    lat.x1 = 1e-3 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)
    lat.x2 = 1e-3 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)
    t_iso = (N * dx / 2 - R0)
    rows = L.run(lat, dt, t_iso, every=3, center=c)
    t = np.array([r["t"] for r in rows])
    E = np.array([r["E_tot"] for r in rows])
    P3 = np.array([3 * r["p"] for r in rows]) * (N * dx) ** 3
    Dm = np.array([r["Dm"] for r in rows])
    cumP3 = np.concatenate([[0], np.cumsum(0.5 * (P3[1:] + P3[:-1]) * np.diff(t))])
    # identity check: Dm(t) - Dm(0) = int 3 P_tot dt
    ident = np.max(np.abs((Dm - Dm[0]) - cumP3)) / (E[0] * t[-1])
    _, w = L.running_wbar(rows)
    bound = 1 / 3 + (R0 + abs(Dm[0]) / E[0]) / (3 * np.maximum(t, 1e-9))
    viol = bool(np.any((w[1:] > bound[1:] + 1e-3)))
    idx = np.unique(np.linspace(0, len(t) - 1, 120).astype(int))
    above = t[w > 1 / 3]
    out = dict(A=A, R0=R0, profile=profile, last_t_wbarT_above_third=float(above.max()) if above.size else None,
               t_isolated=t_iso, E0=float(E[0]), E_drift=float(E.std() / E.mean()),
               Dm0=float(Dm[0]), dilatation_identity_max_rel_err=float(ident), bound_violated=viol,
               t=t[idx].round(3).tolist(), wbar_T=w[idx].round(5).tolist(), bound=bound[idx].round(5).tolist(),
               Dm_over_E=(Dm[idx] / E[0]).round(4).tolist(),
               R0_plus_t=(R0 + t[idx]).round(3).tolist(),
               w_inst=[round(rows[i]["p"] / rows[i]["rho"], 5) for i in idx],
               frac_abs_gt1=[round(rows[i]["frac_abs_gt1"], 6) for i in idx])
    save(f"lumpdil_A{A}_R{R0}" + ("" if profile == "cos2" else "_" + profile), out)
    return {k: out[k] for k in ("E_drift", "Dm0", "dilatation_identity_max_rel_err", "bound_violated", "t_isolated",
                                "last_t_wbarT_above_third")} | \
        dict(wbar_T_end=out["wbar_T"][-1], max_wbar_T_after_5=float(w[t > 5].max()))


def t_small(rng):
    N, dx, dt, T = 64, 1.5, 0.3, 900.0
    res = {}
    for name, A, R0 in (("A0.02_R6", 0.02, 6.0), ("A0.02_R20", 0.02, 20.0)):
        lat = L.Lattice(3, N, dx)
        b = bump3(lat, (N / 2, N / 2, N / 2), R0)
        lat.Pphi = A * b
        lat.x1 = 1e-4 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)     # seed noise confined to the lump
        lat.x2 = 1e-4 * b * L.gaussian_field(rng, 3, N, dx, 0.5, 1.0)
        rows = L.run(lat, dt, T, every=3)
        s = dict(whole=L.summarize_flat(rows), late=L.summarize_flat(rows, t_from=T / 2),
                 first_10=L.summarize_flat(rows, t_to=10.0), series=series(rows))
        t, w = L.running_wbar(rows)
        above = t[w > 1 / 3]
        s["last_t_wbarT_above_third"] = float(above.max()) if above.size else None
        res[name] = s
    save("small", res)
    return {k: {kk: v[kk] for kk in ("whole", "first_10", "last_t_wbarT_above_third")} for k, v in res.items()}


def t_charged(rng):
    N, dx, dt, T = 64, 1.5, 0.3, 900.0
    r = 1.5
    om = math.sqrt(L.LX * r * r)          # circular orbit at phi = 0: omega^2 = 2 lam_X |X|^2 = lam_X r^2
    lat = L.Lattice(3, N, dx)
    lat.phi = 0.01 * rng.normal(size=(N,) * 3)
    lat.x1 = r + 0.01 * rng.normal(size=(N,) * 3)
    lat.x2 = 0.01 * rng.normal(size=(N,) * 3)
    lat.P2 = om * r * np.ones((N,) * 3)
    rows = L.run(lat, dt, T, every=5)
    out = dict(r=r, omega=om, whole=L.summarize_flat(rows), late=L.summarize_flat(rows, t_from=T / 2),
               max_inst_w=max(rw["p"] / rw["rho"] for rw in rows), series=series(rows))
    save("charged", out)
    return {k: out[k] for k in ("whole", "late", "max_inst_w")}


def frw_budget(rows, t0):
    t = np.array([r["t"] for r in rows])
    a = np.array([r["a"] for r in rows])
    H = np.array([r["H"] for r in rows])
    rho = np.array([r["rho"] for r in rows])
    p = np.array([r["p"] for r in rows])
    Phi0 = np.array([r["Phi0"] for r in rows])
    g0 = np.array([r["g0"] for r in rows])
    add = -a / (4 * t * t)                      # a'' for a = (t/t0)^(1/2)
    R = rho * a**4
    tr = lambda y: float(np.trapezoid(y, t))
    lhs = R[-1] - R[0]
    cont = -tr(H * a**4 * (3 * p - rho))
    term_phi = -tr(H * a**4 * Phi0)
    term_b = -(H[-1] * a[-1] ** 4 * g0[-1] - H[0] * a[0] ** 4 * g0[0])
    term_add = tr(add * a**3 * g0)
    return dict(a_end=float(a[-1]), rho_a4_start=float(R[0]), rho_a4_end=float(R[-1]), delta_rho_a4=lhs,
                continuity_rhs=cont, continuity_rel_err=abs(lhs - cont) / abs(R[0]),
                virial_term=term_phi, boundary_term=term_b, addot_term=term_add,
                virial_budget_rel_err=abs(lhs - (term_phi + term_b + term_add)) / abs(R[0]),
                boundary_over_rho_a4=abs(term_b) / R[0], addot_over_rho_a4=abs(term_add) / R[0],
                max_inst_Phi0_over_rho=float((Phi0 / rho).max()),
                rho_a4_monotone_fraction=float(np.mean(np.diff(R) >= -1e-12 * R[0])),
                rho_a4_ratio=float(R[-1] / R[0]))


def t_frw(rng):
    N, dx, dt = 64, 1.5, 0.15
    m = math.sqrt(2 * L.LAM)
    res = {}
    for name in ("network", "one_sign_waves"):
        for Hm in (0.1, 0.02):
            t0 = 0.5 / (Hm * m)
            lat = L.Lattice(3, N, dx, frw=dict(t0=t0))
            if name == "network":
                lat.phi = L.gaussian_field(rng, 3, N, dx, 0.15, 0.1)
                lat.x1 = L.gaussian_field(rng, 3, N, dx, 0.3, 0.2)
                lat.x2 = L.gaussian_field(rng, 3, N, dx, 0.3, 0.2)
            else:
                lat.phi = 1 + L.gaussian_field(rng, 3, N, dx, 0.3, 0.1)
                lat.Pphi = L.gaussian_field(rng, 3, N, dx, 0.3, 0.1 * math.sqrt(2 * L.LAM + 0.09))
                lat.x1 = L.gaussian_field(rng, 3, N, dx, 0.3, 0.1)
                lat.x2 = L.gaussian_field(rng, 3, N, dx, 0.3, 0.1)
            Tspan = 3.0 * t0 if Hm == 0.1 else 1.0 * t0     # a: 1 -> 2 (Hm 0.1); 1 -> 1.41 (Hm 0.02)
            rows = L.run(lat, dt, Tspan, every=2)
            b = frw_budget(rows, t0)
            b["H0_over_m"] = Hm
            b["one_sign"] = bool(min(r["phi_min"] for r in rows) > 0)
            res[f"{name}_Hm{Hm}"] = b
            print("frw", name, Hm, json.dumps({k: b[k] for k in ("rho_a4_ratio", "continuity_rel_err", "virial_budget_rel_err",
                                                                 "boundary_over_rho_a4", "max_inst_Phi0_over_rho", "one_sign")}), flush=True)
    save("frw", res)
    return res


if __name__ == "__main__":
    rng = np.random.default_rng(2026)
    tests = sys.argv[1:]
    fns = dict(net=t_net, lump=t_lump, lumpdil=t_lump_dilatation, small=t_small, charged=t_charged, frw=t_frw,
               tophat=lambda r: t_lump_dilatation(r, A=0.2, R0=40.0, profile="tophat"))
    for tname in tests:
        t0 = time.time()
        r = fns[tname](rng)
        print(tname, "%.0fs" % (time.time() - t0), json.dumps(r, default=str)[:3000], flush=True)
