"""Validation of the charged_exit_sim integrator and diagnostics.  Writes results/validation.json.

V1 charge conservation to round-off in strongly nonlinear evolution (2D and 3D)
V2 energy-budget residual and its dt^2 convergence
V3 exact homogeneous phi = 0 rotating condensate: stays on the branch, n = n_i a^-3, |Y|^2 = (n^2/8B)^(1/3)
V4 linear tachyonic growth of a single phi mode through the n_c crossing vs the mode ODE
V5 conformal-variable integrator vs a direct integrator of f'' + 2 (a'/a) f' - lap f + a^2 dU/df* = 0
V6 periodic labelling / wrapping, core-component classification (wall vs Q-ball), wall length on known shapes
V7 exact Z2 mirror: phi -> -phi initial data gives the mirrored evolution (f+ <-> f-)
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cesim as C  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")


def evolve(sim, tau_i, tau_f, dt, with_X=True, record=None, a_fn=None):
    """Plain KDK loop returning (Q list, H residual list, tau list)."""
    a_of = a_fn or (lambda t: t / tau_i)
    n = int(round((tau_f - tau_i) / dt))
    tau = tau_i
    a = a_of(tau)
    H0 = sim.energy(a)[0]
    Q0 = sim.charge()
    ss = sim.kick(0.5 * dt, a, with_X)
    prev = sim.dHdtau(a, ss)
    budget = 0.0
    out = []
    for k in range(1, n + 1):
        sim.drift(dt, with_X)
        tau = tau_i + k * dt
        a = a_of(tau)
        ss = sim.kick(0.5 * dt, a, with_X)
        dH = sim.dHdtau(a, ss)
        budget += 0.5 * dt * (prev + dH)
        prev = dH
        if record and (k % record == 0 or k == n):
            E = sim.energy(a)
            out.append(dict(tau=tau, a=a, Q=sim.charge(), dQ=(sim.charge() - Q0) / abs(Q0) if Q0 else 0.0,
                            H=E[0], resid=E[0] - H0 - budget, kin_grad=E[1] + E[2], flow=budget))
        sim.kick(0.5 * dt, a, with_X)
    sim.kick(-0.5 * dt, a, with_X)
    return out


def v1_v2():
    res = {}
    for dt_fac in (1.0, 0.5):
        cfg = dict(d=2, N=128, dx=0.25, tau_i=10.0, kind="charged_exit", kappa=2.0, seed=7, eps_S=0.05, eps_X=0.02)
        sim, info = C.setup(cfg)
        dt = 0.05 * dt_fac
        out = evolve(sim, 10.0, 30.0, dt, record=int(round(1.0 / dt)))
        m = C.measure(sim, 3.0, 30.0)
        res[f"2D_dt{dt}"] = dict(max_abs_dQ_rel=max(abs(o["dQ"]) for o in out),
                                 final_resid_over_kin_grad=out[-1]["resid"] / out[-1]["kin_grad"],
                                 max_resid_over_kin_grad=max(abs(o["resid"] / o["kin_grad"]) for o in out),
                                 final_resid_over_H=out[-1]["resid"] / out[-1]["H"],
                                 max_abs_resid=max(abs(o["resid"]) for o in out),
                                 total_flow=out[-1]["flow"], final_f_plus=m["f_plus"], final_f_minus=m["f_minus"],
                                 final_qfrac_lumps=m["qfrac_lumps"])
    r1 = res["2D_dt0.05"]["max_abs_resid"]
    r2 = res["2D_dt0.025"]["max_abs_resid"]
    res["resid_ratio_dt_halving"] = r1 / r2
    # 3D
    cfg = dict(d=3, N=32, dx=0.3, tau_i=10.0, kind="charged_exit", kappa=2.0, seed=3, eps_S=0.05, eps_X=0.02)
    sim, info = C.setup(cfg)
    out = evolve(sim, 10.0, 25.0, 0.06, record=25)
    res["3D_N32"] = dict(max_abs_dQ_rel=max(abs(o["dQ"]) for o in out),
                         final_resid_over_kin_grad=out[-1]["resid"] / out[-1]["kin_grad"])
    return res


def v3():
    out = {}
    for kappa in (2.0, 4.0):
        cfg = dict(d=2, N=8, dx=0.25, tau_i=40.0, kind="charged_exit", kappa=kappa, seed=1, eps_S=0.0, eps_X=0.0)
        sim, info = C.setup(cfg)
        sim.S[...] = 0.0
        sim.PS[...] = 0.0
        sim.ghosts()
        dt = 0.05
        tau_i = 40.0
        n = int(round(80.0 / dt))
        sim.kick(0.5 * dt, 1.0)
        maxdev_mod = 0.0
        maxdev_rel = 0.0
        y0 = info["Y0"]
        nbar_i = info["nbar_i"]
        for k in range(1, n + 1):
            sim.drift(dt)
            a = 1 + k * dt / tau_i
            sim.kick(0.5 * dt, a)
            if k % 100 == 0:
                mod = math.sqrt(sim.R[0, 0] ** 2 + sim.I[0, 0] ** 2)
                maxdev_mod = max(maxdev_mod, abs(mod / y0 - 1))
                nphys = 2 * (sim.R[0, 0] * sim.PI[0, 0] - sim.I[0, 0] * sim.PR[0, 0]) / a ** 3
                y2 = (mod / a) ** 2
                pred = (nphys ** 2 / (8 * C.B_X)) ** (1 / 3)
                maxdev_rel = max(maxdev_rel, abs(y2 / pred - 1), abs(nphys / (nbar_i / a ** 3) - 1))
            sim.kick(0.5 * dt, a)
        out[f"kappa{kappa}"] = dict(max_rel_dev_conformal_modulus=maxdev_mod,
                                    max_rel_dev_branch2_relation_and_dilution=maxdev_rel,
                                    S_max_abs=float(np.abs(sim.S).max()),
                                    a_cross=kappa ** (1 / 3))
    return out


def v4():
    """Single phi mode on the homogeneous condensate: delta chi'' + (k_lat^2 + kappa^(2/3) - a^2) delta chi = 0."""
    kappa, tau_i, N, dx = 2.0, 40.0, 64, 0.25
    L = N * dx
    m = 2
    k = 2 * math.pi * m / L
    klat2 = (2 / dx * math.sin(k * dx / 2)) ** 2
    amp = 1e-9
    cfg = dict(d=2, N=N, dx=dx, tau_i=tau_i, kind="charged_exit", kappa=kappa, seed=1, eps_S=0.0, eps_X=0.0)
    sim, info = C.setup(cfg)
    x = np.arange(N) * dx
    sim.S[...] = amp * np.cos(k * x)[:, None]
    sim.PS[...] = 0.0
    sim.ghosts()
    dt = 0.025
    tau_f = 1.9 * tau_i
    n = int(round((tau_f - tau_i) / dt))
    sim.kick(0.5 * dt, 1.0)
    ts, amps = [], []
    for j in range(1, n + 1):
        sim.drift(dt)
        tau = tau_i + j * dt
        a = tau / tau_i
        sim.kick(0.5 * dt, a)
        if j % 40 == 0:
            ts.append(tau)
            amps.append(float(2 * np.mean(sim.S[:, 0] * np.cos(k * x))))
        sim.kick(0.5 * dt, a)
    c2 = kappa ** (2 / 3)

    def rhs(t, y):
        a = t / tau_i
        return [y[1], -(klat2 + c2 - a * a) * y[0]]
    sol = solve_ivp(rhs, (tau_i, tau_f), [amp, 0.0], t_eval=ts, rtol=1e-11, atol=1e-22)
    ode = sol.y[0]
    lin = [i for i in range(len(ts)) if abs(amps[i]) < 1e-3]
    rel = [abs(amps[i] / ode[i] - 1) for i in lin if abs(ode[i]) > 1e-3 * amp]
    return dict(mode_m=m, k=k, growth_factor_sim=abs(amps[lin[-1]]) / amp, growth_factor_ode=abs(ode[lin[-1]]) / amp,
                tau_last_linear=ts[lin[-1]], max_rel_dev=max(rel), a_cross=kappa ** (1 / 3),
                a_tachyonic_for_mode=math.sqrt(c2 + klat2))


def v5():
    """Direct integration in the physical-field variables (f, f') with the friction form, compared with the
    conformal-variable integrator, on a nonlinear 64^2 configuration."""
    cfg = dict(d=2, N=64, dx=0.25, tau_i=10.0, kind="charged_exit", kappa=2.0, seed=5, eps_S=0.05, eps_X=0.02)
    res = {}
    for dt in (0.02, 0.01):
        sim, info = C.setup(cfg)
        tau_i, tau_f = 10.0, 22.0
        # conformal reference
        evolve(sim, tau_i, tau_f, dt)
        a_f = tau_f / tau_i
        Sc = sim.S / a_f
        Rc = sim.R / a_f
        # direct: fields f = chi / a, f' = (chi' - a' f)/a
        sim2, _ = C.setup(cfg)
        S = sim2.S.copy(); R = sim2.R.copy(); I = sim2.I.copy()
        a = 1.0
        dS = (sim2.PS - S / tau_i) / a
        dR = (sim2.PR - R / tau_i) / a
        dI = (sim2.PI - I / tau_i) / a

        def lap(f):
            return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / 0.25 ** 2

        def acc(S, R, I, dS, dR, dI, tau):
            a = tau / tau_i
            h = 1.0 / tau            # a'/a
            y2 = R * R + I * I
            aS = lap(S) - a * a * S * (S * S - 1 + 2 * C.A_G * y2) - 2 * h * dS
            m2 = C.A_G * S * S + 2 * C.B_X * y2
            aR = lap(R) - a * a * m2 * R - 2 * h * dR
            aI = lap(I) - a * a * m2 * I - 2 * h * dI
            return aS, aR, aI
        # RK4 on the first-order system (independent of the leapfrog used above)
        n = int(round((tau_f - tau_i) / dt))
        tau = tau_i
        for _ in range(n):
            y = (S, R, I, dS, dR, dI)

            def F(y, t):
                aS, aR, aI = acc(*y, t)
                return (y[3], y[4], y[5], aS, aR, aI)
            k1 = F(y, tau)
            k2 = F(tuple(y[i] + 0.5 * dt * k1[i] for i in range(6)), tau + 0.5 * dt)
            k3 = F(tuple(y[i] + 0.5 * dt * k2[i] for i in range(6)), tau + 0.5 * dt)
            k4 = F(tuple(y[i] + dt * k3[i] for i in range(6)), tau + dt)
            S, R, I, dS, dR, dI = (y[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6))
            tau += dt
        res[f"dt{dt}"] = dict(max_abs_diff_S=float(np.abs(S - Sc).max()), max_abs_diff_ReY=float(np.abs(R - Rc).max()),
                             S_range=[float(Sc.min()), float(Sc.max())])
    return res


def v6():
    out = {}
    N = 64
    m = np.zeros((N, N), bool)
    m[10:20, :] = True                       # stripe wrapping along axis 1
    lab, n, sz, wr, sp = C.plabel(m)
    out["stripe"] = dict(n=n, wraps=wr.tolist(), spans=sp.tolist())
    m = np.zeros((N, N), bool)
    m[:5, :5] = True; m[-5:, :5] = True; m[:5, -5:] = True; m[-5:, -5:] = True   # corner blob through boundary
    lab, n, sz, wr, sp = C.plabel(m)
    out["corner_blob"] = dict(n=n, size=sz.tolist(), wraps=wr.tolist())
    i, j = np.indices((N, N))
    m = ((i + j) % N) < 6                    # diagonal stripe: wraps both axes
    lab, n, sz, wr, sp = C.plabel(m)
    out["diagonal_stripe"] = dict(n=n, wraps=wr.tolist())
    # synthetic field: two straight walls (periodic), one Q-ball inside the + domain, one Q-ball on a wall
    N = 256
    dx = 0.25
    sim = C.Sim(2, N, dx, 10.0)
    x = (np.arange(N) + 0.5) * dx
    L = N * dx
    X, Yg = np.meshgrid(x, x, indexing="ij")
    delta = math.sqrt(2.0)
    S = np.tanh((X - L / 4) / delta) * np.tanh((3 * L / 4 - X) / delta)    # + between L/4 and 3L/4
    for (cx, cy) in ((L / 2, L / 2), (L / 4, L / 4)):
        r = np.hypot(X - cx, Yg - cy)
        ball = 0.5 * (1 - np.tanh((r - 3.0) / 0.6))
        S = S * (1 - ball)
        Y = np.sqrt(1.6) * ball
        sim.R[...] += Y
        sim.PI[...] += 0.56 * Y
    sim.S[...] = S
    sim.ghosts()
    mm = C.measure(sim, 1.0, 10.0, qlist=True)
    out["synthetic"] = dict(L_wall_times_area=mm["L_wall"] * L * L, expected_two_walls=2 * L,
                            L_all_times_area=mm["L_all"] * L * L, n_lumps=mm["n_lumps"],
                            n_lumps_attached=mm["n_lumps_attached"], lump_attached_list=mm.get("lump_attached_list"),
                            n_core_walltype=mm["n_core_walltype"], n_core_single=mm["n_core_single"],
                            wraps_plus=mm["wraps_plus"], wraps_minus=mm["wraps_minus"])
    # tanh circle
    for R0 in (8.0, 16.0):
        r = np.hypot(X - L / 2, Yg - L / 2)
        sim.S[...] = np.tanh((r - R0) / delta)
        sim.R[...] = 0; sim.I[...] = 0; sim.PI[...] = 0
        mm = C.measure(sim, 1.0, 10.0)
        out[f"circle_R{R0}"] = dict(L=mm["L_wall"] * L * L, exact=2 * math.pi * R0,
                                    rel_err=mm["L_wall"] * L * L / (2 * math.pi * R0) - 1,
                                    wraps_plus=mm["wraps_plus"], n_comp_minus=mm["n_comp_minus"])
    return out


def v7():
    cfg = dict(d=2, N=128, dx=0.25, tau_i=20.0, kind="charged_exit", kappa=2.0, seed=11)
    sim, _ = C.setup(cfg)
    sim2, _ = C.setup(cfg)
    sim2.S[...] *= -1
    sim2.PS[...] *= -1
    sim2.ghosts()
    evolve(sim, 20.0, 50.0, 0.05)
    evolve(sim2, 20.0, 50.0, 0.05)
    m1 = C.measure(sim, 2.5, 50.0)
    m2 = C.measure(sim2, 2.5, 50.0)
    return dict(max_abs_S_plus_mirror=float(np.abs(sim.S + sim2.S).max()),
                max_abs_ReY_diff=float(np.abs(sim.R - sim2.R).max()),
                f_plus=m1["f_plus"], f_minus=m1["f_minus"], mirror_f_plus=m2["f_plus"], mirror_f_minus=m2["f_minus"])


if __name__ == "__main__":
    t0 = time.time()
    res = {}
    for name, fn in (("V1_V2_charge_energy", v1_v2), ("V3_homogeneous_branch2", v3), ("V4_linear_tachyon", v4),
                     ("V5_direct_vs_conformal", v5), ("V6_labelling_and_estimators", v6), ("V7_Z2_mirror", v7)):
        t = time.time()
        res[name] = fn()
        res[name]["runtime_s"] = time.time() - t
        print(name, json.dumps(res[name])[:600], flush=True)
    res["total_runtime_s"] = time.time() - t0
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "validation.json"), "w") as f:
        json.dump(res, f, indent=1)
