"""Independent checks of what network_sim's wall-speed estimator  sum_{|phi|<1/2} pi^2 / sum_{|phi|<1/2} |grad_c phi|^2
measures, written without reusing audit/estimators/code.

Parts
  A. Offset-averaged continuum-sampled planar kinks in 3D at arbitrary orientation n (quadrature over the normal
     coordinate; lattice points equidistributed in n.x).  Gives the pure-discretisation bias of the estimator as a
     function of dx/delta, gamma and orientation, and its isotropic orientation average.  Also the energy weighting
     test for mixtures (per-wall num and den are additive).
  B. Lattice-evolved planar kinks on the exact lattice reductions of 3D walls:
       normal (1,0,0) -> 1D lattice with h = dx,  (1,1,0) -> h = dx/sqrt2,  (1,1,1) -> h = dx/sqrt3
     (for these normals the 7-point Laplacian and the summed central gradients reduce exactly to 1D with step h).
     Own leapfrog, U = (1/8)(phi^2-1)^2, delta = 2, Courant 1/3 in units of dx as in network_sim.
     Kink speeds from tracked zero crossings.
  C. network_sim 3D PRS series: v_rms_win provenance, window values, max v2_wall, and the dx = 1 vs 0.5 pair.

Usage: python3 -I verify_estimator.py <network_sim_dir> <out_json>
"""
import json
import math
import sys

import numpy as np

DELTA = 2.0
LAM = 0.5


# ------------------------------------------------------------------------------------------------- part A
def plane_num_den(v, n, dx, nu=20001):
    """Per-unit-area numerator and denominator (times dx^3) for one boosted kink with unit normal n."""
    g = 1 / math.sqrt(1 - v * v)
    uc = DELTA * math.atanh(0.5) / g             # |phi| < 1/2  <=>  |u| < uc
    u = np.linspace(-uc, uc, nu)
    f = lambda y: np.tanh(g * y / DELTA)
    pi = -v * g / DELTA / np.cosh(g * u / DELTA) ** 2
    den_int = np.zeros_like(u)
    for ni in n:
        den_int += ((f(u + ni * dx) - f(u - ni * dx)) / (2 * dx)) ** 2
    num = np.trapezoid(pi ** 2, u)
    den = np.trapezoid(den_int, u)
    den_exact = np.trapezoid((g / DELTA / np.cosh(g * u / DELTA) ** 2) ** 2, u)
    return num, den, den_exact


def iso_directions(n_dir=400, seed=3):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n_dir, 3))
    return x / np.linalg.norm(x, axis=1, keepdims=True)


def part_A():
    dirs = iso_directions()
    out = {"orientation_dependence_v0.456": {}, "isotropic_average": {}, "mixtures_isotropic": {}}
    named = {"(1,0,0)": (1, 0, 0), "(1,1,0)": (1 / math.sqrt(2), 1 / math.sqrt(2), 0),
             "(1,1,1)": (1 / math.sqrt(3),) * 3}
    for dx in (1.0, 0.5, 0.25):
        row = {}
        for lab, n in named.items():
            num, den, _ = plane_num_den(0.456, n, dx)
            row[lab] = num / den / 0.456 ** 2 - 1
        out["orientation_dependence_v0.456"][str(dx)] = row
    for v in (0.2, 0.456, 0.7, 0.9):
        row = {}
        for dx in (1.0, 0.5, 0.25):
            nums = dens = 0.0
            for n in dirs:
                a, b, _ = plane_num_den(v, n, dx, nu=4001)
                nums += a
                dens += b
            row[str(dx)] = nums / dens / v ** 2 - 1
        out["isotropic_average"][str(v)] = row
    # mixtures with equal area per wall, isotropic orientation; estimator = sum num / sum den
    cases = {"one_speed_0.456": [0.456] * 4, "rest_and_0.9": [0.0, 0.9], "uniform_0_0.8": list(np.linspace(0, 0.8, 9)),
             "rest_and_0.6": [0.0, 0.6]}
    for name, sp_ in cases.items():
        sp_ = np.array(sp_)
        g = 1 / np.sqrt(1 - sp_ ** 2)
        s_A = float(np.mean(sp_ ** 2))
        s_E = float(np.sum(g * sp_ ** 2) / np.sum(g))
        row = {"s_A": s_A, "s_E": s_E}
        for dx in (1.0, 0.5, 0.25, 0.1):
            N = D = 0.0
            for v in sp_:
                for n in dirs[:100]:
                    a, b, _ = plane_num_den(max(v, 1e-9), n, dx, nu=2001)
                    N += a if v > 0 else 0.0
                    D += b
            est = N / D
            row[str(dx)] = {"est": est, "rel_s_E": est / s_E - 1, "rel_s_A": est / s_A - 1}
        out["mixtures_isotropic"][name] = row
    return out


# ------------------------------------------------------------------------------------------------- part B
def evolve_1d(speeds, h, dt, seg=200.0, t_settle=25.0, t_meas=50.0, seed=11):
    """Own 1D leapfrog for separated kinks; lattice step h along the normal, time step dt."""
    rng = np.random.default_rng(seed)
    nseg = int(round(seg / h))
    phi_l, pi_l = [], []
    for i, v in enumerate(speeds):
        g = 1 / math.sqrt(1 - v * v)
        x = (np.arange(nseg) - nseg / 2 + rng.uniform(-0.5, 0.5)) * h
        sgn = 1.0 if i % 2 == 0 else -1.0
        direc = 1.0 if (i // 2) % 2 == 0 else -1.0
        phi_l.append(sgn * np.tanh(g * x / DELTA))
        pi_l.append(-sgn * direc * v * g / DELTA / np.cosh(g * x / DELTA) ** 2)
    phi = np.concatenate(phi_l)
    pi = np.concatenate(pi_l)

    def F(p):
        return (np.roll(p, 1) + np.roll(p, -1) - 2 * p) / h ** 2 - LAM * p * (p * p - 1)

    # staggered momentum: pi_{1/2} = pi_0 + dt/2 F(phi_0); measure uses the half-step pi as in network_sim
    pi = pi + 0.5 * dt * F(phi)
    ns, nm = int(t_settle / dt), int(t_meas / dt)
    every = max(1, int(round(1.0 / dt)))
    ests, zcs, ts = [], [], []
    for k in range(ns + nm + 1):
        phi = phi + dt * pi
        if k >= ns and (k - ns) % every == 0:
            c = (np.roll(phi, -1) - np.roll(phi, 1)) / (2 * h)
            m = np.abs(phi) < 0.5
            ests.append(float(np.sum(pi[m] ** 2) / np.sum(c[m] ** 2)))
            s = phi > 0
            idx = np.nonzero(s != np.roll(s, -1))[0]
            z = np.array([(j + phi[j] / (phi[j] - phi[(j + 1) % phi.size])) * h for j in idx])
            zcs.append(z)
            ts.append(k * dt)
        pi = pi + dt * F(phi)
    Lt = phi.size * h
    z0, z1 = zcs[0], zcs[-1]
    if len(z0) != len(speeds) or len(z1) != len(speeds):
        return {"error": "crossing count changed", "n0": len(z0), "n1": len(z1)}
    d = (z1 - z0 + Lt / 2) % Lt - Lt / 2
    vm = np.abs(d) / (ts[-1] - ts[0])
    g = 1 / np.sqrt(1 - np.minimum(vm, 0.999999) ** 2)
    s_A = float(np.mean(vm ** 2))
    s_E = float(np.sum(g * vm ** 2) / np.sum(g))
    est = float(np.mean(ests))
    s_set = float(np.sum((1 / np.sqrt(1 - np.array(speeds) ** 2)) * np.array(speeds) ** 2) /
                  np.sum(1 / np.sqrt(1 - np.array(speeds) ** 2)))
    return {"h": h, "dt": dt, "v_measured": [round(float(x), 4) for x in vm], "s_A": s_A, "s_E": s_E,
            "s_E_set": s_set, "est": est, "est_std_t": float(np.std(ests)),
            "rel_s_E_measured": est / s_E - 1 if s_E > 0 else None,
            "rel_s_A_measured": est / s_A - 1 if s_A > 0 else None,
            "rel_s_E_set": est / s_set - 1}


def part_B():
    out = {}
    dt_net = 1.0 / 3.0       # network_sim dt = dx/3 with dx = 1 (3D PRS 256^3 and 384^3)
    for lab, h in (("(1,0,0) h=dx=1", 1.0), ("(1,1,0) h=dx/sqrt2", 1 / math.sqrt(2)),
                   ("(1,1,1) h=dx/sqrt3", 1 / math.sqrt(3)), ("dx=0.5 (1,0,0)", 0.5), ("dx=0.5 (1,1,1)", 0.5 / math.sqrt(3))):
        dt = dt_net if "dx=0.5" not in lab else dt_net / 2
        row = {}
        for name, sp_ in (("one_speed_0.456", [0.456] * 4), ("rest_and_0.9", [0.0, 0.9, 0.0, 0.9]),
                          ("uniform_0_0.8", list(np.linspace(0.0, 0.8, 8)))):
            row[name] = evolve_1d(sp_, h, dt)
        out[lab] = row
    return out


# ------------------------------------------------------------------------------------------------- part C
def part_C(ns_dir):
    out = {}
    runs = {}
    for name in ("prs3d_N256_L256_s1", "prs3d_N256_L256_s2", "prs3d_N256_L256_s3", "prs3d_N384_L384_s1",
                 "prs3d_N256_L128_s1", "prs3d_bias_N128_L128_p0.50_s1", "phys3d_N256_L24_s1", "phys3d_N256_L24_s2",
                 "phys3d_N256_L24_s3", "phys3d_N128_L24_s1", "prs2d_a2_N2048_L2048_s1", "prs2d_a2_N2048_L2048_s2",
                 "prs2d_a2_N2048_L2048_s3", "prs2d_a3_N2048_L2048_s1"):
        r = json.load(open(f"{ns_dir}/results/runs/{name}.json"))
        L = r["config"]["L"]
        ser = r["series"]
        t = np.array([x["tau"] for x in ser])
        v2 = np.array([x["v2_wall"] if x["v2_wall"] is not None else np.nan for x in ser])
        h = t / L
        w = (h >= 0.125 * (1 - 1e-9)) & (h <= 0.25 * (1 + 1e-9))
        runs[name] = (t, v2)
        out[name] = {"dx": r["dx"], "window_tau": [float(t[w][0]), float(t[w][-1])] if w.any() else None,
                     "window_v2": [float(x) for x in v2[w]], "v_rms_win": float(np.sqrt(np.nanmean(v2[w]))) if w.any() else None,
                     "max_v2": float(np.nanmax(v2)), "tau_at_max": float(t[np.nanargmax(v2)]),
                     "n_v2_gt1": int(np.nansum(v2 > 1))}
    vals = [out[f"prs3d_N256_L256_s{i}"]["v_rms_win"] for i in (1, 2, 3)]
    out["prs3d_256_seed_stats"] = {"values": vals, "mean": float(np.mean(vals)), "std_ddof1": float(np.std(vals, ddof=1)),
                                   "mean_of_v2": float(np.mean(np.array(vals) ** 2))}
    t1, a1 = runs["prs3d_bias_N128_L128_p0.50_s1"]
    t2, a2 = runs["prs3d_N256_L128_s1"]
    common = sorted(set(np.round(t1, 6)) & set(np.round(t2, 6)))
    pair = []
    for tc in common:
        i = int(np.argmin(abs(t1 - tc)))
        j = int(np.argmin(abs(t2 - tc)))
        pair.append([float(tc), float(a1[i]), float(a2[j]), float(a1[i] / a2[j] - 1)])
    out["pair_dx1_vs_dx0.5_L128"] = pair
    late = [p for p in pair if p[0] >= 22]
    out["pair_rel_range_tau_ge_22"] = [min(p[3] for p in late), max(p[3] for p in late)] if late else None
    return out


def main(ns_dir, out_path):
    res = {"A_sampled_3D": part_A()}
    print(json.dumps(res["A_sampled_3D"], indent=1), flush=True)
    res["C_network_sim"] = part_C(ns_dir)
    print(json.dumps({k: v for k, v in res["C_network_sim"].items() if k != "pair_dx1_vs_dx0.5_L128"}, indent=1), flush=True)
    print(res["C_network_sim"]["pair_dx1_vs_dx0.5_L128"], flush=True)
    res["B_evolved"] = part_B()
    for lab, row in res["B_evolved"].items():
        for name, r in row.items():
            print(lab, name, {k: r.get(k) for k in ("v_measured", "est", "s_E", "s_A", "s_E_set", "rel_s_E_measured", "rel_s_A_measured", "rel_s_E_set")}, flush=True)
    json.dump(res, open(out_path, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
