"""Lattice-evolved test of network_sim's wall-speed estimator.

Separated phi^4 kinks/antikinks with chosen speeds are evolved in flat space on a periodic 1D lattice with
network_sim's discretisation (3-point Laplacian, leapfrog, U = (1/8)(phi^2-1)^2, delta = 2). After a settling time,
the estimator  sum_{|phi|<1/2} pi^2 / sum_{|phi|<1/2} (central grad phi)^2  is evaluated exactly as in
wallsim.measure (pi is the half-step momentum, as in wallsim.evolve), averaged over many outputs, and compared with
the area-weighted and energy-weighted mean squared speeds of the kinks. Kink speeds are measured from the tracked
zero crossings (not assumed), so radiation shed at initialisation is accounted for.

Usage: python3 -I lattice_kinks.py <out_json>
"""
import json
import math
import sys

import numpy as np

LAM = 0.5
DELTA = 2.0


def dU(phi):
    return LAM * phi * (phi * phi - 1.0)


def init(speeds, dx, seg_len, rng):
    n_seg = int(round(seg_len / dx))
    phis, pis, centers = [], [], []
    for i, v in enumerate(speeds):
        g = 1.0 / math.sqrt(1.0 - v * v)
        off = rng.uniform(-0.5, 0.5)
        x = (np.arange(n_seg) - n_seg / 2 + off) * dx
        s = 1.0 if i % 2 == 0 else -1.0
        # alternate the direction of motion too, so that neighbours do not systematically approach
        vv = v if (i // 2) % 2 == 0 else -v
        u = g * x / DELTA
        phis.append(s * np.tanh(u))
        pis.append(-s * vv * g / DELTA / np.cosh(u) ** 2)
        centers.append((i * n_seg + n_seg / 2 - off) * dx)
    return np.concatenate(phis), np.concatenate(pis), centers, n_seg


def zero_crossings(phi, dx):
    s = phi > 0
    idx = np.nonzero(s != np.roll(s, -1))[0]
    out = []
    for j in idx:
        a, b = phi[j], phi[(j + 1) % phi.size]
        out.append((j + a / (a - b)) * dx)
    return np.array(out)


def run(speeds, dx, courant=1 / 3, seg_len=240.0, t_settle=20.0, t_meas=40.0, seed=1):
    rng = np.random.default_rng(seed)
    phi, pi, centers, n_seg = init(speeds, dx, seg_len, rng)
    N = phi.size
    dt = courant * dx
    # leapfrog: pi at half steps; start pi_{-1/2} ~ pi_0 + (dt/2) F
    def force(p):
        lap = (np.roll(p, 1) + np.roll(p, -1) - 2 * p) / dx ** 2
        return lap - dU(p)
    pi = pi - 0.5 * dt * force(phi)
    n_settle = int(t_settle / dt)
    n_meas = int(t_meas / dt)
    num_t, ests, zc_hist = [], [], []
    for n in range(n_settle + n_meas + 1):
        if n >= n_settle and (n - n_settle) % max(1, int(1.0 / dt)) == 0:
            c = (np.roll(phi, -1) - np.roll(phi, 1)) / (2 * dx)
            core = np.abs(phi) < 0.5
            num = float(np.sum(pi[core] ** 2))
            den = float(np.sum(c[core] ** 2))
            ests.append(num / den)
            zc_hist.append((n * dt, zero_crossings(phi, dx)))
        pi = pi + dt * force(phi)
        phi = phi + dt * pi
    # measured speeds from first and last tracked crossings (same ordering; no collisions by construction)
    (t0, z0), (t1, z1) = zc_hist[0], zc_hist[-1]
    Ltot = N * dx
    vm = []
    if len(z0) == len(z1) == len(speeds):
        d = (z1 - z0 + Ltot / 2) % Ltot - Ltot / 2
        vm = list(np.abs(d) / (t1 - t0))
    vm = np.array(vm) if len(vm) else np.array(speeds)
    g = 1 / np.sqrt(1 - np.minimum(vm, 0.999999) ** 2)
    s_A = float(np.mean(vm ** 2))
    s_E = float(np.sum(g * vm ** 2) / np.sum(g))
    est = float(np.mean(ests))
    return {"dx": dx, "speeds_set": list(map(float, speeds)), "speeds_measured": list(map(float, vm)),
            "s_area": s_A, "s_energy": s_E, "estimator_time_mean": est, "estimator_time_std": float(np.std(ests)),
            "rel_to_s_energy": est / s_E - 1 if s_E > 0 else None, "rel_to_s_area": est / s_A - 1 if s_A > 0 else None}


def main(out):
    cases = {
        "one_speed_0.456": [0.456] * 8,
        "rest_and_0.9": [0.0, 0.9] * 4,
        "uniform_0_to_0.8": list(np.linspace(0.0, 0.8, 8)),
        "rest_and_0.6": [0.0, 0.6] * 4,
    }
    res = {}
    for name, sp in cases.items():
        res[name] = {}
        for dx in (1.0, 0.5, 0.25):
            r = run(sp, dx)
            res[name][str(dx)] = r
            print(f"{name:18s} dx={dx:5.2f} est={r['estimator_time_mean']:.4f}+-{r['estimator_time_std']:.4f} "
                  f"s_E={r['s_energy']:.4f} s_A={r['s_area']:.4f} relE={r['rel_to_s_energy']:+.3f} "
                  f"relA={r['rel_to_s_area'] if r['rel_to_s_area'] is None else round(r['rel_to_s_area'],3)} "
                  f"v_meas={[round(v,3) for v in r['speeds_measured']]}", flush=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
