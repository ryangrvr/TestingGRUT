"""lattice_v.py -- verifier's own 3D PRS wall lattice: small-bias response and unbiased autocorrelation.

Independent implementation (not derived from modules/bias_threshold/code/lattice.py or network_sim/wallsim.py):
  phi'' + (alpha/tau) phi' - lap phi + U'(phi) = 0,  U = (1/8)(phi^2 - 1)^2  (v = m_phi = 1; PRS, beta = 0)
  alpha = 3 (3D radiation era, a ~ tau).  Friction handled exactly with P = tau^alpha phi':
      phi' = P / tau^alpha,   P' = tau^alpha (lap phi - U'(phi))     (kick-drift-kick leapfrog, dt = 0.25, dx = 1)
  IC at tau_i = 2: Gaussian field, white noise filtered with exp(-k^2 s^2/2), s = 1.5, std 0.1, exactly zero mean,
  zero velocity.  Bias: the same field shifted by +-mu (m0 = erf(mu/(sqrt2 * 0.1))).  My own RNG seeds (101...).
Measurements at log-spaced tau:
  nl   = number of sign-changing nearest-neighbour links (wall area up to a constant; K = nl(form)/nl(t))
  m    = mean sign(phi)
  C    = < sign(phi(t)) phi0 > / std(phi0)            (unbiased run only)
  frac = fraction of cells with |phi| < 1/2
Usage: python3 lattice_v.py N seed mu  -> results/lattice_v_N{N}_s{seed}_mu{mu}.json
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "lattice_v")


def ic(N, seed, s=1.5, amp=0.1):
    rng = np.random.default_rng(seed)
    w = rng.standard_normal((N, N, N))
    k = 2 * np.pi * np.fft.fftfreq(N)
    kz = 2 * np.pi * np.fft.rfftfreq(N)
    K2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + kz[None, None, :] ** 2
    f = np.fft.irfftn(np.fft.rfftn(w) * np.exp(-0.5 * K2 * s * s), s=(N, N, N))
    f -= f.mean()
    f *= amp / f.std()
    return f.astype(np.float32)


def lap(phi, out):
    out[...] = -6.0 * phi
    for ax in range(3):
        out += np.roll(phi, 1, axis=ax)
        out += np.roll(phi, -1, axis=ax)
    return out


def force(phi, out):
    lap(phi, out)
    out -= 0.5 * phi * (phi * phi - 1.0)
    return out


def measure(phi, phi0n):
    sg = np.sign(phi)
    nl = 0
    for ax in range(3):
        nl += int(np.count_nonzero(sg != np.roll(sg, 1, axis=ax)))
    return {"nl": nl, "m": float(sg.mean(dtype=np.float64)),
            "C": float(np.mean(sg * phi0n, dtype=np.float64)) if phi0n is not None else None,
            "frac_lt_half": float(np.mean(np.abs(phi) < 0.5))}


def run(N, seed, mu, tau_i=2.0, dt=0.25, alpha=3.0, n_out=40):
    tau_end = N / 2.0
    phi0 = ic(N, seed)
    phi0n = (phi0 / phi0.std()).astype(np.float32)
    phi = (phi0 + np.float32(mu)).astype(np.float32)
    P = np.zeros_like(phi)
    F = np.empty_like(phi)
    outs = set(np.unique(np.round(np.geomspace(tau_i, tau_end, n_out) / dt).astype(int)).tolist())
    nsteps = int(round((tau_end - tau_i) / dt))
    series = []
    t0 = time.time()
    tau = tau_i
    istep0 = int(round(tau_i / dt))
    force(phi, F)
    for n in range(nsteps + 1):
        if (istep0 + n) in outs or n == 0 or n == nsteps:
            mm = measure(phi, phi0n)
            mm["tau"] = tau
            series.append(mm)
        if n == nsteps:
            break
        # kick (half), drift, kick (half)
        P += np.float32(0.5 * dt * tau ** alpha) * F
        tau_mid = tau + 0.5 * dt
        phi += np.float32(dt / tau_mid ** alpha) * P
        tau = tau + dt
        force(phi, F)
        P += np.float32(0.5 * dt * tau ** alpha) * F
    return {"N": N, "seed": seed, "mu": mu, "m0_expected": math.erf(mu / (math.sqrt(2) * 0.1)), "dt": dt,
            "alpha": alpha, "tau_i": tau_i, "wall_s": time.time() - t0, "series": series}


if __name__ == "__main__":
    N, seed, mu = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    os.makedirs(OUT, exist_ok=True)
    r = run(N, seed, mu)
    fn = os.path.join(OUT, f"lattice_v_N{N}_s{seed}_mu{mu:+.5f}.json")
    json.dump(r, open(fn, "w"))
    print(fn, "%.1f s" % r["wall_s"])
