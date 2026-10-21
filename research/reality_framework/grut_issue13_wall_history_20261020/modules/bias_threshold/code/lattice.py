"""lattice.py -- linear response of the wall-network magnetization to a small systematic initial bias.

Module bias_threshold of the issue #13 wall-history packet.  numpy only, deterministic seeds, float32 fields.

Question answered here (computed, [C]): with an initial Gaussian field of short-range correlations and a uniform
shift that gives the expected initial positive-volume fraction p = 1/2 + m0/2 (m0 = 2p - 1, the initial magnetization of
the sign field), how does the sign-field magnetization m(tau) = 2 f_+(tau) - 1 grow with the coarsening factor
K(tau) = xi_ref(tau)/xi_ref(tau_form) of the unbiased network?  Phase-ordering theory (Bray, Adv. Phys. 43 (1994)
357, §4.2 eqs. (91)-(95) and §5.2 "External fields/initial bias") gives, for small m0 and dynamical scaling,

    m(tau) = B m0 K^gamma,   gamma = d - lambda_bar,

with lambda_bar the Fisher-Huse autocorrelation exponent, A(t) = <s(x,t) phi(x,0)> ~ L^-lambda_bar.  The identity
S(k,t,0) = Delta G(k,t) behind gamma = d - lambda_bar only uses the Gaussian initial ensemble, so it holds for
relativistic (second order in time) dynamics too; lambda_bar itself is dynamics-dependent and is measured here.

Dynamics (units v = m_phi = 1, U = (1/8)(phi^2 - 1)^2, planar wall tanh(z/2)):
  'prs'    : phi'' + (alpha/tau) phi' - lap phi + U'(phi) = 0   (Press-Ryden-Spergel, beta = 0; relativistic,
             radiation era; alpha = 3 in 3D, 2 in 2D), PRS leapfrog, dt = dx/3, tau_i = 2.
             Same equation, IC family and discretisation as modules/network_sim (independent implementation).
  'modelA' : phi_t = lap phi - U'(phi)  (overdamped, non-conserved; the Bray/Fisher-Huse reference class),
             explicit Euler, dt = 0.1.
Paired runs from ONE noise realization: the exactly zero-mean field (m0 = 0) and the same field shifted by +-mu,
m0 = erf(mu/sqrt2), so R = (m_+ - m_-)/(2 m0) is a common-random-numbers estimate of the response G(k=0, t).
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

LAM = 0.5
DELTA = 2.0          # planar-wall half width (tanh(z/DELTA)) in 1/m_phi
CORE = 2.0 * math.atanh(0.5) * DELTA   # comoving thickness of the |phi| < 1/2 core of a planar wall


def _erfinv(y):
    from scipy.special import erfinv
    return erfinv(y)


def lap_inplace(phi, out):
    d = phi.ndim
    np.multiply(phi, -2.0 * d, out=out)
    for ax in range(d):
        def s(a, b):
            t = [slice(None)] * d
            t[ax] = slice(a, b)
            return tuple(t)
        out[s(1, None)] += phi[s(None, -1)]
        out[s(0, 1)] += phi[s(-1, None)]
        out[s(None, -1)] += phi[s(1, None)]
        out[s(-1, None)] += phi[s(0, 1)]
    return out


def dU_inplace(phi, out):
    np.multiply(phi, phi, out=out)
    out -= 1.0
    out *= phi
    out *= LAM
    return out


def gaussian_field(d, N, s_c, seed):
    """White noise smoothed by a Gaussian kernel of std s_c (dx = 1), zero mean, unit std (float64)."""
    rng = np.random.default_rng(seed)
    w = rng.standard_normal((N,) * d)
    fk = np.fft.rfftn(w)
    k2 = np.zeros(fk.shape)
    for ax in range(d):
        k = 2 * np.pi * (np.fft.rfftfreq(N) if ax == d - 1 else np.fft.fftfreq(N))
        shp = [1] * d
        shp[ax] = k.size
        k2 = k2 + k.reshape(shp) ** 2
    fk *= np.exp(-0.5 * k2 * s_c ** 2)
    fk.flat[0] = 0.0
    f = np.fft.irfftn(fk, s=(N,) * d, axes=tuple(range(d)))
    f -= f.mean()
    f /= f.std()
    return f


def n_sign_links(pos):
    n = 0
    for ax in range(pos.ndim):
        n += int(np.count_nonzero(pos != np.roll(pos, -1, axis=ax)))
    return n


def run(cfg):
    """One run.  cfg: d, N, dyn ('prs'|'modelA'), alpha, s_c, amp, seed, p, t_end, n_out, name."""
    t0 = time.time()
    d, N = cfg["d"], cfg["N"]
    f0 = gaussian_field(d, N, cfg["s_c"], cfg["seed"])          # unit-std initial field shape (float64)
    s0 = np.sign(f0)
    # bias: uniform shift mu (in units of the initial std) of an exactly zero-mean field; the expected initial
    # sign magnetization is m0 = erf(mu/sqrt2) (Gaussian), i.e. p = (1 + m0)/2.  m0 = 0: no shift at all.
    m0 = cfg["m0"]
    shift = math.sqrt(2.0) * float(_erfinv(m0)) if m0 != 0 else 0.0
    phi = (cfg["amp"] * (f0 + shift)).astype(np.float32)
    f0n = f0.astype(np.float32)
    V = float(N ** d)
    iso = math.pi / 4 if d == 2 else 2.0 / 3.0
    buf = np.empty_like(phi)
    buf2 = np.empty_like(phi)
    if cfg["dyn"] == "prs":
        dt, t_i = 1.0 / 3.0, 2.0
        pi = np.zeros_like(phi)
    else:
        dt, t_i = 0.1, 0.0
        pi = None
    nsteps = int(math.ceil((cfg["t_end"] - t_i) / dt - 1e-9))
    t_out = np.geomspace(max(t_i, 1.0), cfg["t_end"], cfg["n_out"])
    out_steps = sorted(set(int(round((t - t_i) / dt)) for t in t_out) | {0, nsteps})
    out_set = set(out_steps)
    series = []
    for n in range(nsteps + 1):
        t = t_i + n * dt
        if n in out_set:
            pos = phi > 0
            fp = float(np.count_nonzero(pos)) / V
            nl = n_sign_links(pos)
            area = nl * iso                                        # comoving area (dx = 1), isotropic link count
            s = np.where(pos, np.float32(1.0), np.float32(-1.0))
            core = float(np.count_nonzero(np.abs(phi) < 0.5)) / V
            series.append({
                "t": t, "f_plus": fp, "m": 2 * fp - 1, "area_over_V": area / V,
                "xi": (V / area) if area > 0 else None,
                "A_prs": (area / V) * t / 2 if cfg["dyn"] == "prs" else None,
                "C_s_phi0": float(np.mean(s * f0n, dtype=np.float64)),    # <sign phi(t) * phi0/std>
                "C_s_s0": float(np.mean(s * s0, dtype=np.float64)),       # <sign phi(t) * sign phi0>
                "bulk_core_frac": core - CORE * area / V,
                "phi_min": float(phi.min()), "phi_max": float(phi.max()),
            })
            if area == 0 and len(series) > 2 and series[-2]["xi"] is None:
                break
        if n == nsteps:
            break
        lap_inplace(phi, buf)
        dU_inplace(phi, buf2)
        buf -= buf2
        if cfg["dyn"] == "prs":
            delta = 0.5 * cfg["alpha"] * dt / t
            pi *= (1.0 - delta)
            buf *= dt
            pi += buf
            pi *= 1.0 / (1.0 + delta)
            np.multiply(pi, dt, out=buf2)
            phi += buf2
        else:
            buf *= dt
            phi += buf
    return {"config": cfg, "shift_in_std_units": shift, "dt": dt, "wall_seconds": time.time() - t0,
            "series": series}


def jobs():
    """The job list.  Sizes chosen to keep the whole module under ~10 CPU-min on 2 workers."""
    J = []
    def add(tag, d, N, dyn, alpha, t_end, seeds, m0s, n_out, with_p0=True):
        for sd in seeds:
            labs = ([("m0", 0.0)] if with_p0 else []) + [(f"mp{m:g}", m) for m in m0s] + [(f"mm{m:g}", -m) for m in m0s]
            for lab, m in labs:
                J.append({"name": f"{tag}_s{sd}_{lab}", "group": tag, "d": d, "N": N, "dyn": dyn, "alpha": alpha,
                          "s_c": 1.5, "amp": 0.1, "seed": sd, "m0": m, "t_end": t_end, "n_out": n_out})
    # relativistic (PRS) networks in a radiation era
    add("prs3d_N160", 3, 160, "prs", 3.0, 80.0, [11, 12, 13], [0.002], 40)
    add("prs3d_N160", 3, 160, "prs", 3.0, 80.0, [11, 12], [0.02], 40, with_p0=False)    # nonlinear / linearity check
    add("prs2d_N1024", 2, 1024, "prs", 2.0, 512.0, [11, 12, 13], [0.001], 48)
    # overdamped model-A reference (validation of the estimator against the known lambda_bar)
    add("modelA2d_N512", 2, 512, "modelA", None, 1000.0, [21, 22], [0.002], 48)
    add("modelA3d_N96", 3, 96, "modelA", None, 250.0, [21], [0.002], 40)
    return J


def main(argv):
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "lattice_runs")
    os.makedirs(outdir, exist_ok=True)
    J = jobs()
    wid, nw = (int(argv[1]), int(argv[2])) if len(argv) >= 3 else (0, 1)
    for i, cfg in enumerate(J):
        if i % nw != wid:
            continue
        fn = os.path.join(outdir, cfg["name"] + ".json")
        if os.path.exists(fn):
            continue
        r = run(cfg)
        with open(fn, "w") as f:
            json.dump(r, f, indent=0)
        print(f"[{cfg['name']}] done in {r['wall_seconds']:.1f} s, last t = {r['series'][-1]['t']:.1f}, "
              f"m = {r['series'][-1]['m']:.4f}", flush=True)


if __name__ == "__main__":
    main(sys.argv)
