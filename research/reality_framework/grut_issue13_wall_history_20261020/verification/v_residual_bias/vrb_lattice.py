"""Verifier v_residual_bias, part 4: reduced independent check of route D (charge forced into the +v vacuum,
Minkowski limit), 2D classical lattice.  Own implementation and own seeds (not the module's).

Units: S = phi/v, Y = X/v, ell = 1/(sqrt(lambda) v) = 1/(100 GeV).  Ubar = (S^2-1)^2/4 + A S^2|Y|^2 + B|Y|^4,
A = g^2/lambda = 1, B = lambda_X/lambda = 0.1.  Equations (Minkowski):
   S'' = lap S - S (S^2 - 1 + 2 A |Y|^2),   Y'' = lap Y - (A S^2 + 2 B |Y|^2) Y.
n_c = sqrt(8 B) |Y_c|^3 with |Y_c|^2 = 1/(2A).  Preparation: S = 1 + smoothed noise (rms 1e-2, s_c = 1),
Y on its circular orbit in the S = 1 vacuum with 2 omega |Y|^2 = kappa n_c, omega^2 = A + 2B|Y|^2, noise 1e-3.
Velocity Verlet, 5-point Laplacian, dx = 0.5, dt = 0.1, N = 384 (L = 192 ell), tau_max = 300.
Measures: f+ (S > 1/2), f- (S < -1/2); 'formation' = first output after f_core >= 1/2 with both f+- >= 5%;
p_window = mean f+/(f+ + f-) over [tau_f + 5, tau_f + 25]; late p over the last 100.
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "lattice")
os.makedirs(OUT, exist_ok=True)
A_, B_ = 1.0, 0.1
NC = math.sqrt(8 * B_) * (1 / (2 * A_)) ** 1.5


def noise(N, L, sc, amp, rng):
    w = rng.normal(size=(N, N))
    k = 2 * np.pi * np.fft.fftfreq(N, L / N)
    K2 = k[:, None] ** 2 + k[None, :] ** 2
    f = np.fft.ifft2(np.fft.fft2(w) * np.exp(-0.5 * K2 * sc * sc)).real
    f -= f.mean()
    return amp * f / f.std()


def lap(f, dx):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / (dx * dx)


def run(kappa, seed, N=384, dx=0.5, dt=0.1, tmax=300.0, every=5):
    rng = np.random.default_rng(777000 + 31 * seed + int(kappa * 10))
    L = N * dx
    n = kappa * NC
    y2 = brentq(lambda q: 2 * math.sqrt(A_ + 2 * B_ * q) * q - n, 0, 50)
    Y0, om = math.sqrt(y2), math.sqrt(A_ + 2 * B_ * y2)
    S = 1 + noise(N, L, 1.0, 1e-2, rng)
    R = Y0 + noise(N, L, 1.0, 1e-3, rng)
    I = noise(N, L, 1.0, 1e-3, rng)
    vS = np.zeros_like(S)
    vR = -om * I
    vI = om * R

    def acc(S, R, I):
        q = R * R + I * I
        m = A_ * S * S + 2 * B_ * q
        return lap(S, dx) - S * (S * S - 1 + 2 * A_ * q), lap(R, dx) - m * R, lap(I, dx) - m * I
    aS, aR, aI = acc(S, R, I)
    ser = []
    nst = int(round(tmax / dt))
    for it in range(nst + 1):
        if it % every == 0:
            fp = float((S > 0.5).mean())
            fm = float((S < -0.5).mean())
            ser.append((it * dt, float(S.mean()), fp, fm, 1 - fp - fm, float((2 * (R * vI - I * vR)).mean() / NC)))
        if it == nst:
            break
        vS += 0.5 * dt * aS
        vR += 0.5 * dt * aR
        vI += 0.5 * dt * aI
        S += dt * vS
        R += dt * vR
        I += dt * vI
        aS, aR, aI = acc(S, R, I)
        vS += 0.5 * dt * aS
        vR += 0.5 * dt * aR
        vI += 0.5 * dt * aI
    ser = np.array(ser)
    t, mS, fp, fm, fc, nn = ser.T
    restored = np.maximum.accumulate(fc >= 0.5)
    idx = np.where(restored & (fp >= 0.05) & (fm >= 0.05))[0]
    res = dict(kappa=kappa, seed=seed, omega=om, Y0=Y0, n_over_nc_end=float(nn[-1]))
    if idx.size:
        tf = t[idx[0]]
        w = (t >= tf + 5) & (t <= tf + 25) & (fp + fm > 0)
        res.update(tau_form=float(tf), p_window=float((fp[w] / (fp[w] + fm[w])).mean()))
    else:
        res.update(tau_form=None, p_window=None)
    late = t >= t[-1] - 100
    res["late_p_plus"] = float((fp[late] / np.maximum(fp[late] + fm[late], 1e-12)).mean())
    res["late_f_core"] = float(fc[late].mean())
    res["late_mean_S"] = float(mS[late].mean())
    res["end_fp_fm"] = [float(fp[-1]), float(fm[-1])]
    # coherent oscillation frequency of <S> over the first 30 ell (zero crossings of <S> - its mean)
    e = t <= 30
    x = mS[e] - mS[e].mean()
    zc = np.where(np.sign(x[1:]) != np.sign(x[:-1]))[0]
    res["omega_coh_est"] = float(math.pi * (len(zc) - 1) / (t[e][zc[-1]] - t[e][zc[0]])) if len(zc) > 2 else None
    env = []
    for t0 in np.arange(0, t[-1], 10.0):
        ww = (t >= t0) & (t < t0 + 10)
        env.append([float(t0), float(np.abs(mS[ww]).max())])
    res["envelope"] = env
    return res


if __name__ == "__main__":
    jobs = [(1.5, s) for s in (1, 2, 3)] + [(2.0, s) for s in (1, 2, 3)] + [(4.0, s) for s in (1, 2)]
    i, nw = int(sys.argv[1]), int(sys.argv[2])
    for j, (k, s) in enumerate(jobs):
        if j % nw != i:
            continue
        p = os.path.join(OUT, "k%.1f_s%d.json" % (k, s))
        if os.path.exists(p):
            continue
        t0 = time.time()
        r = run(k, s)
        r["wall_s"] = time.time() - t0
        with open(p, "w") as fh:
            json.dump(r, fh)
        print(k, s, r["tau_form"], r["p_window"], r["late_p_plus"], r["late_f_core"], r["omega_coh_est"],
              "%.0fs" % r["wall_s"], flush=True)
