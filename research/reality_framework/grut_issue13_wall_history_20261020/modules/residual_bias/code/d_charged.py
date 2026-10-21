"""Route D: charge forced into a selected +v vacuum above the spinodal n_s1 (outside-action supply [P]).
Small 2D classical lattice ensemble (independent implementation; same equations and preparation as
charged_exit_sim's C2 runs, which are quoted for comparison) to estimate the ensemble bias m0 = 2p(+) - 1 at
vacuum formation as a function of n_i/n_c and of the expansion rate, including the Minkowski (H -> 0) limit that
is physically relevant (physical H/m_X ~ 5e-10, about 1e7 times below the simulated values).

Units (charge_phase / charged_exit_sim conventions): S = phi/v, Y = X/v, lengths and times in ell = 1/(100 GeV),
Ubar = (S^2 - 1)^2/4 + A_G S^2 |Y|^2 + B_X |Y|^4, A_G = 1, B_X = 0.1.  Conformal variables chi = a f in RD
(a = tau/tau_i, a'' = 0) or a = 1 (Minkowski):
    chi_S'' = lap chi_S - chi_S (chi_S^2 - a^2 + 2 A_G |chi_Y|^2)
    chi_Y'' = lap chi_Y - chi_Y (A_G chi_S^2 + 2 B_X |chi_Y|^2)
Leapfrog (kick-drift-kick), periodic 5-point Laplacian, numpy only, deterministic seeds.
Preparation (C2, 'memory_forced'): S = +1 + noise (std 1e-2, Gaussian-smoothed, s_c = 1), X on its own circular
orbit in the fixed +v vacuum at n_i = kappa n_c (omega^2 = A_G + 2 B_X |Y|^2, n = 2 omega |Y|^2), X noise 1e-3.
Formation: first output after restoration (core fraction >= 1/2 seen) with both vacuum fractions >= 5%;
p(+) = f+/(f+ + f-).  If only one sign appears (other < 5% until the end) the outcome is recorded as 'one sign'.
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

A_G, B_X = 1.0, 0.1
NBAR_C = math.sqrt(8 * B_X) * (1.0 / (2 * A_G)) ** 1.5           # 0.316228
NBAR_S1 = 0.427667                                                # [Q charge_phase]
NBAR_STAR = (1.0 / B_X) ** 0.25                                   # 1.778
HERE = os.path.dirname(os.path.abspath(__file__))
RUNDIR = os.path.join(os.path.dirname(HERE), "results", "d_runs")


def fixed_vacuum_state(nbar, S0=1.0):
    f = lambda y2: 2 * math.sqrt(A_G * S0 ** 2 + 2 * B_X * y2) * y2 - nbar
    y2 = brentq(f, 0.0, 100.0)
    return math.sqrt(y2), math.sqrt(A_G * S0 ** 2 + 2 * B_X * y2)


def gaussian_field(N, L, s_c, amp, rng):
    white = rng.standard_normal((N, N))
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    k2 = k1[:, None] ** 2 + k1[None, :] ** 2
    f = np.real(np.fft.ifft2(np.fft.fft2(white) * np.exp(-0.5 * k2 * s_c ** 2)))
    f -= f.mean()
    return f * (amp / f.std())


def lap(f, dx):
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / dx ** 2


def run(kappa, mode, seed, N=384, dx=0.5, tau_i=150.0, a_final=2.5, tau_max=300.0, dt=0.1, out_every=0.5):
    """mode: 'expand' (RD from a = 1 at tau = tau_i) or 'mink' (a = 1 fixed, n fixed)."""
    rng = np.random.default_rng(1000 * seed + int(round(kappa * 100)))
    L = N * dx
    nbar = kappa * NBAR_C
    Y0, om = fixed_vacuum_state(nbar)
    expand = mode == "expand"
    radial = -0.5 / tau_i if expand else 0.0
    sdot = 1.0 / tau_i if expand else 0.0
    cS = 1.0 + gaussian_field(N, L, 1.0, 1e-2, rng)
    pS = np.full((N, N), sdot)
    cR = Y0 + gaussian_field(N, L, 1.0, 1e-3, rng)
    cI = gaussian_field(N, L, 1.0, 1e-3, rng)
    pR = radial * cR - om * cI
    pI = radial * cI + om * cR
    tau = tau_i if expand else 0.0
    tau_end = tau_i * a_final if expand else tau_max
    a_of = (lambda t: t / tau_i) if expand else (lambda t: 1.0)

    def force(a):
        y2 = cR * cR + cI * cI
        s2 = cS * cS
        fS = lap(cS, dx) - cS * (s2 - a * a + 2 * A_G * y2)
        m = A_G * s2 + 2 * B_X * y2
        return fS, lap(cR, dx) - cR * m, lap(cI, dx) - cI * m

    ser = []
    fS, fR, fI = force(a_of(tau))
    nsteps = int(round((tau_end - tau) / dt))
    every = int(round(out_every / dt))
    restored = False
    form = None
    t0 = time.time()
    for it in range(nsteps + 1):
        if it % every == 0:
            a = a_of(tau)
            S = cS / a
            fp = float((S > 0.5).mean())
            fm = float((S < -0.5).mean())
            fc = 1.0 - fp - fm
            n_phys = float((2 * (cR * pI - cI * pR)).mean() / a ** 3)
            ser.append((round(tau, 4), a, float(S.mean()), float(np.abs(S).mean()), fp, fm, fc, n_phys))
            if fc >= 0.5:
                restored = True
            if restored and form is None and fp >= 0.05 and fm >= 0.05:
                form = dict(tau=tau, a=a, p_plus=fp / (fp + fm), f_plus=fp, f_minus=fm, n_over_nc=n_phys / NBAR_C)
        if it == nsteps:
            break
        # KDK
        pS += 0.5 * dt * fS
        pR += 0.5 * dt * fR
        pI += 0.5 * dt * fI
        cS += dt * pS
        cR += dt * pR
        cI += dt * pI
        tau += dt
        fS, fR, fI = force(a_of(tau))
        pS += 0.5 * dt * fS
        pR += 0.5 * dt * fR
        pI += 0.5 * dt * fI
    last = ser[-1]
    end = dict(tau=last[0], f_plus=last[4], f_minus=last[5],
               p_plus=(last[4] / (last[4] + last[5]) if last[4] + last[5] > 0 else None))
    outcome = "two signs" if form is not None else ("one sign (+)" if last[4] > 0.05 and last[5] < 0.05 else
                                                     ("one sign (-)" if last[5] > 0.05 and last[4] < 0.05 else
                                                      "no vacuum formed"))
    return dict(kappa=kappa, n_over_ns1=kappa * NBAR_C / NBAR_S1, n_over_nstar=kappa * NBAR_C / NBAR_STAR,
                mode=mode, seed=seed, N=N, dx=dx, tau_i=tau_i if expand else None, Y0=Y0, omega=om,
                restored=restored, formation=form, end=end, outcome=outcome, wall_s=time.time() - t0,
                series_cols=["tau", "a", "mean_S", "mean_abs_S", "f_plus", "f_minus", "f_core", "n_phys_bar"],
                series=ser)


JOBS = []
for s in (1, 2, 3):
    JOBS += [("expand", 2.0, s, 150.0), ("expand", 4.0, s, 150.0)]
for s in (1, 2, 3):
    JOBS += [("expand", 2.0, s, 600.0), ("expand", 4.0, s, 600.0)]
for s in (1, 2, 3, 4):
    for k in (1.5, 2.0, 3.0, 4.0, 5.0, 7.0):
        JOBS.append(("mink", k, s, None))


def job_name(j):
    mode, k, s, ti = j
    return "%s_k%.1f_s%d%s" % (mode, k, s, ("_t%d" % ti) if ti else "")


def worker(idx, nw):
    os.makedirs(RUNDIR, exist_ok=True)
    for i, j in enumerate(JOBS):
        if i % nw != idx:
            continue
        name = job_name(j)
        path = os.path.join(RUNDIR, name + ".json")
        if os.path.exists(path):
            continue
        mode, k, s, ti = j
        kw = dict(tau_i=ti) if ti else {}
        if mode == "expand" and ti == 600.0:
            kw["a_final"] = 1.8
        r = run(k, mode, s, **kw)
        with open(path, "w") as fh:
            json.dump(r, fh)
        print(name, r["outcome"], r["formation"], "%.0fs" % r["wall_s"], flush=True)


def collect():
    rows = []
    for j in JOBS:
        path = os.path.join(RUNDIR, job_name(j) + ".json")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            r = json.load(fh)
        # coherent envelope: max |mean S| in windows of 10 time units
        ser = np.array(r["series"])
        tau = ser[:, 0]
        env = []
        for t0 in np.arange(tau[0], tau[-1], 10.0):
            w = (tau >= t0) & (tau < t0 + 10.0)
            if w.any():
                env.append((float(t0), float(np.abs(ser[w, 2]).max()), float(ser[w, 3].mean())))
        r2 = {k: v for k, v in r.items() if k not in ("series",)}
        r2["coherent_envelope_[tau0,max|meanS|,mean|S|]"] = env
        r2["m0_form"] = (2 * r["formation"]["p_plus"] - 1) if r["formation"] else None
        rows.append(r2)
    return rows


def validate(N=128, dx=0.5, dt=0.1, tau=60.0, kappa=4.0):
    """Minkowski run: total energy and charge drift (leapfrog, should be O(dt^2) bounded and round-off)."""
    rng = np.random.default_rng(7)
    L = N * dx
    Y0, om = fixed_vacuum_state(kappa * NBAR_C)
    cS = 1.0 + gaussian_field(N, L, 1.0, 1e-2, rng)
    pS = np.zeros((N, N))
    cR = Y0 + gaussian_field(N, L, 1.0, 1e-3, rng)
    cI = gaussian_field(N, L, 1.0, 1e-3, rng)
    pR, pI = -om * cI, om * cR

    def grad2(f):
        return ((np.roll(f, -1, 0) - f) ** 2 + (np.roll(f, -1, 1) - f) ** 2) / dx ** 2

    def energy():
        y2 = cR ** 2 + cI ** 2
        U = (cS ** 2 - 1) ** 2 / 4 + A_G * cS ** 2 * y2 + B_X * y2 ** 2
        return float((0.5 * pS ** 2 + 0.5 * grad2(cS) + pR ** 2 + pI ** 2 + grad2(cR) + grad2(cI) + U).sum() * dx * dx)

    def charge():
        return float((2 * (cR * pI - cI * pR)).sum() * dx * dx)

    def force():
        y2 = cR * cR + cI * cI
        s2 = cS * cS
        m = A_G * s2 + 2 * B_X * y2
        return lap(cS, dx) - cS * (s2 - 1 + 2 * A_G * y2), lap(cR, dx) - cR * m, lap(cI, dx) - cI * m
    # kinetic term |dY|^2 = dR^2 + dI^2: Y'' = lap Y - dU/dY*, dU/dY* = (A S^2 + 2 B |Y|^2) Y; conserved
    # H = sum [pS^2/2 + |grad S|^2/2 + pR^2 + pI^2 + |grad R|^2 + |grad I|^2 + U], Q = sum 2 (R pI - I pR).
    E0, Q0 = energy(), charge()
    fS, fR, fI = force()
    emax = 0.0
    for it in range(int(tau / dt)):
        pS += 0.5 * dt * fS
        pR += 0.5 * dt * fR
        pI += 0.5 * dt * fI
        cS += dt * pS
        cR += dt * pR
        cI += dt * pI
        fS, fR, fI = force()
        pS += 0.5 * dt * fS
        pR += 0.5 * dt * fR
        pI += 0.5 * dt * fI
        if it % 20 == 0:
            emax = max(emax, abs(energy() / E0 - 1))
    return dict(N=N, dx=dx, dt=dt, tau=tau, kappa=kappa, max_rel_energy_drift=emax,
                rel_charge_drift=abs(charge() / Q0 - 1))


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        print(json.dumps(validate()))
    elif sys.argv[1] == "worker":
        worker(int(sys.argv[2]), int(sys.argv[3]))
    elif sys.argv[1] == "time":
        r = run(4.0, "mink", 1, tau_max=20.0)
        print(r["wall_s"], r["series"][-1])
