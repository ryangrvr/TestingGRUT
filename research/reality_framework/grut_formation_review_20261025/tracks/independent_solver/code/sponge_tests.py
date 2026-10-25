#!/usr/bin/env python3
"""Reflection of this track's own absorber (a smooth velocity sponge gamma(r) dU/dt with
gamma = g ((r-R_s)/D)^p, then a hard wall), for momenta p/m_X in {0.02, 0.1, 0.5} and both masses.

Two independent measurements:
  (a) continuum frequency-domain reflection |R|^2 of U'' + (w^2 - m^2 + i w gamma) U = 0 with U(L)=0
      (DOP853, rtol 1e-12), integrated from the wall to the sponge start;
  (b) time-domain free packets in the actual spectral code (exact free propagation, damping applied
      exactly in the kick), compared with the same packet propagated exactly in a box twice as large with
      no sponge.  Error = physical-region (r < R_s) error energy / initial energy, maximised over time.
Also (c) the production "big box" boundary: the same packet in a hard-wall box; error is identically zero
before the causal return time and is reported as a check.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from spectral_radial import Radial, T

SP = dict(width=3000.0, gmax=0.03, power=5.0)


def continuum_reflection(k, m, D=SP['width'], g=SP['gmax'], p=SP['power']):
    w = np.sqrt(k * k + m * m)

    def f(x, y):
        xi = min(max(x / D, 0.0), 1.0)
        return [y[1], -(w * w - m * m + 1j * w * g * xi ** p) * y[0]]
    sol = solve_ivp(f, [D, 0.0], [0j, 1 + 0j], method='DOP853', rtol=1e-12, atol=1e-30)
    U, dU = sol.y[0, -1], sol.y[1, -1]
    A = (1j * k * U + dU) / (2j * k); B = (1j * k * U - dU) / (2j * k)
    return float(abs(B / A) ** 2)


def energy(sim, row_set, rmax):
    """Physical-region free energy (4 pi weighted) of the given rows from nodal fields."""
    q, p = sim.physical()
    dq = sim.derivative(sim.Qh)[:, 1:-1]
    m = sim.r < rmax
    e = 0.0
    for i in row_set:
        e += np.sum(0.5 * p[i, m] ** 2 + 0.5 * dq[i, m] ** 2 + 0.5 * sim.m2[i] * q[i, m] ** 2)
    return 4 * np.pi * sim.h * e


def packet_test(k0, row, h, dt, center, width, Rs, tmax, nchk=40):
    rows = (0, 1) if row == 0 else (2,)
    L = Rs + SP['width']
    sim = Radial(h, L, dt, 'strang', nonlinear=False, sponge=dict(start=Rs, **SP))
    sim.free_packet(k0, width, center, row)
    ref = Radial(h, 2 * L, dt, 'strang', nonlinear=False)  # exact free propagation, no sponge
    ref.free_packet(k0, width, center, row)
    box = Radial(h, L, dt, 'strang', nonlinear=False)  # hard wall at L, no sponge
    box.free_packet(k0, width, center, row)
    E0 = energy(sim, rows, L)
    n = sim.n
    worst = 0.0; worst_box_before = 0.0; t_return = 2 * L - (center + 6 * width) - Rs
    nb = int(round(tmax / nchk / dt))
    hist = []
    for c in range(nchk):
        sim.advance_block(nb)
        # reference: exact free propagation (drift only) for the same time
        for s in (ref, box):
            s.drift(nb * dt); s.tau += nb * dt
        qs, ps = sim.physical(); qr, pr = ref.physical(); qb, pb = box.physical()
        m = sim.r < Rs
        dq_s = sim.derivative(sim.Qh)[:, 1:-1]; dq_r = ref.derivative(ref.Qh)[:, 1:-1][:, :n]
        dq_b = box.derivative(box.Qh)[:, 1:-1]
        err = 0.0; errb = 0.0
        for i in rows:
            err += np.sum(0.5 * (ps[i] - pr[i, :n])[m] ** 2 + 0.5 * (dq_s[i] - dq_r[i])[m] ** 2
                          + 0.5 * sim.m2[i] * (qs[i] - qr[i, :n])[m] ** 2)
            errb += np.sum(0.5 * (pb[i] - pr[i, :n])[m] ** 2 + 0.5 * (dq_b[i] - dq_r[i])[m] ** 2
                           + 0.5 * sim.m2[i] * (qb[i] - qr[i, :n])[m] ** 2)
        err *= 4 * np.pi * h / E0; errb *= 4 * np.pi * h / E0
        worst = max(worst, err)
        if sim.tau < t_return:
            worst_box_before = max(worst_box_before, errb)
        Ein = energy(sim, rows, Rs) / E0
        hist.append([sim.tau, err, errb, Ein, sim.Erec / E0])
    # energy bookkeeping for the sponge run
    Etot = energy(sim, rows, L) / E0
    return dict(k=k0, mass=float(np.sqrt(sim.m2[rows[0]])), h=h, dt=dt, center=center, width=width,
                sponge_start=Rs, L=L, tmax=tmax, max_physical_error_energy_over_initial=worst,
                final_energy_left_in_physical_region=hist[-1][3],
                sponge_energy_receiver_over_initial=sim.Erec / E0,
                field_plus_receiver_minus_one=Etot + sim.Erec / E0 - 1.0,
                bigbox_max_error_before_causal_return=worst_box_before, causal_return_time=t_return,
                history=hist)


def main():
    out = {'sponge': SP, 'continuum': [], 'packets': []}
    for m in ((1.0, np.sqrt(6.0)) if not (len(sys.argv) > 2 and sys.argv[2] == 'soft') else ()):
        for k in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0):
            out['continuum'].append(dict(k=k, mass=m, reflected_energy=continuum_reflection(k, m)))
            print('continuum', k, m, out['continuum'][-1]['reflected_energy'], flush=True)
    cases = [  # k, row, h, dt, center, width, Rs, tmax
        (0.5, 0, 0.5, 0.125, 300.0, 40.0, 600.0, 16000.0),
        (0.1, 0, 1.0, 0.25, 300.0, 60.0, 600.0, 70000.0),
        (0.02, 0, 2.0, 0.5, 600.0, 250.0, 1200.0, 200000.0),
        (0.5, 2, 0.5, 0.125, 300.0, 40.0, 600.0, 35000.0),
        (0.1, 2, 1.0, 0.25, 300.0, 60.0, 600.0, 170000.0),
        (0.02, 2, 2.0, 0.25, 600.0, 250.0, 1200.0, 500000.0),
    ]
    if len(sys.argv) > 2 and sys.argv[2] == 'soft':
        # low-momentum cases with wide, odd-symmetrised packets far from the origin
        out['continuum'] = []
        cases = [(0.02, 0, 2.0, 0.5, 1200.0, 400.0, 2400.0, 230000.0),
                 (0.02, 2, 2.0, 0.25, 1200.0, 400.0, 2400.0, 520000.0)]
    for c in cases:
        t0 = time.time()
        r = packet_test(*c)
        r['seconds'] = time.time() - t0
        out['packets'].append(r)
        print({k: v for k, v in r.items() if k != 'history'}, flush=True)
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('sponge_tests.json')
    dest.write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
