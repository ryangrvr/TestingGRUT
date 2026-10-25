#!/usr/bin/env python3
"""NEW CONTROL using Codex's OLD executable (bb211d1 radial_nodal.cpp, unmodified, built with g++ -O2):
the same nodal Gaussian preparation as refinement.py, but with the ordinary cubic damper moved
from R=270 (L=360) to R=720 (L=810).  Purpose: decide whether the 0.6% gap between the old finest
retention (N00625, 35.153%) and the new CFS finest retention (P00625, 35.363%) comes from the old
boundary.  Initial data are rebuilt here with the same formulas as refinement.py (own code).
Usage: python3 old_far_absorber.py RADIAL_NODAL_BINARY OUTDIR NAME h dt L absorber [tmax]
"""
import hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp


def main():
    binary, outdir, name = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
    h, dt, L, ab = map(float, sys.argv[4:8])
    tmax = float(sys.argv[8]) if len(sys.argv) > 8 else 12000.
    amp, width, probe, eta = .18, 60., 120., .03
    outdir.mkdir(parents=True, exist_ok=True)
    f = lambda r: amp * np.exp(-r * r / (2 * width * width))
    rhs = lambda r, y: np.array([y[1], (3 * (y[0] ** 2 - 1) + f(r) ** 2) * y[0]])
    x = np.linspace(0, L, 601)
    guess = np.array([np.sqrt(1 - f(x) ** 2 / 3), f(x) ** 2 * x / (3 * width * width)])
    s = solve_bvp(rhs, lambda l, r: [l[1], r[0] - 1], x, guess, S=np.diag([0., -2.]), tol=1e-10, bc_tol=1e-11, max_nodes=12000)
    assert s.success
    r = np.arange(1, round(L / h)) * h; F = f(r); z = np.zeros_like(r)
    init = outdir / f'{name}_initial.csv'
    np.savetxt(init, np.column_stack([r, F, z, s.sol(r)[0], z, -F, z]), delimiter=',', header='rho,x,y,s,px,py,ps', comments='')
    prefix = outdir / name
    t0 = time.monotonic()
    proc = subprocess.run([binary, str(init), str(prefix), str(h), str(dt), str(L), str(tmax), str(probe), str(ab), str(eta)],
                          check=True, capture_output=True, text=True)
    hist = np.loadtxt(str(prefix) + '_history.csv', delimiter=',', skiprows=1)
    e0, q0 = hist[0, 1:3]; late = hist[:, 0] > .8 * tmax
    out = dict(label='NEW CONTROL with Codex old executable (radial_nodal.cpp), own driver', name=name,
               parameters=dict(h=h, dt=dt, length=L, absorber=ab, eta=eta, tmax=tmax, probe=probe, amplitude=amp, width=width),
               seconds=time.monotonic() - t0, Qbar_initial=float(q0), Ebar_initial=float(e0),
               max_energy_error=float(np.max(abs((hist[:, 1] + hist[:, 3]) / e0 - 1))),
               max_charge_error=float(np.max(abs((hist[:, 2] + hist[:, 4]) / q0 - 1))),
               core60_late_mean_fraction=float(np.mean(hist[late, 12]) / q0),
               core60_at=dict((str(t), float(np.interp(t, hist[:, 0], hist[:, 12]) / q0)) for t in [2000., 4000., 6000., 8000., 10000., 12000.]),
               Qbar_absorbed_final=float(hist[-1, 4]),
               history_sha256=hashlib.sha256(hist.tobytes()).hexdigest())
    (outdir / f'{name}_RESULTS.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'parameters'}))


if __name__ == '__main__':
    main()
