#!/usr/bin/env python3
"""INDEPENDENT audit (own code): which part of each ff49e5d free packet could have been reflected
back into the measured physical region by the snapshot time tmax?

The packet tests compare a single snapshot at tmax with a large-domain reference.  A mode whose
wall/layer return has not yet re-entered r < R-4h at tmax cannot contribute to that error, so the
snapshot error does not test the layer at that momentum.  Group velocity is the lattice value
v = p p'(k)/omega; inside the CFS layer the real coordinate stretch Re Gamma = 1 + sigma alpha/(alpha^2+omega^2)
lengthens the round trip.  Spectrum: the Gaussian packet's analytic transform, energy weight ~ omega^2 |a(q)|^2.
Usage: python3 indep_packet_coverage.py OUT.json
"""
import json, sys
import numpy as np

BASE = dict(h=.25, dt=.03125, length=450., probe=240., absorber=360., sigma=20., alpha=.1, mass=1.)
CASES = {
    'free_X': dict(h=.5, tmax=6000., k=.1, packet_width=40., center=120.),
    'free_X_time': dict(h=.5, tmax=6000., k=.1, packet_width=40., center=120.),
    'free_s': dict(h=.5, tmax=6500., mass=np.sqrt(6.), k=.2, packet_width=30., center=120.),
    'free_soft': dict(h=2., length=1090., absorber=1000., probe=800., tmax=60000., k=.02, packet_width=100., center=300.),
}


def main(out):
    rows = []
    for name, extra in CASES.items():
        c = {**BASE, **extra}
        h, m, W = c['h'], c['mass'], c['length'] - c['absorber']
        q = np.linspace(-np.pi / h, np.pi / h, 400001)[1:-1]
        s = np.sin(q * h / 2); p = 2 * s / h * (1 + s * s / 6); om = np.sqrt(m * m + p * p)
        v = np.abs(p * np.cos(q * h / 2) * (1 + .5 * s * s) / om)
        amp = np.exp(-(q - c['k']) ** 2 * c['packet_width'] ** 2 / 2)
        wgt = om ** 2 * amp ** 2; wgt /= wgt.sum()
        x = np.linspace(0, W, 2001)
        stretch = np.array([np.trapezoid(1 + c['sigma'] * (x / W) ** 3 * c['alpha'] / (c['alpha'] ** 2 + w * w), x) for w in om[::1000]])
        stretch = np.interp(om, om[::1000][np.argsort(om[::1000])], stretch[np.argsort(om[::1000])])
        R = c['absorber']; edge = R - 4 * h; c0 = c['center']
        path = np.where(q > 0, R - c0, R + c0) + 2 * stretch + (R - edge)
        tret = path / np.maximum(v, 1e-300)
        covered = tret < c['tmax']
        # momentum below which the round trip is incomplete (q>0 branch)
        qpos = q[(q > 0) & ~covered]
        kmin = float(qpos.max()) if len(qpos) else 0.
        rows.append(dict(name=name, k0=c['k'], sigma_k=1 / c['packet_width'], tmax=c['tmax'],
                         energy_fraction_whose_return_can_reenter_by_tmax=float(wgt[covered].sum()),
                         energy_fraction_not_tested=float(wgt[~covered].sum()),
                         smallest_tested_k_on_outgoing_branch=float(q[(q > 0) & covered].min()),
                         energy_fraction_with_k_below_0p01=float(wgt[np.abs(q) < .01].sum())))
        print(rows[-1])
    json.dump(dict(label='INDEPENDENT (own code) snapshot-coverage audit of the free-packet tests', rows=rows), open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
