#!/usr/bin/env python3
"""INDEPENDENT check (own code, own method) of the old ordinary-damper reflection (9636950).

Model (taken from the old radial_nodal.cpp, not from Codex's boundary.py): the momentum of every
field is multiplied by exp(-dt*eta*z^3/2) before and after each Yoshida step, z=(r-ab)/(L-ab),
with a Dirichlet wall at r=L.  For a time-harmonic charged wave e^{-i omega t} in the linear
exterior (S=1, |X| small) this is U'' = -(k^2 + i omega eta z^3) U, U(D)=0, omega^2 = 1 + k^2.

Method (different from Codex's DOP853 shooting and its 2nd-order FD BVP): exact piecewise-constant
transfer matrices on N slabs with midpoint eta, then Richardson extrapolation in N (error ~ N^-2).
The returned |R|^2 is the energy reflection of a unit incident wave at the layer entrance.
Usage: python3 indep_old_damper.py OUT.json
"""
import json, sys
import numpy as np


def reflection(k, width, eta, nslab, power=3):
    omega = np.hypot(1.0, k)
    d = width / nslab
    # start at wall: U=0, U'=1; march toward the entrance x=0
    U, P = 0j, 1 + 0j
    for j in range(nslab - 1, -1, -1):
        xm = (j + 0.5) * d
        q = np.sqrt(k * k + 1j * omega * eta * (xm / width) ** power + 0j)
        c, s = np.cos(q * d), np.sin(q * d)
        # moving from x to x-d:  U(x-d)=c U - s/q U',  U'(x-d)= q s U + c U'
        U, P = c * U - s / q * P, q * s * U + c * P
    inc = (U + P / (1j * k)) / 2
    ref = (U - P / (1j * k)) / 2
    return ref / inc


def converged(k, width, eta):
    n = 4000
    r1 = reflection(k, width, eta, n)
    r2 = reflection(k, width, eta, 2 * n)
    r3 = reflection(k, width, eta, 4 * n)
    rich = (4 * r3 - r2) / 3
    return rich, abs(rich - (4 * r2 - r1) / 3)


def main(out):
    cases = [('original', 90., .03), ('half_strength', 90., .015), ('old_larger_box', 120., .03),
             ('candidate_wide', 360., .1), ('candidate_very_wide', 720., .1)]
    ks = [.01, .02, .03, .05, .1, .2, .5, 1.]
    rows = []
    for name, width, eta in cases:
        for k in ks:
            r, err = converged(k, width, eta)
            rows.append(dict(case=name, width=width, eta=eta, k_over_m=k, reflected_energy=float(abs(r) ** 2),
                             richardson_amplitude_uncertainty=float(err)))
    # zero-damping control: exact -exp(2ikD)
    ctrl = []
    for k in [.01, .1, 1.]:
        r, _ = converged(k, 90., 0.)
        ctrl.append(dict(k=k, zero_damping_error=float(abs(r + np.exp(2j * k * 90.)))))
    json.dump(dict(label='INDEPENDENT (own transfer-matrix code)', model='U\'\'=-(k^2+i omega eta (x/D)^3)U, U(D)=0',
                   rows=rows, zero_damping_controls=ctrl), open(out, 'w'), indent=1)
    for r in rows:
        if r['case'] in ('original', 'old_larger_box', 'candidate_very_wide'):
            print(r['case'], r['k_over_m'], '%.6g%%' % (100 * r['reflected_energy']))
    print(ctrl)


if __name__ == '__main__':
    main(sys.argv[1])
