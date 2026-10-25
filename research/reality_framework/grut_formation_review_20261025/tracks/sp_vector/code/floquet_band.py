#!/usr/bin/env python3
"""Track sp_vector (own code): fine search of the first transverse-vector resonance band in the late X orbit
(k near omega_radial/2), to resolve narrow bands missed by the coarse grid of floquet_vector.py.
Usage: python3 floquet_band.py ../results/floquet_vector.json"""
import json, sys
import numpy as np
from scipy.optimize import minimize_scalar
import floquet_vector as fv

def main(path):
    sol, T, rmax, rmin = fv.orbit()
    omega = 2 * np.pi / T
    ts = np.linspace(0, T, 4001)
    r2 = np.array([sol.sol(t)[0]**2 + sol.sol(t)[1]**2 for t in ts])
    out = {}
    for e in [1e-7, 1e-5, 1e-4, 1e-3, 3e-3, 5e-3]:
        m2bar = 2 * e * e * r2.mean()
        kc = np.sqrt(max((omega / 2)**2 - m2bar, 1e-12))
        best = (0.0, kc)
        for n in (1, 2):  # first two bands, k ~ n omega/2
            kcn = np.sqrt(max((n * omega / 2)**2 - m2bar, 1e-12))
            ks = np.linspace(kcn * (1 - 2e-2), kcn * (1 + 2e-2), 801)
            mus = np.array([fv.floquet(sol, T, e, k) for k in ks])
            i = int(np.argmax(mus))
            if mus[i] > best[0]:
                best = (float(mus[i]), float(ks[i]))
        # analytic narrow-resonance estimate for the first band: mu ~ |c_1| / (4 k), c_1 = Fourier amplitude of m^2 at omega
        m2 = 2 * e * e * r2
        c1 = 2 * abs(np.mean(m2 * np.exp(-1j * omega * ts[:len(m2)]))) 
        out[str(e)] = dict(mu_max=best[0], k=best[1], mu_narrow_estimate=float(c1 / (4 * kc)),
                           growth_2mu_x400=2 * best[0] * 400)
        print(e, out[str(e)], flush=True)
    d = json.load(open(path)); d['band_search'] = out; d['omega_radial'] = omega
    json.dump(d, open(path, 'w'), indent=1)

if __name__ == '__main__':
    main(sys.argv[1])
