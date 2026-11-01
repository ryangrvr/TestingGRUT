"""Skeptic part F: slaved-mass control (X mass g M(s) imposed, phi absent) in the independent DOP853 code."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sk_common import *
out = []
for orb in ['tau20', 'late']:
    y4 = tau20_state(); q = charge(y4)
    if orb == 'late': y4 = np.array([31.494634874218747, 0.0, 0.0, q / 31.494634874218747])
    info = orbit_info(y4, nsamp=1024)
    Mstart = 0.5 * np.sqrt(G2 / LP) * info['rmin']
    for eps in [1e-3, 1e-4]:
        devs = []
        for ph in np.arange(8) / 8:
            st = phase_state(y4, info['T_r'], ph)
            m2 = G2 * Mstart ** 2
            J0 = radial_quad(energy4(st, LX, m2), charge(st), LX, m2)['J_r'] / charge(st)
            r = run_measure(st, info, eps, 0.0, mode=2, rtol=1e-12)
            devs.append(r['nminus'] - J0)
        out.append(dict(orbit=orb, eps=eps, max_abs_dev=float(np.max(np.abs(devs))), mean_dev=float(np.mean(devs))))
        print(out[-1])
save_json('skF_controls.json', dict(slaved_mass=out))
