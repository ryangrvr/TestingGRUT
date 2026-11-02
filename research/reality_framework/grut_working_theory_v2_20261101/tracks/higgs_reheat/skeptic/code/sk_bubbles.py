"""Thermal bubbles in our past light cone along the engine history, re-done with the canonical S3/T of sk_thermal.py
(same counting as the track's thermal.py: N_b = e^{3N} int dtau a^3 Gamma/H_i^4, Gamma = T^4 (S3/2piT)^{3/2} e^{-S3/T}).
Writes ../results/sk_bubbles.json."""
import json
import os
import sys

import numpy as np
from scipy.interpolate import interp1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'code')))
from background import T_GeV  # noqa: E402

bg = json.load(open(os.path.join(HERE, '..', '..', 'results', 'background.json')))['table_downsampled']
tau, a, rr = np.array(bg['tau']), np.array(bg['a']), np.array(bg['rho_r'])
th = json.load(open(os.path.join(HERE, '..', 'results', 'sk_thermal.json')))
out = {}
for Mt in ('172.87', '172.6'):
    rows = th[f'Mt={Mt}']
    Ts = np.array([float(k) for k in rows])
    res = {}
    for norm in ('track', 'canon'):
        Ss = np.array([rows[k][norm]['S3_over_T'] for k in rows])
        Sint = interp1d(np.log(Ts), Ss, bounds_error=False, fill_value=(Ss[0], Ss[-1]))
        for H_GeV in (1.29e11, 1.51e11):
            T = T_GeV(rr, H_GeV)
            m = T > 1e12
            S3 = Sint(np.log(T[m]))
            lnG = 4 * np.log(T[m] / H_GeV) + 1.5 * np.log(S3 / (2 * np.pi)) - S3
            integ = np.log(np.trapezoid(np.exp(lnG - lnG.max()) * a[m] ** 3, tau[m])) + lnG.max()
            res[f'{norm} H={H_GeV:.3g}'] = dict(log10_Nb_N60=float((integ + 180) / np.log(10)),
                                                min_S3_over_T=float(S3.min()))
    out[f'Mt={Mt}'] = res
    print(Mt, res)
json.dump(out, open(os.path.join(HERE, '..', 'results', 'sk_bubbles.json'), 'w'), indent=1)
