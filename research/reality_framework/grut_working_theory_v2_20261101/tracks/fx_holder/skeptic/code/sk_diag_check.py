#!/usr/bin/env python3
"""fx_holder_skeptic (Claude; own code): the forced-Y point used by the free-Y diagnostic (NYL/NYH) of rpholder.c and
rptri.c is yp = (-A J + B conj(J))/D with B = 2 mu_eff Sigma.  The Y equation is dV/dY* = M^2 Y + kap X^3 + 2 conj(u) S conj(Y)
+ ..., so the static forced solution of M^2 Y + c conj(Y) = -kap X^3 has c = 2 conj(u) S.  The two agree only for real u.
This evaluates both at the prepared point (u = 10 i, lambda = 1) and at the late vacuum value u = mu.  Output: printed +
../results/sk_diag_check.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_common as sc

z = np.array([87.19556213287578, -5.75424796328559, -28.387181944361622j])      # sk_point.py lam = 1 prepared point
out = {}
for nm, u in (('prep_u_10i', 10j), ('late_u_mu', sc.MU + 0j)):
    J = sc.KAP * z[0] ** 3; A = sc.M ** 2
    for lab, c in (('track_B_2uS', 2 * u * z[2]), ('correct_c_2conj(u)S', 2 * np.conj(u) * z[2])):
        yp = (-A * J + c * np.conj(J)) / (A * A - abs(c) ** 2)
        out[f'{nm}_{lab}'] = dict(Y_forced=[yp.real, yp.imag], actual_Y=[z[1].real, z[1].imag], displacement=float(abs(z[1] - yp)),
                                  rotation_phase=float(0.5 * np.angle(c)))
print(json.dumps(out, indent=1))
json.dump(out, open(Path(__file__).resolve().parent.parent / 'results' / 'sk_diag_check.json', 'w'), indent=1)
