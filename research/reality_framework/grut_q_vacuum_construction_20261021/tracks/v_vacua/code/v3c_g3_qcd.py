"""Verifier v_vacua, part 3c: g3(mu = Mt) from alpha_s^(5)(MZ) by pure-QCD running (1-4 loops, nf = 5).

Checks the Buttazzo et al. fit g3(Mt) = 1.1666 + ... (1307.3536 sec. 3; their TeX carries a commented-out 1.1645)
used by track `vacua`, against an independent pure-QCD evaluation, SMDR's full fit and Bednyakov et al. 1507.08833
(g_s = 1.16540 at alpha_s = 0.1185, Mt = 173.21). Top decoupling at mu = Mt changes g3 by O(1e-4) only (range shown).
Beta coefficients: van Ritbergen-Vermaseren-Larin / Czakon [R]. Run: python3 v3c_g3_qcd.py -> ../results/v3c_g3_qcd.json
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
Z3 = 1.2020569031595943


def beta(nf):
    b0 = 11 - 2 * nf / 3
    b1 = 102 - 38 * nf / 3
    b2 = 2857 / 2 - 5033 * nf / 18 + 325 * nf ** 2 / 54
    b3 = (149753 / 6 + 3564 * Z3) - (1078361 / 162 + 6508 / 27 * Z3) * nf + (50065 / 162 + 6472 / 81 * Z3) * nf ** 2 \
        + 1093 / 729 * nf ** 3
    return np.array([b0 / 4, b1 / 16, b2 / 64, b3 / 256])     # da/dln mu^2 = -sum b_i a^{i+2}, a = alpha_s/pi


def run(als, mu0, mu1, nf=5, L=4):
    b = beta(nf)[:L]
    f = lambda t, a: [-sum(b[i] * a[0] ** (i + 2) for i in range(L))]
    s = solve_ivp(f, (0, np.log(mu1 ** 2 / mu0 ** 2)), [als / np.pi], rtol=1e-12, atol=1e-15)
    return s.y[0, -1]


out = {}
for als, Mt in ((0.1184, 173.34), (0.1180, 172.56), (0.1185, 173.21)):
    row = {}
    for L in (1, 2, 3, 4):
        a = run(als, 91.1876, Mt, 5, L)
        row[f'{L}-loop'] = float(np.sqrt(4 * np.pi ** 2 * a))
    a = run(als, 91.1876, Mt, 5, 4)
    row['4-loop, decoupling band (c2 in [-11/72, 7/24])'] = [float(np.sqrt(4 * np.pi ** 2 * a * (1 + c * a * a)))
                                                            for c in (-11 / 72, 7 / 24)]
    row['Buttazzo fit'] = 1.1666 + 0.00314 * (als - 0.1184) / 0.0007 - 0.00046 * (Mt - 173.34)
    out[f'alpha_s={als}, Mt={Mt}'] = row
out['SMDR fit (173.34, 125.15, 0.1184)'] = 1.164856
out['SMDR fit (172.56, 125.20, 0.1180)'] = 1.163427
out['Bednyakov et al. 1507.08833 (0.1185, 173.21)'] = 1.16540
json.dump(out, open(os.path.join(HERE, '..', 'results', 'v3c_g3_qcd.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
