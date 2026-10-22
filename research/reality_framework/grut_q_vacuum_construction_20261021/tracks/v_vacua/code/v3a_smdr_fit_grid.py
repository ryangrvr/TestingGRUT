"""SMDR full-precision MSbar fits (pole masses -> couplings at Q = Mt) on a grid of (Mt, Mh, alpha_s).

Verifier v_vacua. Each SMDR fit takes ~25 s, so the grid is computed once (4 processes) and cached in
../results/v3a_smdr_fit_grid.json; v3b_sm_check.py interpolates in it.  Run: python3 v3a_smdr_fit_grid.py
"""
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import smdr_common as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'v3a_smdr_fit_grid.json')

GRID = [
    # Mt scan at the PDG 2025 point (Mh = 125.20, alpha_s = 0.1180)
    (170.4, 125.20, 0.1180), (170.8, 125.20, 0.1180), (171.2, 125.20, 0.1180), (171.6, 125.20, 0.1180),
    (172.0, 125.20, 0.1180), (172.56, 125.20, 0.1180),
    # Mh and alpha_s derivatives at Mt = 171.2
    (171.2, 124.70, 0.1180), (171.2, 125.70, 0.1180), (171.2, 125.20, 0.1171), (171.2, 125.20, 0.1189),
    # Buttazzo et al. reference point and Bednyakov et al. inputs
    (173.34, 125.15, 0.1184), (171.4, 125.70, 0.1185),
    # far points for the Froggatt-Nielsen double condition (Mh extrapolation)
    (172.6, 128.60, 0.1180), (174.2, 132.00, 0.1180), (171.7, 126.80, 0.1180), (171.2, 128.6, 0.1180),
]


def one(p):
    t = time.time()
    r = S.fit(*p)
    r['inputs'] = list(p)
    r['seconds'] = time.time() - t
    return r


if __name__ == '__main__':
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(one, GRID))
    with open(OUT, 'w') as f:
        json.dump({'source': 'SMDR v1.3 SMDR_Fit_Inputs (Martin-Robertson 1907.02500), PDG-2025 reference values for '
                             'all other inputs (GF, MZ, alpha, light masses, Delta alpha_had)',
                   'fits': res}, f, indent=1)
    print(json.dumps(res, indent=1))
