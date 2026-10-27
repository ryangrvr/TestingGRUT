"""SMDR full-precision MSbar fits (pole masses -> couplings at Q = Mt) on a grid around the PDG 2026 point.

Track cf_vacuum (closure-forcing pass, 27 Oct 2026). Reuses the SMDR driver and wrapper written by verifier v_vacua
(q-theory construction, tracks/v_vacua/code/smdr_driver.c, smdr_common.py; copied here unchanged) and the SMDR v1.3
library of Martin and Robertson (arXiv:1907.02500), already built in the session scratchpad by v_vacua.

PDG 2026 inputs (pdglive, 2026 edition, read 27 Oct 2026):
  M_h = 125.13 +- 0.11 GeV (S = 1.5)             [S126M]
  M_t = 172.60 +- 0.27 GeV (direct, S = 1.5)     [Q007TP]
  M_t = 172.1  +- 0.6  GeV (pole from xsec)      [Q007TP4]
  alpha_s(M_Z) = 0.1180 +- 0.0009                 [QCD review rev. Aug 2025, eq. (9.25)]
Each fit takes ~28 s; run once (4 processes), cached in ../results/c1_fit_grid.json.
"""
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import smdr_common as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'c1_fit_grid.json')
MH0, AS0 = 125.13, 0.1180

GRID = [(mt, MH0, AS0) for mt in (170.6, 171.0, 171.4, 171.8, 172.1, 172.6, 173.0)]
GRID += [(171.0, MH0 - 0.5, AS0), (171.0, MH0 + 0.5, AS0), (171.0, MH0, AS0 - 0.0009), (171.0, MH0, AS0 + 0.0009)]


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
        json.dump({'source': 'SMDR v1.3 SMDR_Fit_Inputs (Martin-Robertson 1907.02500); other inputs = SMDR PDG-2025 '
                             'reference values (GF, MZ, alpha, light masses, Delta alpha_had)',
                   'pdg2026_point': [MH0, AS0], 'fits': res}, f, indent=1)
    print(json.dumps(res, indent=1))
