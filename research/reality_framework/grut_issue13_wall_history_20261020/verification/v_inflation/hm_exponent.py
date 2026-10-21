"""Hawking-Moss / stochastic exponent check (verifier v_inflation): Lambda_1 from the 1D solver vs
B_HM = (8 pi^2/3) DeltaV/H^4 and the overdamped Kramers form Lambda_1 = sqrt2 a/(3 pi) exp(-B), a = |w''(0)|."""
import json, math
import numpy as np
import fp1d_verify as F
rows = []
for H in (100.0, 120.0, 141.42, 170.0, 200.0, 230.0, 251.19, 298.5, 335.0, 398.1):
    hg = 0.002 if H > 160 else 0.001
    L1, br, info = F.lambda1(H, hgoal=hg)
    L1b = F.lambda1(H, hgoal=hg / 2)[0]
    B = F.B_HM(H)
    K = F.kramers_lambda1(H)
    rows.append(dict(H=H, B_HM=B, Lambda1=L1, Lambda1_half_h=L1b, Lambda1_Richardson=(4 * L1b - L1) / 3,
                     ln_Lambda1_plus_B=math.log((4 * L1b - L1) / 3) + B, ratio_to_Kramers=((4 * L1b - L1) / 3) / K,
                     B_noxi=F.B_HM(H, xi_on=False)))
    print(rows[-1], flush=True)
# effective exponent: slope of ln Lambda1 vs B between neighbouring points
sl = []
for a, b in zip(rows[:-1], rows[1:]):
    sl.append(dict(H_pair=[a['H'], b['H']], dlnL_dB=(math.log(a['Lambda1_Richardson']) - math.log(b['Lambda1_Richardson'])) / (a['B_HM'] - b['B_HM'])))
json.dump(dict(rows=rows, slopes=sl), open('results_hm.json', 'w'), indent=1)
for s in sl:
    print(s)
