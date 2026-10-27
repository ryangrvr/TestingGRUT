#!/usr/bin/env python3
"""cf_scale_skeptic: small closed-form numbers used in SKEPTIC.md (own code).
(i) P4-S: e(h) implied by e(M_P) = e* (one loop, b_e = 46/3); (ii) f_g and e* when hypercharge is calibrated at the
reduced Planck mass (own 2-loop g1 at 2.435e18 GeV, PDG 2026 direct m_t) instead of Buttazzo's M_Pl = 1.22e19 GeV.
Usage: python3 -B summary_numbers.py <out.json>"""
import json, sys
from pathlib import Path
import numpy as np
K = 16 * np.pi**2; B_E = 46 / 3; MP, H = 2.435e18, 5.33e12
sm = json.loads((Path(__file__).resolve().parents[1] / 'results' / 'sm_two_loop.json').read_text())
out = {}
for tag, g1 in (('Buttazzo_MPl_1.22e19', 0.6154), ('own_2loop_reduced_MP', sm['pdg2026']['mt_direct']['couplings_MPred_2loop']['g1'])):
    gY = g1 * np.sqrt(3 / 5); fg = gY**2 * (41 / 6) / K; es = np.sqrt(K * fg / B_E)
    eh = 1 / np.sqrt(1 / es**2 + 2 * B_E * np.log(MP / H) / K)
    out[tag] = dict(gY=gY, f_g=fg, e_star=es, A_crit=(432 + np.sqrt(80 * 7776)) * es**2 / K, e_at_h=eh)
    print(tag, {k: round(v, 5) for k, v in out[tag].items()})
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
