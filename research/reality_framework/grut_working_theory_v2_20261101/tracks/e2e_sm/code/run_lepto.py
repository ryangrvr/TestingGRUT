#!/usr/bin/env python3
"""e2e_sm driver: leptogenesis cases -> results/lepto.json (one core, ~20 s)."""
import json, os, sys
import sm_lepto as L

HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(HERE, '..', '..', 'e2e_kinetic', 'results')
RES = os.path.join(HERE, '..', 'results')
cases = [
    ('S1_HFa_1.3e11', 'fluid_HFa.json', 'seesaw_H1.3e11.json', 'S1_mt050', 1.3e11, {}),
    ('S1_HFb_1.3e11', 'fluid_HFb.json', 'seesaw_H1.3e11.json', 'S1_mt050', 1.3e11, {}),
    ('S1_HFa_washout_x0.5', 'fluid_HFa.json', 'seesaw_H1.3e11.json', 'S1_mt050', 1.3e11, dict(washout_scale=0.5)),
    ('S1_HFa_washout_x2', 'fluid_HFa.json', 'seesaw_H1.3e11.json', 'S1_mt050', 1.3e11, dict(washout_scale=2.0)),
    ('S1_HFa_no_nonthermal', 'fluid_HFa.json', 'seesaw_H1.3e11.json', 'S1_mt050', 1.3e11, dict(nt_scale=0.0)),
    ('S1_HFa_rescaled_8.4e10', 'fluid_HFa.json', 'seesaw_H1.3e11.json', 'S1_mt050', 8.4e10, {}),
    ('S2_HFa_1.3e11_history_inconsistent', 'fluid_HFa.json', 'seesaw_H1.3e11.json', 'S2_mt030', 1.3e11, {}),
]
out = dict(label='C on [Q]/[E] rates; model choices; see sm_lepto.py docstring', cases={})
for name, h, s, sc, H, kw in cases:
    r = L.run(os.path.join(KIN, h), os.path.join(RES, s), sc, H, **kw)
    out['cases'][name] = r
    print('%-36s etaB/delta %.3e  thermal-only %.3e  delta_needed %.2e  surv(3000) %.1e surv(1e4) %.2e' % (
        name, r['etaB_per_unit_delta'], r['etaB_thermal_only_per_unit_delta'], r['delta_needed'],
        r['washout_survival_from_t'].get('3000', float('nan')), r['washout_survival_from_t'].get('10000', float('nan'))))
json.dump(out, open(os.path.join(RES, 'lepto.json'), 'w'), indent=1)
