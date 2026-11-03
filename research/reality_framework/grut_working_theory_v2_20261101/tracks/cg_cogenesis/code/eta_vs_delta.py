#!/usr/bin/env python3
"""cg_cogenesis (Claude): eta_B at order-one lepton CP phase on the M2 histories, for candidates K1/K4.
Reuses e2e_card's corrected leptogenesis (excess-over-equilibrium source) and e2e_sm's sm_lepto.run UNCHANGED (read-only).
sm_lepto.run takes the seesaw row at delta_ref (eps_j and rates at that delta, so m-tilde shifts and washout at large
delta are included) and returns eta_B per unit delta; eta_B(delta) = that x delta_ref. [C on the model choices]"""
import json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
TR = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(TR, 'e2e_card', 'code'))
import lepto_fix as F  # noqa: E402  (read-only; its own __main__ is not run)
ETA_OBS = 6.12e-10
out = {}
for name, p, Hg in F.cases:
    if not os.path.exists(p) or name.startswith('final_GMg03'):
        continue
    hp, _ = F.patch(p, Hg, name)
    row = {}
    for d in (1e-3, 0.1, 0.3, 1.0):
        r = F.L.run(hp, F.SEE, 'S1_mt050', Hg, delta_ref=d)
        eta = r['etaB_per_unit_delta'] * d
        row['%g' % d] = dict(etaB=eta, ratio_to_obs=eta / ETA_OBS)
    out[name] = row
    print(name, ' '.join('d=%s: eta %.2e (x%.0f)' % (k, v['etaB'], v['ratio_to_obs']) for k, v in row.items()), flush=True)
json.dump(out, open(os.path.join(HERE, '..', 'results', 'eta_vs_delta.json'), 'w'), indent=1)
