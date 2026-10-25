#!/usr/bin/env python3
"""Rest-mass accounting (Claude review, track `conversion`), own arithmetic.

Quoted inputs [Q] (bb211d1 parents, FIXED_YIELD_RESULTS.json joint_screen_population
and PARAMETER_POINT.json): late X per object 12521.802854958796, object mass
17064.537454786467 GeV, late total mass per object 18221.721841256607 GeV,
m_X = 0.0924089637375362 GeV.  Omega from my own profile solve (profile.py).
Exact (evolving ground-state) release from my E(Q) - E(Q - 2N) (profile_omega.json).
"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / 'results'
p = json.loads((RES / 'profile_omega.json').read_text())
OM = p['Q200000']['Omega']
MX = 0.0924089637375362
NX = 12521.802854958796
MOBJ = 17064.537454786467
MTOT = 18221.721841256607
pairs = NX / 2
anti_rest = pairs * MX
rel_frozen = 2 * OM * MX * pairs
rel_exact = p['exact_core_energy_release_GeV']
out = {}
for tag, rel in [('frozen_profile', rel_frozen), ('evolving_ground_state', rel_exact)]:
    rem = MTOT - rel
    out[tag] = dict(release_GeV=rel, release_fraction=rel / MTOT, remaining_ratio=rem / MTOT,
                    converted_X_fraction_of_remaining=anti_rest / rem,
                    all_free_X_fraction_of_remaining=2 * anti_rest / rem)
out['initial_anti_rest_fraction'] = anti_rest / MTOT
out['core_charge_fraction_transferred'] = 2 * pairs / 2e5
out['budget_check_GeV'] = dict(object=MOBJ, free_X=NX * MX, D=MTOT - MOBJ - NX * MX)
out['kinematics'] = dict(Omega=OM, E_out_cold=1 + 2 * OM, u_esc=((1 + 2 * OM) ** 2 - 1) ** 0.5,
                         K_out_warm=1.5 ** 0.5 - 1 + 2 * OM, p_out_warm=((1.5 ** 0.5 + 2 * OM) ** 2 - 1) ** 0.5,
                         E_out_cold_MeV=(1 + 2 * OM) * MX * 1e3,
                         Omega_after_complete_conversion=p['omega_over_m_after_complete_conversion'],
                         u_esc_after_complete_conversion=((1 + 2 * p['omega_over_m_after_complete_conversion']) ** 2 - 1) ** 0.5)
print(json.dumps(out, indent=1))
(RES / 'mass_fractions.json').write_text(json.dumps(out, indent=1) + '\n')
