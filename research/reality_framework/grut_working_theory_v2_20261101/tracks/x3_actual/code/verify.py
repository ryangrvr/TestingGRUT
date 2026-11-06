#!/usr/bin/env python3
"""Audit saved evidence. Green execution checks do not mean a viable model."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from moments import TRACK, TRACKS, RES, dump


def main():
    checks=[]
    def check(name,ok,detail=None):
        checks.append(dict(name=name,passed=bool(ok),detail=detail))
    moments=json.loads((RES/'moment_feasibility.json').read_text())
    quantum=json.loads((RES/'quantum_controls.json').read_text())
    pops=json.loads((RES/'population_invariants.json').read_text())
    runs={p.stem:json.loads(p.read_text()) for p in (RES/'kinetic').glob('*.json')}
    check('all_seven_runs_present',len(runs)==7)
    base=runs['actual_ref_6400'];proxy=runs['proxy_ref_6400']
    ratios={k:v['summary']['K']/base['summary']['K']-1 for k,v in runs.items()}
    numerical=['actual_ref_3200','actual_ref_12800','actual_ref_dt_half']
    check('K_numerical_controls_within_registered_five_percent',max(abs(ratios[k]) for k in numerical)<.05,ratios)
    for name,r in runs.items():
        prefix=RES/'kinetic'/name
        input_bytes=gzip.decompress(Path(str(prefix)+'.in.gz').read_bytes())
        check(name+'_input_hash',hashlib.sha256(input_bytes).hexdigest()==r['input_sha256'])
        rec=np.loadtxt(Path(str(prefix)+'_rec.txt.gz'))
        check(name+'_records_finite',bool(np.isfinite(rec).all()))
        K=rec[:,1]**2*rec[:,2]
        check(name+'_K_from_raw',math.isclose(K[-1],r['summary']['K'],rel_tol=1e-13))
        err=3*r['state']['mP']**2*rec[0,2]**2/r['splice_ledger']['represented_energy']-1
        check(name+'_active_energy_continuity',abs(err)<1e-6,err)
        check(name+'_late_drift',r['checks']['last_decade_K_relative_drift']<1e-3)
        check(name+'_full_energy_failure_disclosed',not r['checks']['full_source_energy_continuity_pass'])
        check(name+'_no_point_promoted',r['accepted_cosmological_point'] is False)
    last=moments['rows'][-1]['species']
    check('robust_N1_moment_failure',last[0]['caps']['1.0']['E_over_E_min']<.5)
    check('independent_Fermi_minimization',moments['checks']['relative_error']<2e-6)
    check('all_half_occupation_native_moments_infeasible',all(not s['caps']['0.5']['feasible'] for s in last))
    check('line_unitarity',max(x['max_norm_error'] for x in quantum['line'])<1e-9)
    check('3d_mode_unitarity',max(x['norm_error'] for x in quantum['three_dimensional_mode'])<1e-9)
    full=[r for r in quantum['three_dimensional_mode'] if r['profile']=='full']
    check('3d_spatial_convergence',abs(full[3]['positive_probability']/full[4]['positive_probability']-1)<1e-8)
    check('3d_time_convergence',abs(full[1]['positive_probability']/full[2]['positive_probability']-1)<1e-3)
    check('3d_is_not_sea_certification',all(x['filled_sea'] is False for x in quantum['three_dimensional_mode']))
    hom=[r for r in quantum['three_dimensional_mode'] if r['profile']=='hom' and not r['static']]
    check('3d_independent_ODE',hom[-1]['independent_absolute_error']<3e-6,hom[-1]['independent_absolute_error'])
    vac=[r for r in quantum['three_dimensional_mode'] if r['static']]
    check('3d_static_vacuum_error_decreases',vac[-1]['positive_probability']<vac[0]['positive_probability']/10)
    check('reaction_charge_balance',len(pops['reactions'])==23 and all(r['charge_change']==0 for r in pops['reactions']))
    check('reaction_energy_balance',all(r['energy_error']==0 for r in pops['exact_energy_examples']))
    check('no_population_inferred',pops['cosmological_population_computed'] is False)
    change=base['summary']['K']/proxy['summary']['K']
    summary=dict(K_ratio_actual_to_proxy=change,K_change_percent=100*(change-1),
        conditional_fixed_charge_yield_ratio=change**(-1.5),
        numerical_control_relative_changes={k:ratios[k] for k in numerical},
        spectrum_control_relative_change=ratios['actual_dilute_6400'],
        species_partition_relative_change=ratios['actual_native_6400'],
        maximum_represented_energy_error=max(abs(x['checks']['represented_network_energy_error']) for x in runs.values()),
        source_scalar_budget_absmax=moments['source_scalar_budget_absmax'],
        scientific_acceptance=False,
        verdict='faithful native-moment splice fails; conditional network converges; cosmological closure unresolved',
        reasons=['native N1 energy and number violate the physical free-Fermi bound',
                 'full source energy splice exceeds 1e-6: omitted spectator and printed-column residual',
                 'H55 drain still transferred from H84; no full H55 quantum backreaction',
                 'no evolved core population or new accepted relic root'])
    dump(RES/'summary.json',summary)
    dump(RES/'verification.json',dict(passed=all(c['passed'] for c in checks),checks=checks,
         scientific_acceptance=False,scope='internal evidence and arithmetic checks, not outside scientific review'))
    print(json.dumps(summary,indent=2))
    print(sum(c['passed'] for c in checks),'/',len(checks),'evidence checks passed')
    for c in checks:
        if not c['passed']:print('FAIL',c)
    if not all(c['passed'] for c in checks):raise SystemExit(1)


if __name__=='__main__':main()
